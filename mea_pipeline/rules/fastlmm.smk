import pathlib

include: "utils.smk"

current_dataset = pathlib.Path.cwd().parent.name + '/' + pathlib.Path.cwd().name

complete_region = config.get('complete_region')
chromtranslation_file = config.get('chromtranslation_file')
covars_file = config.get('covars_file')
#phenotype_columns = config.get('phenotype_columns')
covars_columns = config.get('covars_columns', [])
source_dir = config.get('source_dir')
out_dir = config.get('out_dir')


vcf_file = f'{source_dir}/{{test_notation}}/concat/{complete_region}.vcf.gz'
fws_file = f'{source_dir}/{{test_notation}}/concat/fws-{complete_region}/fws.tsv'
ibd_file = f'{source_dir}/{{test_notation}}/concat/ibd-{complete_region}/ibd.tsv'
dis_file = f'{source_dir}/dist/concat/{complete_region}.vcf.gz'

set_het_as_missing = config['set_het_as_missing']
use_genetic_relation_matrix = config['use_genetic_relation_matrix']



# TODO:
# - handle missing values in phenotype and covariate files, currently set to "" (empty string)
#   see prepare_phenotype_and_covars_file rule
# - use hmmIBD to calculate proportion of shared-IBD, convert it to
#   similarity matrix (eg. 1 - shared_proportion), and save it as
#   .npz file for consumption of K0 argument of single_snp()


rule prepare_samples:
    input:
        vcf_fn = vcf_file,
    output:
       samples = f'{out_dir}/P@{{pheno}}/{{test_notation}}/vcf-samples.txt',
    shell:
        'bcftools query -l {input.vcf_fn} > {output.samples}'


rule vcf_to_ped:
    threads: 4
    input:
        vcf_fn = vcf_file,
    output:
        ped_fn = '{pfx}/_geno/variants.ped',
        map_fn = '{pfx}/_geno/variants.map',
    params:
        outprefix = lambda w, output: str(output.ped_fn).removesuffix('.ped'),
        #het2miss = '--set-hh-missing' if set_het_as_missing else '',
    shell:
        'plink2 --out {params.outprefix} --recode ped --allow-extra-chr '
        ' --threads {threads} '
        ' --double-id --vcf {input.vcf_fn} && '
        'spcli $PYS/gen/translate-chrom-name.py '
        ' --translation-file {chromtranslation_file} '
        ' -o {output.map_fn} {output.map_fn}'


use rule vcf_to_ped as vcf_to_ped2 with:
    input:
        vcf_fn = dis_file
    output:
        ped_fn = '{pfx}/_geno/distance.ped',
        map_fn = '{pfx}/_geno/distance.map',


rule vcf_to_bed:
    threads: 4
    input:
        vcf_fn = vcf_file,
        sample_fn = '{pfx}/P@{pheno}/{test_notation}/samples.txt',
        phenotype_fn = '{pfx}/P@{pheno}/{test_notation}/phenotype.txt',
    output:
        bed_fn = '{pfx}/P@{pheno}/{test_notation}/variants.bed',
        bim_fn = '{pfx}/P@{pheno}/{test_notation}/variants.bim',
        fam_fn = '{pfx}/P@{pheno}/{test_notation}/variants.fam',
    params:
        outprefix = lambda w, output: str(output.bed_fn).removesuffix('.bed'),
        het2miss = '--set-hh-missing' if set_het_as_missing else '',
    shell:
        'plink2 --out {params.outprefix} --make-bed --allow-extra-chr'
        ' --threads {threads}'
        ' {params.het2miss}'
        #'  --geno 0.5'
        ' --keep {input.sample_fn}'
        ' --double-id --vcf {input.vcf_fn}'
        ' && '
        'mea-pl translate-chrom-name'
        ' --translation-file {chromtranslation_file}'
        ' -o {output.bim_fn} {output.bim_fn}'
        ' && '
        'mea-pl tab-modify-fam'
        ' --phenotype-file {input.phenotype_fn}'
        ' --outfile {output.fam_fn} {output.fam_fn}'


rule vcf_to_bed2:
    input:
        vcf_fn = dis_file,
        sample_fn = '{pfx}/P@{pheno}/{vcf_notation}/samples.txt',
        #fws_fn = fws_file,
    output:
        bed_fn = '{pfx}/P@{pheno}/{vcf_notation}/distance.bed',
        bim_fn = '{pfx}/P@{pheno}/{vcf_notation}/distance.bim',
        fam_fn = '{pfx}/P@{pheno}/{vcf_notation}/distance.fam',
    params:
        outprefix = lambda w, output: str(output.bed_fn).removesuffix('.bed'),
        het2miss = '--set-hh-missing' if set_het_as_missing else '',
    shell:
        'plink2 --out {params.outprefix} --make-bed --allow-extra-chr'
        ' --threads {threads}'
        ' {params.het2miss}'
        #'  --geno 0.5'
        ' --keep {input.sample_fn}'
        ' --double-id --vcf {input.vcf_fn}'
        ' && '
        'mea-pl translate-chrom-name'
        ' --translation-file {chromtranslation_file}'
        ' -o {output.bim_fn} {output.bim_fn}'


rule prepare_phenotype_and_covars_file:
    localrule: True
    input:
        covars_fn = covars_file,
        sample_fn = '{pfx}/P@{pheno}/{test_notation}/vcf-samples.txt',
    output:
        phenotype_fn = '{pfx}/P@{pheno}/{test_notation}/phenotype.txt',
        covars_fn = '{pfx}/P@{pheno}/{test_notation}/covars.txt',
        sample_fn = '{pfx}/P@{pheno}/{test_notation}/samples.txt',

    run:
        import pandas as pd
        import numpy as np

        # prepare covariate + phenotype
        if input.covars_fn.endswith('.feather'):
            covars_df = pd.read_feather(input.covars_fn)
        else:
            covars_df = pd.read_table(input.covars_fn)

        if 'FWS' in covars_columns:
            # get FWS f
            pass

        covars_df = covars_df.loc[:, ['WGSID', wildcards.pheno] + covars_columns]

        sample_df = pd.read_table(input.sample_fn, header=None)
        sample_df.rename(columns={0: 'WGSID'}, inplace=True)
        merged_df = sample_df.merge(covars_df)
        #import IPython; IPython.embed()
        #merged_df['log_CQ_EC50'] = np.log(merged_df.CQ_EC50)
        #merged_df.loc[merged_df.CQ_EC50 > 300, 'CQ_EC50'] = 300
        merged_df = merged_df.dropna()

        # write phenotype
        phenotype_df = merged_df.loc[:, ['WGSID', 'WGSID', wildcards.pheno]]
        phenotype_df.to_csv(output.phenotype_fn, header=False, index=False, sep='\t', na_rep='')

        # write covariates
        if any(covars_columns):
            print('Preparing covariates: ', covars_columns)
        covariates_df = merged_df.loc[:, ['WGSID', 'WGSID'] + covars_columns]
        covariates_df.to_csv(output.covars_fn, header=False, index=False, sep='\t', na_rep='')

        # write samples
        #with open(output.sample_fn, 'wt') as sample_out:
        #    sample_out.write('#IID\n')
        #    sample_out.write('\n'.join(merged_df.WGSID))
        merged_df['#FID'] = merged_df.WGSID
        merged_df['IID'] = merged_df.WGSID
        merged_df[['#FID', 'IID']].to_csv(output.sample_fn, index=False, sep='\t')


rule plot_phenotype:
    threads: 1
    input:
        fn = '{pfx}/{pheno}/phenotype.txt'
    output:
        fn = '{pfx}/{pheno}/phenotype.png'
    run:
        import pandas as pd
        from matplotlib import pyplot as plt
        import seaborn as sns

        phenotype_df = pd.read_table(input.fn, header=None)
        #import IPython; IPython.embed()
        sns.scatterplot(x=range(len(phenotype_df)), y=sorted(phenotype_df[2]))
        plt.savefig(output.fn)


rule run_fastlmm:
    threads: 12
    input:
        ped_fn = '{pfx}/P@{pheno}/{vcf_notation}/variants.bed',
        dist_fn = '{pfx}/P@{pheno}/{vcf_notation}/distance.bed',
        phenotype_fn = '{pfx}/P@{pheno}/{vcf_notation}/phenotype.txt',
        covars_fn = '{pfx}/P@{pheno}/{vcf_notation}/covars.txt'
    output:
        fn = '{pfx}/P@{pheno}/{vcf_notation}/normal/gwas_results.feather',
    params:
        covars = (lambda w, input: f'--covars-file {input.covars_fn}') if any(covars_columns) else '',
        similarity_file = (lambda w, input: f'--similarity-file {input.dist_fn}') if use_genetic_relation_matrix else '',
    shell:
        "{cli} gwas-run-fastlmm --phenotype-file {input.phenotype_fn} "
        " {params.covars} "
        " {params.similarity_file} "
        " -o {output.fn} {input.ped_fn}"


rule plot_results:
    localrule: True
    input:
        fn = '{pfx}/P@{pheno}/{vcf_notation}/{lmm_type}/gwas_results.feather',
    output:
        qq = '{pfx}/P@{pheno}/{vcf_notation}/{lmm_type}/qq.png',
        mht = '{pfx}/P@{pheno}/{vcf_notation}/{lmm_type}/mht.png',
    params:
        line1 = f'Dataset: {current_dataset}',
        line2 = lambda w: f'Phenotype: {w.pheno} {"[warped]" if w.lmm_type == "warped" else ""}',
        line3 = f'Covars: {" ".join(covars_columns)}',
    shell:
        '{cli} tab-plot-gwas-results --outqq {output.qq}'
        '  --outmht {output.mht}'
        '  --add-title-line "{params.line1}"'
        '  --add-title-line "{params.line2}"'
        '  --add-title-line "{params.line3}"'
        '  --use-bonferroni 0.05'
        '  {input.fn}'


rule run_warpedlmm:
    threads: 12
    input:
        snp_fn = '{pfx}/{pheno}/variants.bed',
        dist_fn = '{pfx}/{pheno}/distance.bed',
        phenotype_fn = '{pfx}/{pheno}/phenotype.txt',
        covars_fn = '{pfx}/{pheno}/covars.txt',
    output:
        fn = '{pfx}/{pheno}/warped/gwas_results.feather',
    params:
        covars = (lambda w, input: f'--covariates {input.covars_fn}') if any(covars_columns) else '',
        k0_snp_file = (lambda w, input: f'--k0-snp-file {input.dist_fn}') if use_genetic_relation_matrix else '',
    shell:
        'python3 -m warpedlmm {input.snp_fn} {input.phenotype_fn} '
        '  --output_directory {wildcards.pfx}/{wildcards.pheno}/warped --save '
        '  {params.covars} '
        '  {params.k0_snp_file}'


# EOF
