# this snakemake file will perform GWAS using gemma

include: "fastlmm.smk"

rule calculate_blink_pca:
    threads: 1
    input:
        bed_fn = '{pfx}/P@{pheno}/{vcf_notation}/distance.bed',
        bim_fn = '{pfx}/P@{pheno}/{vcf_notation}/distance.bim',
        fam_fn = '{pfx}/P@{pheno}/{vcf_notation}/distance.fam',
    output:
        pca_fn = '{pfx}/P@{pheno}/{vcf_notation}/distance.eigenvec',
    params:
        prefix = lambda w, output: str(output.pca_fn).removesuffix('.eigenvec'),
    shell:
        "blink_linux --pca --file {params.prefix} --plink --out {params.prefix}"


rule prepare_phenotype:
    localrule: True
    input:
        phenotype_fn = '{pfx}/P@{pheno}/{vcf_notation}/phenotype.txt',
    output:
        phenotype_fn = '{pfx}/P@{pheno}/{vcf_notation}/variants.txt',
    run:
        import sys
        import pandas as pd

        phenotype_df = pd.read_table(input.phenotype_fn, header=None)
        phenotype_df = phenotype_df.drop(1, axis=1).rename(columns={0: 'taxa', 2: wildcards.pheno})

        phenotype_df.to_csv(output.phenotype_fn, sep='\t', index=False)


rule prepare_covars:
    localrule: True
    input:
        covars_fn = '{pfx}/P@{pheno}/{vcf_notation}/covars.txt',
        pca_fn = '{pfx}/P@{pheno}/{vcf_notation}/distance.eigenvec'
    output:
        covars_fn = '{pfx}/P@{pheno}/{vcf_notation}/variants.cov',
    run:
        import sys
        import pandas as pd

        covars_df = pd.read_table(input.covars_fn, header=None)
        covars_df = covars_df.drop(covars_df.columns[1], axis=1)
        pca_df = pd.read_table(input.pca_fn)

        final_df = pca_df[['FID', 'PC1', 'PC2', 'PC3']].rename(columns={'FID': 'taxa'})
        final_df = final_df.merge(covars_df, left_on='taxa', right_on=0).drop(0, axis=1)

        final_df.to_csv(output.covars_fn, sep='\t', index=False, header=False)


rule run_blink:
    threads: 1
    input:
        bed_fn = '{pfx}/P@{pheno}/{vcf_notation}/variants.bed',
        bim_fn = '{pfx}/P@{pheno}/{vcf_notation}/variants.bim',
        fam_fn = '{pfx}/P@{pheno}/{vcf_notation}/variants.fam',
        phenotype_fn = '{pfx}/P@{pheno}/{vcf_notation}/variants.txt',
        covars_fn = '{pfx}/P@{pheno}/{vcf_notation}/variants.cov',
    output:
        assoc_fn = '{pfx}/P@{pheno}/{vcf_notation}/variants_{pheno}_GWAS_result.txt',
    params:
        prefix = lambda w, input: input.phenotype_fn.removesuffix('.txt'),
    shell:
        "blink_linux --gwas --file {params.prefix}  --out {params.prefix} --plink"


rule annotate_result:
    input:
        assoc_fn = '{pfx}/P@{pheno}/{vcf_notation}/variants_{pheno}_GWAS_result.txt'
    output:
        table_fn = '{pfx}/P@{pheno}/{vcf_notation}/gwas.tsv',
    shell:
        "mea-pl tab-annotate-signals -o {output.table_fn}"
        "  --column p_value --neg-log10"
        "  --quantile 0.999 --minvariants 1"
        "  {input.assoc_fn}:chr,pos"


rule plot_manhattan:
    input:
        table_fn = '{pfx}/P@{pheno}/{vcf_notation}/gwas.tsv'
    output:
        plot_fn = '{pfx}/P@{pheno}/{vcf_notation}/mht.png',
    shell:
        "mea-pl tab-plot-manhattan -o {output.plot_fn} --column=\"-log10(p_value)\" {input.table_fn}"


# EOF