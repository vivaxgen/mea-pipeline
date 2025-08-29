# this snakemake file will perform GWAS using gemma

include: "fastlmm.smk"
include: "blink.smk"

# GEMMA here runs with PLINK binary-PED (BED) format (ie. *.bed, *.bim, *.fam files)
# GEMMA automatically recognize -9 or NA as missing values in the PLINK binary-PED (BED) format
# (ie. in *.fam file)
# For covariates, GEMMA uses BIMBAM format and missing values should be coded as NA

rule prepare_covar_file:
    input:
        fn = '{pfx}/P@{pheno}/{vcf_notation}/covars.txt',
        pca = '{pfx}/P@{pheno}/{vcf_notation}/distance.eigenvec'
    output:
        fn = '{pfx}/P@{pheno}/{vcf_notation}/covars-gemma.txt',
    run:
        import pandas as pd

        covars = pd.read_table(input.fn, header=None)
        pca = pd.read_table(input.pca)

        # use covars to get the taxa order
        combined_df = pd.DataFrame({'taxa': covars[0], 'Intercept': [1] * len(covars)})
        combined_df = combined_df.merge(pca[['IID', 'PC1', 'PC2', 'PC3']], left_on='taxa', right_on='IID')
        combined_df.drop(['taxa', 'IID'], axis=1, inplace=True)
        covars_gemma = pd.concat([combined_df, covars.iloc[:, 2:]], axis=1)
        covars_gemma.to_csv(output.fn, sep='\t', header=False, index=False, na_rep='NA')


rule run_gemma:
    threads: 1
    input:
        test_fn = '{pfx}/P@{pheno}/{vcf_notation}/variants.bed',
        dist_fn = '{pfx}/P@{pheno}/{vcf_notation}/distance.bed',
        phenotype_fn = '{pfx}/P@{pheno}/{vcf_notation}/phenotype.txt',
        covars_fn = '{pfx}/P@{pheno}/{vcf_notation}/covars-gemma.txt'
    output:
        fn = '{pfx}/P@{pheno}/{vcf_notation}/gemma/linreg.assoc.txt',
    params:
        covars = (lambda w, input: f'-c {input.covars_fn}') if any(covars_columns) else '',
        similarity_file = (lambda w, input: f'--similarity-file {input.dist_fn}') if use_genetic_relation_matrix else '',
        bedprefix = lambda w, input: str(input.test_fn).removesuffix('.bed'),
        outprefix = lambda w, output: str(output.fn).rsplit('/', 1)[0],
    shell:
        "gemma -bfile {params.bedprefix} -lm 4"
        "  {params.covars}"
        "  -outdir {params.outprefix} -o linreg"


rule run_gemma_bslmm:
    threads: 1
    input:
        test_fn = '{pfx}/P@{pheno}/{vcf_notation}/variants.bed',
        dist_fn = '{pfx}/P@{pheno}/{vcf_notation}/distance.bed',
        phenotype_fn = '{pfx}/P@{pheno}/{vcf_notation}/phenotype.txt',
        covars_fn = '{pfx}/P@{pheno}/{vcf_notation}/covars.txt'
    output:
        fn = '{pfx}/P@{pheno}/{vcf_notation}/gemma/bslmm.param.txt',
    params:
        covars = (lambda w, input: f'-c {input.covars_fn}') if any(covars_columns) else '',
        similarity_file = (lambda w, input: f'--similarity-file {input.dist_fn}') if use_genetic_relation_matrix else '',
        bedprefix = lambda w, input: str(input.test_fn).removesuffix('.bed'),
        outprefix = lambda w, output: str(output.fn).rsplit('/', 1)[0],
    shell:
        "gemma -bfile {params.bedprefix} -bslmm 1"
        "  -outdir {params.outprefix} -o bslmm"


rule plot_gemma_result:
    threads: 1
    input:
        fn = '{pfx}/P@{pheno}/{vcf_notation}/gemma/linreg.assoc.txt',
    output:
        mht = '{pfx}/P@{pheno}/{vcf_notation}/gemma/linreg-mht.png',
        qq = '{pfx}/P@{pheno}/{vcf_notation}/gemma/linreg-qq.png',
    params:
        line1 = f'Dataset: {current_dataset}',
        line2 = lambda w: f'Phenotype: {w.pheno}',
        line3 = f'Covars: {" ".join(covars_columns)}',
    shell:
        '{cli} tab-plot-gwas-results --outqq {output.qq}'
        '  --outmht {output.mht}'
        '  --add-title-line "{params.line1}"'
        '  --add-title-line "{params.line2}"'
        '  --add-title-line "{params.line3}"'
        '  --columns chr,ps,p_wald'
        '  --use-bonferroni 0.05'
        '  {input.fn}'


rule annotate_gemma_bslmm_result:
    localrule: True
    input:
        fn = '{pfx}/P@{pheno}/{vcf_notation}/gemma/bslmm.param.txt',
    output:
        fn = '{pfx}/P@{pheno}/{vcf_notation}/gemma/bslmm.param.anno.txt',
    run:
        import pandas as pd

        df = pd.read_table(input.fn)
        df['eff'] = (df.beta * df.gamma).abs()
        df.rename(columns={'chr': 'CHROM_NO', 'ps': 'POS'}, inplace=True)
        df['CHROM'] = df.CHROM_NO
        df['SIGNAL'] = 0
        df.loc[df.gamma > 0.5, 'SIGNAL'] = 1
        df['size'] = (df.eff.abs() * 100).clip(lower=0.1)
        df.to_csv(output.fn, sep='\t', index=False)


rule plot_gemma_bslmm_result:
    threads: 1
    input:
        tsv = '{pfx}/P@{pheno}/{vcf_notation}/gemma/bslmm.param.anno.txt',
    output:
        plot = '{pfx}/P@{pheno}/{vcf_notation}/gemma/bslmm-mht.png',
    params:
        y_label = lambda w: f'PIP {w.pheno}',
        bedfile = f'--bedfile {fn}' if (fn := config['bedfile_highlight_manhattan']) else '',
    shell:
        "{cli} tab-plot-manhattan  -o {output.plot}  --column gamma"
        "  --hline  --size-column size --y-label '{params.y_label}'  {params.bedfile}"
        "  {input.tsv}"


# EOF
