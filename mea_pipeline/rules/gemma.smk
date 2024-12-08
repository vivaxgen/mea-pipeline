# this snakemake file will perform GWAS using gemma

include: "fastlmm.smk"


rule run_gemma:
    threads: 1
    input:
        test_fn = '{pfx}/P@{pheno}/{vcf_notation}/variants.bed',
        dist_fn = '{pfx}/P@{pheno}/{vcf_notation}/distance.bed',
        phenotype_fn = '{pfx}/P@{pheno}/{vcf_notation}/phenotype.txt',
        covars_fn = '{pfx}/P@{pheno}/{vcf_notation}/covars.txt'
    output:
        fn = '{pfx}/P@{pheno}/{vcf_notation}/gemma/linreg.assoc.txt',
    params:
        covars = (lambda w, input: f'--covars-file {input.covars_fn}') if any(covars_columns) else '',
        similarity_file = (lambda w, input: f'--similarity-file {input.dist_fn}') if use_genetic_relation_matrix else '',
        bedprefix = lambda w, input: str(input.test_fn).removesuffix('.bed'),
        outprefix = lambda w, output: str(output.fn).rsplit('/', 1)[0],
    shell:
        "gemma -bfile {params.bedprefix} -lm 4"
        "  -outdir {params.outprefix} -o linreg"


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
        '  {input.fn}'

# EOF
