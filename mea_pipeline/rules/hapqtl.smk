
include: "fastlmm.smk"


rule prepare_hapQTL_input:
    # this rule prepares bimbam files for hapQTL inputs
    threads: 1
    input:
        vcf_fn = vcf_file,
        covars_fn = covars_file,
    output:
        geno_fn = expand(
            '{{pfx}}/P@{{pheno}}/{{test_notation}}/hapQTL/{region}.geno.txt',
            region=REGIONS
        ),
        pos_fn = expand(
            '{{pfx}}/P@{{pheno}}/{{test_notation}}/hapQTL/{region}.pos.txt',
            region=REGIONS
        ),
        pheno_fn = '{pfx}/P@{pheno}/{test_notation}/hapQTL/test.pheno.txt',
    params:
        outdir = lambda w, output: str(output.pheno_fn).removesuffix('/test.pheno.txt'),
        regions = ' '.join(f"-r {region}" for region in REGIONS),
        covars_fn = lambda w, input: ("--covariant-file " + input.covars_fn + ":WGSID," + ",".join(covars_columns)) if any(covars_columns) else "",
        translation_file = chromtranslation_file,
    shell:
        "mea-pl vcf-convert-to-bimbam"
        "  --phenotype-file {input.covars_fn}:WGSID,{wildcards.pheno}"
        "  {params.covars_fn}"
        "  {params.regions}"
        "  --translation-file {params.translation_file}"
        #"  --outprefix {wildcards.region}"
        "  --outdir {params.outdir}"
        "  {input.vcf_fn}"
        #"  ../lmm/P@CQ_RELAX_CODE/test/atomized.vcf.gz"


rule prepare_hapQTL_input_complete_region:
    input:
        vcf_fn = vcf_file,
        covars_fn = covars_file,
    output:
        geno_fn = f'{{pfx}}/P@{{pheno}}/{{test_notation}}/hapQTL/{complete_region}.geno.txt',
        pos_fn = f'{{pfx}}/P@{{pheno}}/{{test_notation}}/hapQTL/{complete_region}.pos.txt',
        pheno_fn = '{pfx}/P@{pheno}/{test_notation}/hapQTL/test.pheno.txt',
    params:
        outdir = lambda w, output: str(output.pheno_fn).removesuffix('/test.pheno.txt'),
        regions = "",
        covars_fn = lambda w, input: ("--covariant-file " + input.covars_fn + ":WGSID," + ",".join(covars_columns)) if any(covars_columns) else "",
        translation_file = chromtranslation_file,
    shell:
        "mea-pl vcf-convert-to-bimbam"
        "  --phenotype-file {input.covars_fn}:WGSID,{wildcards.pheno}"
        "  {params.covars_fn}"
        "  {params.regions}"
        "  --translation-file {params.translation_file}"
        "  --outprefix {complete_region}"
        "  --outdir {params.outdir}"
        "  {input.vcf_fn}"
        #"  ../lmm/P@CQ_RELAX_CODE/test/atomized.vcf.gz"

ruleorder: prepare_hapQTL_input_complete_region > prepare_hapQTL_input

rule prepare_random_seed:
    threads: 1
    output:
        seed_fn = '{pfx}/P@{pheno}/{test_notation}/hapQTL/seed.txt',
    run:
        import os
        
        random_data = os.urandom(3)
        seed = int.from_bytes(random_data, byteorder="big")
        with open(output.seed_fn, 'w') as f:
            f.write(str(seed))


def get_seed(wildcards, input):
    with open(input.seed_fn) as f:
        return f.read().strip()


rule run_hapQTL:
    threads: 1
    input:
        geno_fn = '{pfx}/P@{pheno}/{test_notation}/hapQTL/{region}.geno.txt',
        pos_fn = '{pfx}/P@{pheno}/{test_notation}/hapQTL/{region}.pos.txt',
        pheno_fn = '{pfx}/P@{pheno}/{test_notation}/hapQTL/test.pheno.txt',
        seed_fn = '{pfx}/P@{pheno}/{test_notation}/hapQTL/seed.txt',
    output:
        fn = '{pfx}/P@{pheno}/{test_notation}/hapQTL/output/{region}.bf.txt',
        snpinfo = '{pfx}/P@{pheno}/{test_notation}/hapQTL/output/{region}.snpinfo.txt',
    params:
        upper_cluster = config.get('hapqtl_upper_cluster', 3),
        lower_cluster = config.get('hapqtl_lower_cluster', 8),
        runs = config.get('hapqtl_runs', 1),
        linear_step = config['hapqtl_linear_step'],   # might need to bump to 80
        quadratic_step = config['hapqtl_quadratic_step'],  # default is 30
        morgan = "-morgan {m}" if (m := config['hapqtl_morgan']) else "",
        seed = get_seed,
        #outprefix = lambda w, output: str(output.fn).removesuffix('.bf.txt'),
    shell:
        "cd {wildcards.pfx}/P@{wildcards.pheno}/{wildcards.test_notation}/hapQTL"
        " && "
        "hapQTL-lin"
        "  -g {wildcards.region}.geno.txt"
        "  -p 1 -pos {wildcards.region}.pos.txt"
        "  -FILE test.pheno.txt"
        "  -C {params.upper_cluster} -c {params.lower_cluster}"
        "  -e {params.runs} -w {params.linear_step} -s {params.quadratic_step}"
        "  --pk1"
        "  -R {params.seed}"
        "  -o {wildcards.region}"
        " || true"

        #"hapQTL-lin -g out.geno.txt -p 1 -pos out.pos.txt -FILE out.pheno.txt -C 3 -c 8 -o p-C3-c8-e3-s50 -e 3 -s 50"


rule gather_hapQTL_bf:
    localrule: True
    input:
        bf_fn = expand(
            '{{pfx}}/P@{{pheno}}/{{test_notation}}/hapQTL/output/{region}.bf.txt',
            region=REGIONS
        ),
        pos_fn = expand(
            '{{pfx}}/P@{{pheno}}/{{test_notation}}/hapQTL/output/{region}.snpinfo.txt',
            region=REGIONS
        )
    output:
        bf_fn = '{pfx}/P@{pheno}/{test_notation}/hapQTL/all.bf.txt',
    params:
        use_id = False,
    run:

        import pandas as pd

        dfs = []

        for bf_fn, pos_fn in zip(input.bf_fn, input.pos_fn):
            bf_df = pd.read_csv(bf_fn, sep='\t')
            pos_df = pd.read_csv(pos_fn, sep='\t')

            df = pd.concat([pos_df, bf_df], axis=1)

            dfs.append(df)

        bf_df = pd.concat(dfs)

        rename_columns = {}
        for col in bf_df.columns:
            rename_columns[col] = col.strip()
        bf_df.rename(columns=rename_columns, inplace=True)
        bf_df.rename(columns={'chr': 'CHROM_NO', 'pos': 'POS'}, inplace=True)
        # from column rs, split the value by ':' and take the 1st part
        # import IPython; IPython.embed()
        if params.use_id:
            bf_df[['CHROM', 'POS']] = bf_df['rs'].str.split(':', n=1, expand=True)
        else:
            bf_df['CHROM'] = bf_df['rs'].str.split(':').str[0]
        bf_df['SIGNAL'] = 0
        bf_df.loc[(bf_df.bf2 >= 2) | (bf_df.bf1 >= 2), 'SIGNAL'] = 1
        bf_df['UPPER_THRESHOLD'] = 2

        bf_df.to_csv(output.bf_fn, sep='\t', index=False)


use rule gather_hapQTL_bf as gather_hapQTL_bf_complete_region with:
    localrule: True
    input:
        bf_fn = [f'{{pfx}}/P@{{pheno}}/{{test_notation}}/hapQTL/output/{complete_region}.bf.txt',],
        pos_fn = [f'{{pfx}}/P@{{pheno}}/{{test_notation}}/hapQTL/output/{complete_region}.snpinfo.txt',],
    output:
        bf_fn = f'{{pfx}}/P@{{pheno}}/{{test_notation}}/hapQTL/{complete_region}.bf.txt',
    params:
        use_id = True,


rule plot_hapQTL_localhap:
    input:
        bf_fn = '{pfx}/P@{pheno}/{test_notation}/hapQTL/all.bf.txt',
    output:
        plot_fn = '{pfx}/P@{pheno}/{test_notation}/hapQTL/all.bf-hap.png',
    params:
        y_label = lambda w: f'log (BF_Haplotype) {w.pheno}',
        bedfile = f'--bedfile {fn}' if (fn := config['bedfile_highlight_manhattan']) else '',
    shell:
        "{cli} tab-plot-manhattan  -o {output.plot_fn}  --column bf2"
        "  --hline --size-column pv2 --size-factor 10"
        "  --y-label '{params.y_label}'  {params.bedfile}"
        "  {input.bf_fn}"


use rule plot_hapQTL_localhap as plot_hapQTL_localhap_complete_region with:
    input:
        bf_fn = f'{{pfx}}/P@{{pheno}}/{{test_notation}}/hapQTL/{complete_region}.bf.txt',
    output:
        plot_fn = f'{{pfx}}/P@{{pheno}}/{{test_notation}}/hapQTL/{complete_region}.bf-hap.png',


rule plot_hapQTL_pval_localhap:
    input:
        bf_fn = '{pfx}/P@{pheno}/{test_notation}/hapQTL/all.bf.txt',
    output:
        plot_fn = '{pfx}/P@{pheno}/{test_notation}/hapQTL/all.pval-hap.png',
    params:
        y_label = lambda w: f'-log(PVal_Haplotype) {w.pheno}',
        bedfile = f'--bedfile {fn}' if (fn := config['bedfile_highlight_manhattan']) else '',
    shell:
        "{cli} tab-plot-manhattan  -o {output.plot_fn}  --column pv2"
        "  --size-column bf2 --size-factor 10"
        "  --y-label '{params.y_label}'  {params.bedfile}"
        "  {input.bf_fn}"


rule plot_hapQTL_SNP:
    input:
        bf_fn = '{pfx}/P@{pheno}/{test_notation}/hapQTL/all.bf.txt',
    output:
        plot_fn = '{pfx}/P@{pheno}/{test_notation}/hapQTL/all.bf-snp.png',
    params:
        y_label = lambda w: f'log (BF_SNP) {w.pheno}',
        bedfile = f'--bedfile {fn}' if (fn := config['bedfile_highlight_manhattan']) else '',
    shell:
        "{cli} tab-plot-manhattan  -o {output.plot_fn}  --column bf1"
        "  --hline  --size-column pv1 --size-factor 10"
        "  --y-label '{params.y_label}'  {params.bedfile}"
        "  {input.bf_fn}"


use rule plot_hapQTL_SNP as plot_hapQTL_SNP_complete_region with:
    input:
        bf_fn = f'{{pfx}}/P@{{pheno}}/{{test_notation}}/hapQTL/{complete_region}.bf.txt',
    output:
        plot_fn = f'{{pfx}}/P@{{pheno}}/{{test_notation}}/hapQTL/{complete_region}.bf-snp.png',


rule plot_hapQTL_pval_SNP:
    input:
        bf_fn = '{pfx}/P@{pheno}/{test_notation}/hapQTL/all.bf.txt',
    output:
        plot_fn = '{pfx}/P@{pheno}/{test_notation}/hapQTL/all.pval-snp.png',
    params:
        y_label = lambda w: f'-log(PVal_SNP) {w.pheno}',
        bedfile = f'--bedfile {fn}' if (fn := config['bedfile_highlight_manhattan']) else '',
    shell:
        "{cli} tab-plot-manhattan  -o {output.plot_fn}  --column pv1"
        "  --size-column bf1 --size-factor 10"
        "  --y-label '{params.y_label}'  {params.bedfile}"
        "  {input.bf_fn}"


# EOF
