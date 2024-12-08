from mea_pipeline import cerr


include: "utils.smk"


refseq_file = config["refseq_file"]
sample_file = config["sample_file"]
variant_dir = config["variant_dir"]
variant_file = config["variant_file"]

snpEff_config_file = config["snpEff_config_file"]
snpEff_data_dir = config["snpEff_data_dir"]
snpEff_db = config["snpEff_db"]

gff_file = config["gff_file"]
POI_slope = config["POI_slope"]

query_list = config.get("queries", [])
target_list = config.get("targets", [])
outfile_list = config.get("outfiles", [])

complete_region = config["complete_region"]


strict_samples = config["strict_samples"]

# strict variants indicate whether we are using -T (strict only for the position)
# or using -R (check for overlaps as well in case of indels) during variant filtering
# with bcftools view
strict_variants = config["strict_variants"]

cerr(f"query_list: {query_list}")


def out_files(queries=None, outfiles=None):
    if queries is None:
        queries = query_list
    if outfiles is None:
        outfiles = outfile_list
    file_list = []
    for q in queries:
        if q.endswith('/'):
            q = q.removesuffix('/')
        file_list += [f"{q}/{fn}" for fn in outfiles]
    print("File list: ", file_list)
    return file_list


def vcf_files(queries=None):
    if queries is None:
        queries = query_list
    file_list = []
    for q in queries:
        file_list += [f"{q}/vcfs/{reg}.vcf.gz.tbi" for reg in REGIONS]
    return file_list


def qc_files(queries=None):
    if queries is None:
        queries = query_list
    file_list = []
    for q in queries:
        file_list += [f"{q}/QC/{fn}" for fn in ["samples.png", "variants.png"]]
    return file_list


def concatenated_vcf_files(queries=None):
    if queries is None:
        queries = query_list
    file_list = []
    for q in queries:
        file_list.append(f"{q}/concat/{complete_region}.vcf.gz.tbi")
    return file_list


def concatenated_qc_files(queries=None):
    if queries is None:
        queries = query_list
    file_list = []
    for q in queries:
        file_list += [f"{q}/concat/{complete_region}.{fn}-qc.png" for fn in ["samples", "variants"]]
    return file_list


def region_file(w):
    if variant_dir:
        return f"{variant_dir}/{w.reg}.bed.gz"
    elif variant_file:
        return variant_file
    return "NO-VARIANTS-DIR-NOR-FILE"


rule files:
    input:
        out_files(),


rule VCF:
    input:
        out_files(outfiles=[f'vcfs/{reg}.vcf.gz.tbi' for reg in REGIONS]),
        #vcf_files(),


rule QC:
    input:
        qc_files(),


rule concatdVCF:
    input:
        concatenated_vcf_files(),


rule concatdQC:
    input:
        concatenated_qc_files(),


# variant manipulation

rule set_hets_depth:
    threads: 2
    input:
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
        idx="{pfx}/vcfs/{reg}.vcf.gz.tbi",
    output:
        vcf="{pfx}-HETd{dhets}r{rhets}-d{mindepth}/vcfs/{reg}.vcf.gz",
    log:
        log1="{pfx}-HETd{dhets}r{rhets}-d{mindepth}/logs/mea_vcf-set-GT-{reg}.log",
        log2="{pfx}-HETd{dhets}r{rhets}-d{mindepth}/logs/bcftools-fill-tags-{reg}.log",
    shell:
        "{cli} vcf-set-GT --min-minor-depth {wildcards.dhets} "
        " --min-minor-ratio {wildcards.rhets}  --minimum-depth {wildcards.mindepth}"
        " {input.vcf} 2> {log.log1}"
        " | bcftools +fill-tags -o {output} -- -t TYPE,F_MISSING,AC,AF,AN"


rule set_depth:
    threads: 2
    input:
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
        idx="{pfx}/vcfs/{reg}.vcf.gz.tbi",
    output:
        vcf="{pfx}-d{mindepth}/vcfs/{reg}.vcf.gz",
    log:
        log1="{pfx}-d{mindepth}/logs/bcftools-setGT-{reg}.log",
    shell:
        'bcftools +setGT {input.vcf} -- -t q -n . -i "FORMAT/DP<{wildcards.mindepth}" 2> {log.log1}'
        " | bcftools view -o {output}"


ruleorder: set_hets_depth > set_depth

HETOPTS = {
    'REF': '--set-het-to-ref',
    'ALT': '--set-het-to-alt',
    'MISS': '--set-het-to-missing',
    'MAJ': '--min-minor-depth 2 --min-minor-ratio 0.49 --set-het-to-missing',
}

rule convert_hets:
    threads: 1
    input:
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
        idx="{pfx}/vcfs/{reg}.vcf.gz.tbi",
    output:
        vcf="{pfx}-HET2{allele}/vcfs/{reg}.vcf.gz",
    params:
        flag=lambda w: HETOPTS[w.allele]
    shell:
        "{cli} vcf-set-GT -o {output.vcf} {params.flag} {input.vcf}"


# variant annotation


rule snpEff_annotate:
    threads: 2
    input:
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
        idx="{pfx}/vcfs/{reg}.vcf.gz.tbi",
        config_file=snpEff_config_file,
        data_dir=snpEff_data_dir,
    output:
        vcf="{pfx}-ANN/vcfs/{reg}.vcf.gz",
    log:
        log1="{pfx}-ANN/logs/snpEff-ann-{reg}.log"
    shell:
        "snpEff ann -c {input.config_file} -dataDir {input.data_dir} {snpEff_db} {input.vcf} 2> {log.log1}"
        " | bcftools view -o {output.vcf}"


rule bcsq_annotate:
    threads: 2
    input:
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
        idx="{pfx}/vcfs/{reg}.vcf.gz.tbi",
        refseq_file=refseq_file,
        gff_file=gff_file,
    output:
        vcf="{pfx}-BCSQ/vcfs/{reg}.vcf.gz",
    params:
        flags="--local-csq"
    log:
        log1="{pfx}-BCSQ/logs/bcftools-csq-{reg}.log"
    shell:
        "bcftools csq -f {input.refseq_file} -g {input.gff_file} {params.flags} "
        " -o {output.vcf} {input.vcf} 2> {log.log1}"


# variant and sample filtering


rule flt_pass:
    threads: 2
    input:
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
        idx="{pfx}/vcfs/{reg}.vcf.gz.tbi",
    output:
        vcf="{pfx}-PASS/vcfs/{reg}.vcf.gz",
    log:
        log1="{pfx}-PASS/logs/bcftool-f-PASS-{reg}.log",
    shell:
        "bcftools view -f PASS {input.vcf} 2> {log.log1}"
        " | bcftools view -o {output}"


rule flt_variants_samples:
    threads: 3
    input:
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
        vcf_idx="{pfx}/vcfs/{reg}.vcf.gz.csi",
        variant_file=variant_file,
        sample_file=sample_file,
    output:
        vcf="{pfx}-variants-samples/vcfs/{reg}.vcf.gz",
    log:
        log1="{pfx}-variants-samples/logs/bcftools-view-R-{reg}.log",
        log2="{pfx}-variants-samples/logs/bcftools-view-S-{reg}.log",
    shell:
        "bcftools view -R {input.variant_file} {input.vcf} 2> {log.log1}"
        " | bcftools view -S {input.sample_file} -o {output} 2> {log.log2}"


rule flt_samples:
    # filter VCFs by samples
    # requires:
    #   source_dir
    #   samples_file [a text file containing sample ids]
    threads: 2
    input:
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
        sample_file=sample_file,
    output:
        vcf="{pfx}-samples/vcfs/{reg}.vcf.gz",
    log:
        log1="{pfx}-samples/logs/bcftools-view-S-{reg}.log",
    params:
        strict_flag="" if strict_samples else "--force-samples",
    shell:
        "bcftools view -S {input.sample_file} {params.strict_flag}"
        " -o {output.vcf} {input.vcf} 2> {log.log1}"


use rule flt_samples as flt_samples_fn with:
    input:
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
        sample_file="{fn}",
    output:
        vcf="{pfx}-S{{fn}}/vcfs/{reg}.vcf.gz",
    log:
        log1="{pfx}-S{{fn}}/logs/bcftools-view-S-{reg}.log",


rule flt_variants:
    # filter VCFs by variants
    threads: 2
    input:
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
        vcf_idx="{pfx}/vcfs/{reg}.vcf.gz.tbi",
        region_file=region_file,
    output:
        vcf="{pfx}-variants/vcfs/{reg}.vcf.gz",
    params:
        region_flag="-T" if strict_variants else "-R",
    log:
        log1="{pfx}-variants/logs/bcftools-view-R-{reg}.log",
    shell:
        "bcftools view {params.region_flag} {input.region_file} -o {output.vcf} {input.vcf} 2> {log.log1}"


use rule flt_variants as flt_variants_fn with:
    input:
    input:
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
        vcf_idx="{pfx}/vcfs/{reg}.vcf.gz.tbi",
        region_file="{fn}",
    output:
        vcf="{pfx}-V{{fn}}/vcfs/{reg}.vcf.gz",
    log:
        log1="{pfx}-V{{fn}}/logs/bcftools-view-R-{reg}.log",


ruleorder: flt_variants_samples > flt_variants > flt_samples

# allele-specific transformation

rule vcf_atomize:
    # this rule decompose complex variants and multi-allele SNVs to individual
    # VCF lines, adding TYPE and F_MISSING to INFO field
    threads: 2
    input:
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
        ref=refseq_file,
    output:
        "{pfx}-atom/vcfs/{reg}.vcf.gz",
    log:
        log1="{pfx}-atom/logs/norm-a-m-snps-f-{reg}.log",
        log2="{pfx}-atom/logs/fill-tags-TYPE-F_MISSING-{reg}.log",
    shell:
        "bcftools norm -a -m -snps -f {input.ref} {input.vcf} 2> {log.log1}"
        " | bcftools +fill-tags -o {output} -- -t TYPE,F_MISSING,AC,AF,AN 2> {log.log2}"

rule vcf_split_snv_join:
    # this rule split (not decompose) complex variants and multi-allele SNVs to individual
    # VCF lines, adding TYPE and F_MISSING to INFO field
    threads: 3
    input:
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
        ref=refseq_file,
    output:
        vcf="{pfx}-split-snv-join/vcfs/{reg}.vcf.gz",
    log:
        log1="{pfx}-split-snv-join/logs/norm-m-any-f-{reg}.log",
        log2="{pfx}-split-snv-join/logs/bcftools-view-v-snps-{reg}.log",
        log3="{pfx}-split-snv-join/logs/mea_vcf-process-duplicate-AC-reorder-{reg}.log",
        log4="{pfx}-split-snv-join/logs/join-{reg}.log",
        log5="{pfx}-split-snv-join/logs/fill-tags-{reg}.log",

    shell:
        "bcftools norm -m -any {input.vcf} 2> {log.log1}"
        " | "
        "bcftools view -v snps  2> {log.log2}"
        " | "
        "{cli} vcf-process-duplicate --key AC --action reorder  2> {log.log3}"
        " | "
        "bcftools norm -m +any  2> {log.log4}"
        " | "
        "bcftools +fill-tags -- -t TYPE,F_MISSING,AC,AF,AN  2> {log.log5}"
        " | "
        "bcftools view -o {output.vcf}"


rule vcf_split_snv_dedup:
    # this rule split (not decompose) complex variants and multi-allele SNVs to individual
    # VCF lines
    threads: 3
    input:
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
        ref=refseq_file,
    output:
        vcf="{pfx}-split-snv-dedup/vcfs/{reg}.vcf.gz",
    log:
        log1="{pfx}-split-snv-dedup/logs/norm-m-any-f-{reg}.log",
        log2="{pfx}-split-snv-dedup/logs/bcftools-view-v-snps-{reg}.log",
        log3="{pfx}-split-snv-dedup/logs/mea_vcf-process-duplicate-AC-deduplicate-{reg}.log",
        log4="{pfx}-split-snv-dedup/logs/fill-tags-{reg}.log",

    shell:
        "bcftools norm -m -any {input.vcf}  2> {log.log1}"
        " | "
        "bcftools view -v snps  2> {log.log2}"
        " | "
        "{cli} vcf-process-duplicate  --key AC  --action deduplicate  2> {log.log3}"
        " | "
        "bcftools +fill-tags -- -t TYPE,F_MISSING,AC,AF,AN  2> {log.log4}"
        " | "
        "bcftools view -o {output.vcf}"


rule vcf_split:
    threads: 2
    input:
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
        ref=refseq_file,
    output:
        "{pfx}-split/vcfs/{reg}.vcf.gz",
    log:
        log1="{pfx}-split/logs/norm-m-any-f-{reg}.log",
        log2="{pfx}-split/logs/fill-tags-TYPE-F_MISSING-{reg}.log",
    shell:
        "bcftools norm -m -any -f {input.ref} {input.vcf}  2> {log.log1}"
        " | bcftools +fill-tags -o {output} -- -t TYPE,F_MISSING,AC,AF,AN 2> {log.log2}"


rule vcf_join:
    # reorder the duplicate variants so that alternate variants with higher AC comes first
    threads: 3
    input:
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
    output:
        vcf="{pfx}-join/vcfs/{reg}.vcf.gz",
    log:
        log1="{pfx}-join/logs/mea_vcf-process-duplicate-AC-reorder-{reg}.log",
        log2="{pfx}-join/logs/join-{reg}.log",
        log3="{pfx}-join/logs/fill-tags-{reg}.log",
    shell:
        "{cli} vcf-process-duplicate --key AC --action reorder {input.vcf}  2> {log.log1}"
        " | bcftools norm -m +any  2> {log.log2}"
        " | bcftools +fill-tags -- -t TYPE,F_MISSING,AC,AF,AN  2> {log.log3}"
        " | bcftools view -o {output.vcf}"


rule flt_snv:
    # this rule filter-in only single nucleotide variants
    threads: 2
    input:
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
    output:
        vcf="{pfx}-snv/vcfs/{reg}.vcf.gz",
    log:
        log1="{pfx}-snv/logs/bcftools-view-v-snps-{reg}.log",
    shell:
        "bcftools view -v snps -o {output.vcf} {input.vcf}  2> {log.log1}"


rule vcf_dedup:
    # this rule remove duplicated positions and keeping the variant that
    # has bigger minor allele count
    # combination of -atom-dedup is identical to mark indels and other alternate
    # alleles (eg. non-biallelic) as missing values
    threads: 3
    input:
        vcf = "{pfx}/vcfs/{reg}.vcf.gz"
    output:
        vcf = "{pfx}-dedup/vcfs/{reg}.vcf.gz"
    log:
        log1 = "{pfx}-dedup/logs/mea_vcf-process-duplicate-AC-deduplicate-{reg}.log"
    shell:
        "{cli} vcf-process-duplicate  -o {output.vcf}  --key AC  --action deduplicate"
        "  {input.vcf} 2> {log.log1}"


rule vcf_split_dedup:
    threads: 3
    input:
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
        ref=refseq_file,
    output:
        vcf="{pfx}-split-dedup/vcfs/{reg}.vcf.gz",
    log:
        log1="{pfx}-split-dedup/logs/norm-m-any-f-{reg}.log",
        log2="{pfx}-split-dedup/logs/bcftools-view-v-snps-{reg}.log",
        log3="{pfx}-split-dedup/logs/mea_vcf-process-duplicate-AC-deduplicate-{reg}.log",
        log4="{pfx}-split-dedup/logs/fill-tags-{reg}.log",

    shell:
        "bcftools norm -m -any {input.vcf}  2> {log.log1}"
        " | "
        "{cli} vcf-process-duplicate --key AC --action deduplicate  2> {log.log3}"
        " | "
        "bcftools +fill-tags -- -t TYPE,F_MISSING,AC,AF,AN  2> {log.log4}"
        " | "
        "bcftools view -o {output.vcf}"


ruleorder: vcf_split_snv_join > vcf_split_snv_dedup > vcf_split_dedup > vcf_split > vcf_join > flt_snv > vcf_dedup

# Minor alelele-based filtering


rule flt_MAC:
    threads: 3
    input:
        "{pfx}/vcfs/{reg}.vcf.gz",
    output:
        "{pfx}-MAC{MAC}/vcfs/{reg}.vcf.gz",
    shell:
        'bcftools view -e "MAC<{wildcards.MAC}" -o {output} {input}'


rule flt_MAF:
    # this rule filters variants based on minimum Minor Allele Frequency (MAF)
    threads: 3
    input:
        "{pfx}/vcfs/{reg}.vcf.gz",
    output:
        "{pfx}-MAF{MAF}/vcfs/{reg}.vcf.gz",
    shell:
        # 'bcftools +fill-tags {input} -- -t AF '
        # '| bcftools view -e "MAF<{wildcards.MAF}" -o {output} {input}'
        "bcftools +fill-tags {input} -- -t AF "
        "| bcftools view -q {wildcards.MAF}:minor  -o {output}"


# QC-based filtering


rule flt_V:
    # filter variants based on missingness threshold
    threads: 3
    input:
        variant="{pfx}/variants_{vt}.txt",
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
        vcf_idx="{pfx}/vcfs/{reg}.vcf.gz.csi",
    output:
        "{pfx}-V{vt}/vcfs/{reg}.vcf.gz",
    shell:
        "bcftools view -R {input.variant} {input.vcf} "
        '| bcftools view -e "F_MISSING>{wildcards.vt}" -o {output} '


rule prep_V:
    # prepare variant list for certain threshold
    localrule: True
    input:
        variant="{pfx}/QC/variants.tsv",
    output:
        variant="{pfx}/variants_{vt}.txt",
    run:
        import pandas as pd

        df_variant = pd.read_table(input.variant)
        df_variant[df_variant.F_MISSING <= float(wildcards.vt)][["CHROM", "POS"]].to_csv(
            output.variant, sep="\t", index=False, header=False
        )


rule flt_S:
    # filter samples based on missingness threshold
    threads: 3
    input:
        sample="{pfx}/samples_{st}.txt",
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
        vcf_idx="{pfx}/vcfs/{reg}.vcf.gz.csi",
    output:
        "{pfx}-S{st}/vcfs/{reg}.vcf.gz",
    shell:
        "bcftools view -S {input.sample} -o {output} {input.vcf}"


rule prep_S:
    # prepare sample list for certain threshold
    localrule: True
    input:
        sample="{pfx}/QC/samples.tsv",
    output:
        sample="{pfx}/samples_{st}.txt",
    run:
        import pandas as pd

        df_sample = pd.read_table(input.sample)
        df_sample[df_sample.F_MISSING <= float(wildcards.st)][["SAMPLE"]].to_csv(
            output.sample, index=False, header=False
        )


rule flt_VS:
    # filter variants and samples at the same time
    threads: 3
    input:
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
        vcf_idx="{pfx}/vcfs/{reg}.vcf.gz.csi",
        variants="{pfx}/variants_{vt}.txt",
        samples="{pfx}/samples_{st}.txt",
    output:
        vcf="{pfx}-V{vt}_S{st}/vcfs/{reg}.vcf.gz",
    shell:
        "bcftools view -R {input.variants} {input.vcf} "
        "| bcftools view -S {input.samples} -o {output.vcf}"


rule prep_VPOI:
    threads: 1
    input:
        variants="{pfx}/QC/variants.tsv",
    output:
        variants="{pfx}/variants_POI.txt",
        plot="{pfx}/variants_POI.png",
        yaml="{pfx}/variants_POI.yaml",
    shell:
        "{cli} select-by-POI  --slope {POI_slope}"
        "  -o {output.variants}  --outplot {output.plot}  --outyaml {output.yaml}"
        "  {input.variants}:CHROM,POS,MISSING"


rule flt_VPOI:
    threads: 3
    input:
        variants="{pfx}/variants_POI.txt",
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
        vcf_idx="{pfx}/vcfs/{reg}.vcf.gz.csi",
    output:
        vcf="{pfx}-VPOI/vcfs/{reg}.vcf.gz",
    shell:
        "bcftools view -R {input.variants} -o {output.vcf} {input.vcf}"


rule prep_SPOI:
    threads: 1
    input:
        samples="{pfx}/QC/samples.tsv",
    output:
        samples="{pfx}/samples_POI.txt",
        plot="{pfx}/samples_POI.png",
        yaml="{pfx}/samples_POI.yaml",
    shell:
        "{cli} select-by-POI  --slope {POI_slope}"
        "  -o {output.samples}  --outplot {output.plot}  --outyaml {output.yaml}"
        "  {input.samples}:SAMPLE,MISSING"


rule flt_SPOI:
    threads: 3
    input:
        samples="{pfx}/samples_POI.txt",
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
        vcf_idx="{pfx}/vcfs/{reg}.vcf.gz.csi",
    output:
        "{pfx}-SPOI/vcfs/{reg}.vcf.gz",
    shell:
        "bcftools view -S {input.samples} -o {output} {input.vcf}"


# -- QC related --


rule check_QC:
    threads: 1
    input:
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
    output:
        sample="{pfx}/QC/{reg}-samples.tsv",
        variant="{pfx}/QC/{reg}-variants.tsv",
    shell:
        "{cli} vcf-generate-QC-metrics --outsample {output.sample} --outvariant {output.variant} {input.vcf}"


rule check_complete_QC:
    threads: 1
    input:
        vcf="{pfx}/{reg}.vcf.gz",
    output:
        sample="{pfx}/{reg}.samples-qc.tsv",
        variant="{pfx}/{reg}.variants-qc.tsv",
    log:
        "{pfx}/logs/check_single_QC-{reg}.log",
    benchmark:
        "{pfx}/logs/benchmark-check_single_QC-{reg}.txt"
    shell:
        "{cli} vcf-generate-QC-metrics --fraction --outsample {output.sample} --outvariant {output.variant} {input.vcf} 2> {log}"


rule plot_complete_QC_png:
    threads: 1
    input:
        "{pfx}/{fname}.{var}-qc.tsv",
    output:
        "{pfx}/{fname}.{var}-qc.png",
    run:
        import pandas as pd
        import seaborn as sns
        from matplotlib import pyplot as plt

        df = pd.read_table(input[0], sep="\t")
        # plot
        ax = sns.scatterplot(
            x=range(len(df)), y=df.F_MISSING.sort_values(), s=3, linewidth=0, alpha=1.0
        )
        ax.set_ylabel(f"F_MISSING")
        ax.set_xlabel(f"{wildcards.var} index (n={len(df)})")
        plt.savefig(output[0])
        plt.close()


use rule plot_complete_QC_png as plot_complete_QC_pdf with:
    output:
        "{pfx}/{reg}.{var}-qc.pdf",


rule gather_sample_QC_png:
    localrule: True
    input:
        expand("{{pfx}}/QC/{reg}-samples.tsv", reg=REGIONS),
    output:
        outfile="{pfx}/QC/samples.tsv",
        outplot="{pfx}/QC/samples.png",
    run:
        import pandas as pd
        import seaborn as sns
        from matplotlib import pyplot as plt

        # gather
        dfs = [pd.read_table(infile, sep="\t") for infile in input]
        df = pd.concat(dfs).groupby("SAMPLE").agg("sum").reset_index()
        df["F_MISSING"] = df.MISSING / df.N_VARIANTS
        df.to_csv(output.outfile, sep="\t", index=False)

        # plot
        ax = sns.scatterplot(
            x=range(len(df)), y=df.F_MISSING.sort_values(), s=3, linewidth=0, alpha=1.0
        )
        ax.set_ylabel("Variant F_MISSING")
        ax.set_xlabel("Sample Index")
        plt.savefig(output.outplot)
        plt.close()


use rule gather_sample_QC_png as rule_gather_sample_QC_pdf with:
    output:
        outplot="{pfx}/QC/samples.pdf",


rule gather_variant_QC_png:
    localrule: True
    input:
        expand("{{pfx}}/QC/{reg}-variants.tsv", reg=REGIONS),
    output:
        outfile="{pfx}/QC/variants.tsv",
        outplot="{pfx}/QC/variants.png",
    run:
        import pandas as pd
        import seaborn as sns
        from matplotlib import pyplot as plt

        #import IPython; IPython.embed()
        dfs = [pd.read_table(infile, sep="\t") for infile in input]
        df = pd.concat(dfs)
        df["F_MISSING"] = df.MISSING / df.N_SAMPLES
        df.to_csv(output.outfile, sep="\t", index=False)

        # plot
        ax = sns.scatterplot(
            x=range(len(df)), y=df.F_MISSING.sort_values(), s=3, linewidth=0, alpha=1.0
        )
        ax.set_ylabel("Sample F_MISSING")
        ax.set_xlabel("Variant Index")
        plt.savefig(output.outplot)
        plt.close()


use rule gather_variant_QC_png as gather_variant_QC_pdf with:
    output:
        outplot="{pfx}/QC/variants.pdf",


# additional utilities


rule concat_vcfs:
    # concatenate VCF files into a single VCF file
    # requires:
    #    REGIONS
    #    complete_region [the name to assign for the new VCF file]
    threads: 2
    input:
        expand("{{pfx}}/vcfs/{reg}.vcf.gz", reg=REGIONS)
    output:
        f"{{pfx}}/concat/{complete_region}.vcf.gz"
    shell:
        "bcftools concat -o {output} {input}"


# EOF
