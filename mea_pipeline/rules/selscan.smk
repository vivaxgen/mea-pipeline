

# -- signal of selection scan soss_ --

# selscan can only work with non-missing biallelic variants only
# hence need to convert the other alternate alleles and missing variants
# into ref alleles

rule soss_prep:
    localrule: True
    # prepare the file so the steps can be run locally within selscan-{complete_region} directory
    input:
        src = '{pfx}/hmmIBD-{complete_region}/{population_notation}.independent_samples.txt'
    output:
        dst = '{pfx}/selscan-{complete_region}/{population_notation}.sample-list.txt'
    shell:
        "ln -sr {input.src} {output.dst}"

rule soss_miss2ref_vcf:
    threads: 2
    input:
        vcf = '{pfx}/{complete_region}.vcf.gz',
        idx = '{pfx}/{complete_region}.vcf.gz.tbi',
    output:
        vcf = '{pfx}/selscan-{complete_region}/{complete_region}-miss2ref.vcf.gz',
    log:
        log1 = '{pfx}/selscan-{complete_region}/logs/bcftools-norm-m-{complete_region}.log',
        log2 = '{pfx}/selscan-{complete_region}/logs/mea_vcf-process-duplicate-{complete_region}.log',
        log3 = '{pfx}/selscan-{complete_region}/logs/mea_vcf-set-GT-{complete_region}.log',
    shell:
        "{cli} vcf-set-GT  -o {output.vcf}  --set-id"
        #"  --minimum-depth 2  --min-minor-depth 2  --min-minor-ratio 0.495"
        "  --set-missing-to-ref"
        "  {input.vcf} 2> {log.log3}"


rule soss_prep_map:
    threads: 1
    input:
        vcf = '{pfx}/{complete_region}-miss2ref.vcf.gz',
    output:
        mapfile = '{pfx}/{complete_region}-maps/{reg}.map.txt',
    shell:
        "{cli} vcf-generate-mapfile  --region {wildcards.reg}"
        "  --translation-file {chromtranslation_file}  -o {output.mapfile}  {input.vcf}"


rule soss_prep_vcf:
    threads: 2
    input:
        vcf = '{pfx}/selscan-{complete_region}/{complete_region}-miss2ref.vcf.gz',
        idx = '{pfx}/selscan-{complete_region}/{complete_region}-miss2ref.vcf.gz.tbi',
        sample = '{pfx}/selscan-{complete_region}/{population_notation}.sample-list.txt',
    output:
        vcf = '{pfx}/selscan-{complete_region}/{population_notation}.vcf.gz'
    shell:
        'bcftools view  -o {output.vcf}  -S {input.sample}  {input.vcf}'
        #' | mea-pl vcf-set-GT -o {output.vcf}'
        #' --min-alt-count 2 --min-alt-ratio 0.45'
        #' --set-het-to-ref --set-missing-to-ref'
        #' --set-missing-to-ref --set-alt2-to-ref'


rule index_vcf_region:
    threads: 1
    input:
        vcf = '{pfx}/selscan-{complete_region}/{population_notation}.vcf.gz',
    output:
        idx = '{pfx}/selscan-{complete_region}/{population_notation}.vcf.gz.tbi'
    shell:
        'bcftools index --tbi {input.vcf}'

ruleorder: index_vcf_region > index_tbi

rule soss_prep_vcf_region:
    localrule: True
    #threads: 1
    input:
        vcf = '{pfx}/{population_notation}.vcf.gz',
        idx = '{pfx}/{population_notation}.vcf.gz.tbi',
    output:
        vcf = '{pfx}/{population_notation}-vcfs/{reg}.vcf.gz'
        #vcf = '{pfx}/{pop1}-+-{pop2}-vcfs/{population_notation}-{reg}.vcf.gz'
    shell:
        'bcftools view -o {output.vcf} -r {wildcards.reg} {input.vcf}'


rule soss_vcf2tped:
    # convert vcf to tpepd based on chromosomes
    threads: 1
    input:
        vcf = '{pfx}/{population_notation}.vcf.gz'
    output:
        tped = '{pfx}/{pop1}-+-{pop2}-xpehh/{population_notation}-{reg}.tped'
    shell:
        'spcli $PYS/wgs/vcf2tped.py'
        ' --translation-file {chromtranslation_file}'
        ' --region {wildcards.reg}'
        ' -o {output.tped}'
        ' {input.vcf}'


rule soss_selscan_xpehh:
    # options to use --wagh
    threads: 16
    input:
        pop1 = "{pfx}/{pop1}-vcfs/{reg}.vcf.gz",
        pop2 = "{pfx}/{pop2}-vcfs/{reg}.vcf.gz",
        #pop1 = '{pfx}/{pop1}-+-{pop2}-vcfs/{pop1}-{reg}.vcf.gz',
        #pop2 = '{pfx}/{pop1}-+-{pop2}-vcfs/{pop2}-{reg}.vcf.gz',
        mapfile = f'{{pfx}}/{complete_region}-maps/{{reg}}.map.txt',
    output:
        outfile = '{pfx}/{pop1}-+-{pop2}-xpehh/{reg}.xpehh.out',
    params:
        flags = config['selscan_xpehh_flags'],
    shell:
        'selscan --xpehh  --threads {threads}  {params.flags}'
        '  --out {wildcards.pfx}/{wildcards.pop1}-+-{wildcards.pop2}-xpehh/{wildcards.reg}'
        '  --vcf {input.pop1}  --vcf-ref {input.pop2}'
        #'  --map {input.mapfile}'
        '  --pmap'


rule soss_concat_normalize_xpehh:
    localrule: True
    input:
        expand('{{pfx}}/{{pop1}}-+-{{pop2}}-xpehh/{reg}.xpehh.out', reg=REGIONS)
    output:
        outfile = f'{{pfx}}/{{pop1}}-+-{{pop2}}-xpehh/{complete_region}.xpehh.tsv'
    params:
        normfiles = lambda w, input: ' '.join(f'{infile}.norm' for infile in input)
    shell:
        'norm --xpehh --files {input}'
        ' && '
        'mea-pl concat-tables -o {output.outfile} {params.normfiles}'


rule soss_selscan_xpnsl:
    threads: 16
    input:
        pop1 = "{pfx}/{pop1}-vcfs/{reg}.vcf.gz",
        pop2 = "{pfx}/{pop2}-vcfs/{reg}.vcf.gz",
        #pop1 = '{pfx}/{pop1}-+-{pop2}-vcfs/{pop1}-{reg}.vcf.gz',
        #pop2 = '{pfx}/{pop1}-+-{pop2}-vcfs/{pop2}-{reg}.vcf.gz',
    output:
        outfile = '{pfx}/{pop1}-+-{pop2}-xpnsl/{reg}.xpnsl.out',
    params:
        flags = config['selscan_xpnsl_flags']
    shell:
        'selscan --xpnsl  --threads {threads} {params.flags}'
        ' --out {wildcards.pfx}/{wildcards.pop1}-+-{wildcards.pop2}-xpnsl/{wildcards.reg} '
        ' --vcf {input.pop1} --vcf-ref {input.pop2}'


rule soss_concat_normalize_xpnsl:
    threads: 1
    input:
        expand('{{pfx}}/{{pop1}}-+-{{pop2}}-xpnsl/{reg}.xpnsl.out', reg=REGIONS)
    output:
        outfile = f'{{pfx}}/{{pop1}}-+-{{pop2}}-xpnsl/{complete_region}.xpnsl.tsv'
    params:
        normfiles = lambda w, input: ' '.join(f'{infile}.norm' for infile in input)
    shell:
        'norm --xpnsl --files {input}'
        ' && '
        '{cli} concat-tables -o {output.outfile} {params.normfiles}'


rule soss_annotate_signals:
    threads: 1
    input:
        tsv =  f'{{pfx}}/{complete_region}.{{score}}.tsv'
    output:
        tsv =  f'{{pfx}}/{complete_region}.{{score}}.anno.tsv'
    params:
        flags = config["mea-pl_tab-annotate-signals"]
    shell:
        "{cli} tab-annotate-signals  -o {output.tsv}"
        "  --column normxpehh  --use-id id  {params.flags}"
        "  --translation-file {chromtranslation_file}"
        "  {input.tsv}"


def get_population_label(population, sample_list_file):
    if (label := config['population_label'].get(population, "")):
        return f'{label} (N={count_file_lines(sample_list_file)})'
    return ""


rule soss_plot_xpehh:
    localrule: True
    input:
        tsv =  f'{{pfx}}/{{pop1}}-+-{{pop2}}-xpehh/{complete_region}.xpehh.anno.tsv',
        pop1 = '{pfx}/{pop1}.sample-list.txt',
        pop2 = '{pfx}/{pop2}.sample-list.txt',
    output:
        plot = f'{{pfx}}/{{pop1}}-+-{{pop2}}-xpehh/{complete_region}.xpehh.anno.png'
    params:
        y_label = 'xp-EHH Score',
        toplabel = lambda w, input: get_population_label(w.pop1, input.pop1),
        bottomlabel = lambda w, input: get_population_label(w.pop2, input.pop2),
        bedfile = f'--bedfile {fn}' if (fn := config['bedfile_highlight_manhattan']) else '',
    shell:
        "{cli} tab-plot-manhattan  -o {output.plot}  --column normxpehh"
        "  --hline  --y-label '{params.y_label}'  {params.bedfile}"
        "  --toptext '{params.toplabel}'  --bottomtext '{params.bottomlabel}'"
        "  {input.tsv}"


#import IPython; IPython.embed()

use rule soss_plot_xpehh as soss_plot_xpnsl with:
    input:
        tsv = f'{{pfx}}/{{pop1}}-+-{{pop2}}-xpnsl/{complete_region}.xpnsl.anno.tsv',
        pop1 = '{pfx}/{pop1}.sample-list.txt',
        pop2 = '{pfx}/{pop2}.sample-list.txt',
    output:
        plot = f'{{pfx}}/{{pop1}}-+-{{pop2}}-xpnsl/{complete_region}.xpnsl.anno.png'
    params:
        y_label = 'xp-nSL Score',
        toplabel = lambda w, input: get_population_label(w.pop1, input.pop1),
        bottomlabel = lambda w, input: get_population_label(w.pop2, input.pop2),
        bedfile = f'--bedfile {fn}' if (fn := config['bedfile_highlight_manhattan']) else '',


# EOF
