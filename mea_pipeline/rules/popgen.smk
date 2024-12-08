meta_file = config["meta_file"]

# -- distance-based rules (NJ and PCoA) --


wildcard_constraints:
    alleletype="major|GT",


rule calc_gendist:
    # this rule generate proportional genetic distance matrix from zarr file
    threads: 24
    input:
        "{pfx}/{complete_region}.vcz",
    output:
        "{pfx}/gendist-{complete_region}/distance.tsv",
    shell:
        "{cli} vcz-calculate-distance -o {output} {input}"


rule calc_gendist_genotype:
    # this rule generate proportional genetic distance matrix from zarr file
    threads: 16
    input:
        vcz="{pfx}/{complete_region}.vcz",
    output:
        tsv="{pfx}/gendist-{complete_region}/distance-A@{alleletype}.tsv",
    shell:
        "{cli} vcz-calculate-distance --genotype {wildcards.alleletype} -o {output.tsv} {input.vcz}"


rule plot_njtree:
    # this rule generate NJ tree plot from genetic distance
    # requires:
    # - metafile
    threads: 1
    input:
        dist="{pfx}/gendist-{complete_region}/distance.tsv",
    output:
        tree="{pfx}/gendist-{complete_region}/njtree.pdf",
    shell:
        "Rscript {MEA_PIPELINE_BIN}/dist2treeplot.R -o {output} {input}"


rule plot_njtree_genotype:
    # this rule generate NJ tree plot from genetic distance
    # requires:
    # - metafile
    threads: 1
    input:
        dist="{pfx}/gendist-{complete_region}/distance-A@{alleletype}.tsv",
    output:
        tree="{pfx}/gendist-{complete_region}/njtree-A@{alleletype}.pdf",
    shell:
        "Rscript {MEA_PIPELINE_BIN}/dist2treeplot.R -o {output} {input}"


rule prep_dist_anno:
    threads: 1
    input:
        distmat="{pfx}/gendist-{complete_region}/distance-A@{alleletype}.tsv",
        meta_file=meta_file,
    output:
        indv="{pfx}/gendist-{complete_region}/dist-A@{alleletype}-G@{grp}.indv.tsv",
        group="{pfx}/gendist-{complete_region}/dist-A@{alleletype}-G@{grp}.group.tsv",
    shell:
        "{cli} tab-generate-annotation --metafile {input.meta_file}:SAMPLE,{wildcards.grp} "
        "--outsample {output.indv} --outgroup {output.group} {input.distmat}"


rule plot_njtree_meta:
    threads: 1
    input:
        dist="{pfx}/gendist-{complete_region}/distance-A@{alleletype}.tsv",
        indv="{pfx}/gendist-{complete_region}/dist-A@{alleletype}-G@{grp}.indv.tsv",
        group="{pfx}/gendist-{complete_region}/dist-A@{alleletype}-G@{grp}.group.tsv",
    output:
        tree="{pfx}/gendist-{complete_region}/njtree-A@{alleletype}-G@{grp}.pdf",
    shell:
        "Rscript {MEA_PIPELINE_BIN}/dist2treeplot.R -c {input.indv} -l {input.group} -o {output} {input.dist}"


rule plot_njtree_meta_type:
    threads: 1
    input:
        dist="{pfx}/gendist-{complete_region}/distance.tsv",
        indv="{pfx}/gendist-{complete_region}/dist-A@{alleletype}-G@{grp}.indv.tsv",
        group="{pfx}/gendist-{complete_region}/dist-A@{alleletype}-G@{grp}.group.tsv",
    output:
        "{pfx}/gendist-{complete_region}/njtree-A@{alleletype}-G@{grp}-T@{type}.pdf",
    shell:
        "Rscript {MEA_PIPELINE_BIN}/dist2treeplot.R -c {input.indv} -l {input.group} -o {output} "
        "-t {wildcards.type} -L {input.dist}"


ruleorder: plot_njtree_meta_type > plot_njtree_meta


rule calculate_pcoa:
    threads: 1
    input:
        dist="{pfx}/gendist-{complete_region}/distance-A@{alleletype}.tsv",
        meta_file=meta_file,
    output:
        pcoa="{pfx}/gendist-{complete_region}/pcoa-A@{alleletype}.tsv",
        txt="{pfx}/gendist-{complete_region}/pcoa-A@{alleletype}.out.txt",
    params:
        meta_file = lambda w, input: f'--metafile {input.meta_file}' if input.meta_file else ''
    shell:
        "{cli} dis-calculate-pcoa"
        "  -o {output.pcoa}  --outlog {output.txt}  {params.meta_file}"
        "  {input.dist}"

rule plot_pcoa:
    threads: 1
    input:
        pcoa="{pfx}/gendist-{complete_region}/pcoa-A@{alleletype}.tsv",
    output:
        png="{pfx}/gendist-{complete_region}/pcoa-A@{alleletype}-G@{grp}.png"
    shell:
        "{cli} tab-plot-scatter"
        "  --outplot {output.png} --hue {wildcards.grp}"
        "  --axis PC1:PC2  --axis PC1:PC3"
        "  {input.pcoa}"


# -- FWS related rules --


rule vcf2gds:
    # convert VCF to GDS using R SeqArray library driven by vcf2gds.py
    threads: 2
    input:
        f"{{pfx}}/{complete_region}.vcf.gz",
    output:
        f"{{pfx}}/{complete_region}.gds",
    shell:
        "{cli} vcf-convert-to-gds -o {output} {input}"


rule calc_fws:
    # calculate FWS from GDS file using R moimix driven by gds2fws.py
    threads: 2
    input:
        "{pfx}/{complete_region}.gds",
    output:
        "{pfx}/moimix-{complete_region}/fws.tsv",
    shell:
        "{cli} gds-calculate-fws -o {output} {input}"


rule group_fws:
    threads: 1
    input:
        table = "{pfx}/moimix-{complete_region}/fws.tsv",
        meta = meta_file,
    output:
        table = "{pfx}/moimix-{complete_region}/fws-G@{grp}.tsv",
    run:
        import pandas as pd

        fws_df = pd.read_table(input.table)
        meta_df = pd.read_table(input.meta)
        fws_meta_df = fws_df.merge(meta_df)

        fws_meta_df.to_csv(output.table, sep='\t', index=False) 


rule plot_fws:
    threads: 1
    input:
        "{pfx}/moimix-{complete_region}/fws.tsv",
    output:
        "{pfx}/moimix-{complete_region}/fws.png",
    run:
        import pandas as pd
        import seaborn as sns
        from matplotlib import pyplot as plt

        fws_df = pd.read_table(input[0])

        ax = sns.catplot(
            data=fws_meta_df,
            x=wildcards.grp,
            order=sorted(fws_meta_df[wildcards.grp].unique()),
            y="FWS",
            palette=sns.color_palette(color_palettes["xgfs_normal12"]),
        )
        plt.xticks(rotation=60, ha="right", rotation_mode="anchor")
        plt.tight_layout()
        plt.savefig(output[0])
        plt.close()


rule plot_fws_group:
    threads: 1
    input:
        "{pfx}/moimix-{complete_region}/fws.tsv",
    output:
        "{pfx}/moimix-{complete_region}/fws-G@{grp}.png",
    run:
        import pandas as pd
        import seaborn as sns
        from matplotlib import pyplot as plt

        fws_df = pd.read_table(input[0])
        meta_df = pd.read_table(meta_file)
        fws_meta_df = fws_df.merge(meta_df)

        ax = sns.catplot(
            data=fws_meta_df,
            x=wildcards.grp,
            order=sorted(fws_meta_df[wildcards.grp].unique()),
            y="FWS",
            palette=sns.color_palette(color_palettes["xgfs_normal12"]),
        )
        plt.xticks(rotation=60, ha="right", rotation_mode="anchor")
        plt.tight_layout()
        plt.savefig(output[0], dpi=600)
        plt.close()


ruleorder: plot_fws_group > plot_fws


rule fws_ecdf_group:
    threads: 1
    input:
        tsv="{pfx}/moimix-{complete_region}/fws-G@{grp}.tsv",
    output:
        plot="{pfx}/moimix-{complete_region}/fws-G@{grp}.ecdf.png"
    shell:
        "{cli} plt-distribution -o {output.plot}"
        "  --kind ecdf  --use-y-axis --hue {wildcards.grp}  --add-hline 0.95  --add-Q1"
        "  {input.tsv}:FWS"


rule fws_cat_group:
    threads: 1
    input:
        tsv="{pfx}/moimix-{complete_region}/fws-G@{grp}.tsv",
    output:
        plot="{pfx}/moimix-{complete_region}/fws-G@{grp}.cat.png"
    shell:
        "{cli} tab-plot-categories  -o {output.plot}"
        "  --hue {wildcards.grp}  --add-hline 0.95  --add-means"
        "  {input.tsv}:FWS"


rule prep_FWS:
    localrule: True
    input:
        fws=f"{{pfx}}/concat/moimix-{complete_region}/fws.tsv",
    output:
        fws="{pfx}/fws_{fws}.txt",
    run:
        import pandas as pd

        df_fws = pd.read_table(input.fws)
        df_fws[df_fws.FWS >= float(wildcards.fws)].SAMPLE.to_csv(
            output.fws,
            index=False,
            header=False,
        )


rule flt_FWS:
    # filter for Fws
    threads: 2
    input:
        vcf="{pfx}/vcfs/{reg}.vcf.gz",
        vcf_idx="{pfx}/vcfs/{reg}.vcf.gz.csi",
        samples="{pfx}/fws_{fws}.txt",
    output:
        vcf="{pfx}-FWS{fws}/vcfs/{reg}.vcf.gz",
    shell:
        "bcftools view -S {input.samples} -o {output} {input.vcf}"


# -- hmmIBD --


rule select_sample_by_population:
    localrule: True
    input:
        #"{pfx}/fws-{complete_region}/monoclonal-samples.txt"
        sample="{pfx}/{complete_region}.samples-qc.tsv",
    output:
        sample="{pfx}/hmmIBD-{complete_region}/P#{population}-G@{grp}.sample.txt",
    run:
        import pandas as pd

        df = pd.read_table(input.sample, header=None)
        df.rename(columns={0: "SAMPLE"}, inplace=True)
        meta_df = df.merge(pd.read_table(meta_file))
        pops = wildcards.population.split("+")
        selected_df = meta_df.loc[meta_df[wildcards.grp].isin(pops)]
        if len(selected_df) == 0:
            raise ValueError(f'ERR: column {wildcards.grp} does not have items {pops}')
        selected_df.SAMPLE.to_csv(output.sample, header=False, index=False)


rule vcz2hmmibd:
    threads: 1
    input:
        vcz="{pfx}/{complete_region}.vcz",
        sample="{pfx}/hmmIBD-{complete_region}/P#{population}-G@{grp}.sample.txt",
    output:
        tsv="{pfx}/hmmIBD-{complete_region}/P#{population}-G@{grp}.hmmIBD-input.tsv",
    shell:
        "{cli} vcz-convert-to-hmmibd"
        " -s {input.sample}"
        " -o {output}"
        " -d {mindepth}"
        " -t {chromtranslation_file}"
        " {input.vcz}"


checkpoint hmmibd:
    threads: 1
    input:
        tsv="{pfx}/P#{population}-G@{grp}.hmmIBD-input.tsv",
    output:
        txt="{pfx}/P#{population}-G@{grp}.hmm_fract.txt",
    log:
        log1="{pfx}/logs/hmmIBD-P#{population}-G@{grp}.log",
    shell:
        "hmmIBD -i {input.tsv} -m 1000 -n 1000 -o {wildcards.pfx}/P#{wildcards.population}-G@{wildcards.grp} 2> {log.log1}"


def split_population_notation(wildcards):
    return wildcards.population_notation.split('-+-')
    return ['abc', 'def']
    #return expand("{{pfx}}/P#{population}-G@{{grp}}.hmm_fract.txt", population=['abc', 'def'])


rule gather_ibd:
    localrule: True
    input:
        #gather_ibd_input,
        expand("{{pfx}}/{population}.hmm_fract.txt", population=split_population_notation),
    output:
        tsv="{pfx}/{population_notation}.hmm_fract.combined.tsv",
    run:
        import pandas as pd

        dfs = []
        for infile in input:
            df = pd.read_table(infile)
            population = infile.removeprefix(wildcards.pfx + '/').removesuffix('.hmm_fract.txt')
            if (label := config['population_label'].get(population, "")):
                df['Label'] = label
            else:
                raise ValueError(population)
            dfs.append(df)
        df = pd.concat(dfs)
        df.to_csv(output.tsv, sep='\t', index=False)


rule ibd_ecdf_group:
    threads: 1
    input:
        tsv="{pfx}/{population_notation}.hmm_fract.combined.tsv",
    output:
        plot="{pfx}/{population_notation}.hmm_fract.combined.ecdf.png"
    shell:
        "{cli} plt-distribution -o {output.plot}"
        "  --kind ecdf  --use-y-axis --hue Label  --add-hline 0.5  --add-Q3  --set-ymin -0.05"
        "  {input.tsv}:fract_sites_IBD"


rule ibd_cat_group:
    localrule: True
    input:
        tsv="{pfx}/{population_notation}.hmm_fract.combined.tsv",
    output:
        plot="{pfx}/{population_notation}.hmm_fract.combined.cat.png"
    shell:
        "{cli} tab-plot-categories  -o {output.plot}"
        "  --hue Label --add-hline 0.5 --add-means --set-ymin -0.05"
        "  {input.tsv}:fract_sites_IBD"


rule ibd_ecdf:
    threads: 1
    input:
        txt="{pfx}/P#{population}-G@{grp}.hmm_fract.txt",
    output:
        img="{pfx}/P#{population}-G@{grp}-T@{plotype}.png",
    shell:
        "{cli} plt-distribution -o {output} --title {wildcards.population}"
        " --kind {wildcards.plotype} --use-y-axis  --set-ymin -0.05"
        " {input}:fract_sites_IBD"


rule plot_ibd:
    threads: 1
    input:
        txt="{pfx}/{population_notation}.hmm_fract.txt",
    output:
        plot="{pfx}/{population_notation}.ibd_plot.pdf",
    shell:
        "{cli} plt-igraph -o {output.plot}"
        " {input.txt}:sample1,sample2,fract_sites_IBD"


rule plot_ibd_meta:
    threads: 1
    input:
        txt="{pfx}/{population_notation}.hmm_fract.txt",
    output:
        plot="{pfx}/{population_notation}-C@{color_grp}.ibd_plot.pdf",
    shell:
        "{cli} plt-igraph -o {output.plot}"
        " --metafile {meta_file}:{wildcards.color_grp}"
        " {input}:sample1,sample2,fract_sites_IBD"


rule generate_cluster:
    threads: 1
    input:
        txt="{pfx}/{population_notation}.hmm_fract.txt",
    output:
        cluster="{pfx}/{population_notation}.clusters.yaml",
    shell:
        "{cli} plt-igraph --outcluster {output.cluster} {input}:sample1,sample2,fract_sites_IBD"


rule independent_samples:
    localrule: True
    input:
        cluster="{pfx}/hmmIBD-{complete_region}/{population_notation}.clusters.yaml",
        sample_qc="{pfx}/{complete_region}.samples-qc.tsv",
    output:
        txt="{pfx}/hmmIBD-{complete_region}/{population_notation}.independent_samples.txt",
    shell:
        "{cli} select-independent-samples -o {output} --sample-qc {input.sample_qc} {input.cluster}"


rule calc_RoH:
    threads: 1
    input:
        vcf = f"{{pfx}}/{complete_region}.vcf.gz",
    output:
        tsv = f"{{pfx}}/generics/RoH.tsv",
    params:
        threshold = '--threshold 0.0005',
        winsize = '--winsize 100',
    shell:
        "{cli} vcf-calculate-RoH -o {output.tsv}  {params.threshold}  {params.winsize}  {input.vcf}"


use rule group_fws as group_roh with:
    input:
        table = "{pfx}/generics/RoH.tsv",
        meta = meta_file,
    output:
        table = "{pfx}/generics/RoH-G@{grp}.tsv",


rule roh_ecdf_group:
    threads: 1
    input:
        tsv="{pfx}/generics/RoH-G@{grp}.tsv",
    output:
        plot="{pfx}/generics/RoH-G@{grp}.ecdf.png"
    shell:
        "{cli} plt-distribution -o {output.plot}"
        "  --kind ecdf  --use-y-axis --hue {wildcards.grp}"
        "  {input.tsv}:ROH"


rule roh_cat_group:
    threads: 1
    input:
        tsv="{pfx}/generics/RoH-G@{grp}.tsv",
    output:
        plot="{pfx}/generics/RoH-G@{grp}.cat.png"
    shell:
        "{cli} tab-plot-categories  -o {output.plot}"
        "  --hue {wildcards.grp}  --add-means"
        "  {input.tsv}:ROH"

# EOF
