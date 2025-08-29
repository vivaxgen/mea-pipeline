#!/usr/bin/env mea-pl
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

import sys
from mea_pipeline import arg_parser, cerr


def init_argparser():
    p = arg_parser("vcf-plot-phenotype")
    p.add_argument(
        "--phenotype-file",
        help="a TSV file containing phenotype data, "
        "use the filename:sample_column,phenotype_column to specify columns to use",
    )
    p.add_argument(
        "--size-file",
        default=None,
        help="a TSV file containing normalized values for representing dot sizes, "
        "use the filename:column specification to specify the column to use",
    )
    p.add_argument("--threads", type=int, default=1)
    p.add_argument(
        "-t",
        "--targets",
        default=None,
        help="comma-separated chromosom position, in the format of CHROM:POS or CHROM:START-END",
    )
    p.add_argument(
        "-T",
        "--target-file",
        default=None,
        help="a TSV file containing target positions, "
        "use the filename:column specification to specify the column to use",
    )
    p.add_argument(
        "--contig-pattern",
        default=None,
        help="string pattern for contig name after chromosome number, eg. 'PvP01_{:02}_v2 ",
    )
    p.add_argument("--count", type=int, default=10)
    p.add_argument("--outprefix", default="outplot")
    p.add_argument("infile", help="VCF file")
    return p


def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def vcf_plot_phenotype(args):

    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    from cyvcf2 import VCF
    from mea_pipeline.utils import tabutils

    vcf = VCF(args.infile, threads=args.threads, gts012=True)

    samples = vcf.samples

    # merge samples data with phenotype data using tabutils join_metafile
    pheno_df, diff = tabutils.join_metafile(samples, args.phenotype_file)
    pheno_column = pheno_df.columns[1]

    if args.size_file:
        size_df, diff = tabutils.join_metafile(samples, args.size_file)
        size_column = size_df.columns[1]

        if pheno_df.iloc[:, 0].equals(size_df.iloc[:, 0]):
            cerr("Phenotype and size files have different sample order")
            sys.exit(1)

    targets = []
    sort_order = 1
    if args.target_file:
        # split filename and column specification
        if ":" in args.target_file:
            target_file, target_columns = args.target_file.split(":")
            target_columns = target_columns.split(",")
            if target_columns[-1].startswith("-"):
                target_columns[-1] = target_columns[-1][1:]
                sort_order = -1

        else:
            target_file = args.target_file
            target_columns = [0, 1, 2]

        target_df = tabutils.read_file(target_file)[target_columns]
        target_df.sort_values(
            by=target_df.columns[-1],
            ascending=True if sort_order > 0 else False,
            inplace=True,
        )
        print(target_df[:5])
        if args.contig_pattern:
            contig_pattern = args.contig_pattern + ":{}"
        else:
            contig_pattern = "{}:{}"
        targets += [
            contig_pattern.format(
                int(row[0]) if is_float(row[0]) else row[0], int(row[1])
            )
            for row in target_df[: args.count].itertuples(index=False)
        ]

    if args.targets:
        targets += args.targets.split(",")

    counter = 0
    for target in targets:

        contig, pos = target.split(":")
        if "-" not in pos:
            pos = f"{pos}-{pos}"
            target = f"{contig}:{pos}"

        for v in vcf(target):
            plot_df = pd.DataFrame(
                {
                    "sample": pheno_df.iloc[:, 0],
                    # "genotype": np.array(["0/0", "0/1", "1/1", "./."])[v.gt_types],
                    "genotype": v.gt_bases,
                    pheno_column: pheno_df.loc[:, pheno_column],
                    "alt_ratio": v.gt_alt_depths / (v.gt_alt_depths + v.gt_ref_depths),
                }
            )
            plot_df = plot_df[plot_df["genotype"] != "./."]
            plot_df = plot_df.dropna()
            plot_df.sort_values(by="genotype", inplace=True)

            fg = sns.catplot(
                x="genotype",
                y=pheno_column,
                data=plot_df,
                legend=False,
                kind="strip",
                hue="alt_ratio",
            )
            fg.ax.set_yscale("log")

            for i, label in enumerate(plot_df.genotype.unique()):
                data = plot_df.loc[plot_df.genotype == label][pheno_column]
                mean = data.mean()
                q1, q2, q3 = data.quantile([0.25, 0.50, 0.75])
                for x, c in zip([mean, q2, q1, q3], ["black", "red", "green", "green"]):
                    fg.ax.plot([i - 0.2, i + 0.2], [x, x], color=c, ls="--", lw=0.5)

            if plot_df.genotype.str.len().max() > 8:
                plt.xticks(rotation=60, ha="right", rotation_mode="anchor")

            plt.title(f"{v.CHROM}:{v.POS}")
            plt.tight_layout()

            if args.outprefix:
                outplot = f"{args.outprefix}{counter:02}-{v.CHROM}_{v.POS}.png"
                outtable = f"{args.outprefix}{counter:02}-{v.CHROM}_{v.POS}.tsv"
                plot_df.to_csv(outtable, sep="\t", index=False)
                cerr(f"Saving table to {outtable} with {len(plot_df)} data points")
                cerr(f"Saving plot to {outplot}")
                plt.savefig(outplot, dpi=300)
                plt.close()

        counter += 1


def main(args):
    vcf_plot_phenotype(args)


# EOF
