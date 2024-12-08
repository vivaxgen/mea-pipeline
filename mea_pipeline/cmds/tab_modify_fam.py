#!/usr/bin/env mea-pl
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


import sys
from mea_pipeline import arg_parser, cerr, cexit, gzopen
from mea_pipeline.utils import tabutils
import sys


def init_argparser():
    p = arg_parser("Modify plink fam file")

    p.add_argument(
        "--chrom-translation-file", default=None, help="chromosome translation file"
    )

    p.add_argument(
        "--phenotype-file", default=None, help="phenotype file in TSV format"
    )

    p.add_argument("-o", "--outfile", help="output filename, can be the same as infile")

    p.add_argument(
        "infile",
        help="plink fam file",
    )

    return p


# .fam file format
# family_id individual_id paternal_id maternal_id sex phenotype

# phenotype file
# family_id individual_id phenotype


def tab_modify_fam(args):
    import pandas as pd

    df = pd.read_table(args.infile, header=None)
    if args.chrom_translation_file:
        chrom_translation_df = pd.read_table(args.chrom_translation_file, header=None)
        merge_df = df.merge(chrom_translation_df, on=0)
        df = merge_df[["1_y", "1_x"] + list(range(2, len(df.columns)))]

    if args.phenotype_file:
        phenotype_df = pd.read_table(args.phenotype_file, header=None)
        merge_df = df[[0, 1]].merge(phenotype_df, on=[0, 1])
        df[5] = merge_df[2]

    if args.outfile:
        df.to_csv(args.outfile, sep="\t", header=False, index=False)


def main(args):
    return tab_modify_fam(args)


# EOF
