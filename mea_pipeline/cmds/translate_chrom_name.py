#!/usr/bin/env spcli

from mea_pipeline import cerr, cexit, arg_parser


def init_argparser():
    p = arg_parser("transform chromosome name (1st column) to numeric values")

    p.add_argument("--translation-file", required=True, help="translation file")
    p.add_argument("-o", "--outfile", default="outdata.tsv", help="output filename")
    p.add_argument("infile", help="input text file")

    return p


def translate_chrom_name(args):

    import pandas as pd

    df = pd.read_table(args.infile, header=None)
    translation_df = pd.read_table(args.translation_file, header=None)
    merge_df = df.merge(translation_df, on=0)

    # import IPython; IPython.embed()

    merge_df = merge_df[["1_y", "1_x"] + list(range(2, len(df.columns)))]
    merge_df.to_csv(args.outfile, header=False, index=False, sep="\t")


def main(args):
    translate_chrom_name(args)


# EOF
