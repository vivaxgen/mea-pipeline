#!/usr/bin/env mea-pl

from mea_pipeline import cout, cerr, arg_parser


def init_argparser():
    p = arg_parser("Create multiple manhattan plot from columnar data")
    p.add_argument("--column", required=True)
    p.add_argument("--neg-log10", default=False, action="store_true")
    p.add_argument("--neg-log", default=False, action="store_true")
    p.add_argument("--use-id", default=None)
    p.add_argument("--quantile", type=float, default=-1)
    p.add_argument("--two-ways", default=False, action="store_true")
    p.add_argument("--minvariants", type=int, default=10)
    p.add_argument("--window-size", type=int, default=10000)
    p.add_argument(
        "--translation-file", default=None, help="translation file for chromosome"
    )
    p.add_argument("-o", "--outfile", default="outplot.png")

    p.add_argument("infile")

    return p


def tab_annotate_signals(args):

    import numpy as np
    import pandas as pd

    chrom_translation = None
    if args.translation_file:
        chrom_translation = {}
        cerr(f"Reading chromosome translation file")
        with open(args.translation_file) as f:
            for line in f:
                tokens = line.split()
                chrom_translation[tokens[0]] = tokens[1]

    cerr(f"Reading {args.infile}")
    if args.infile.endswith(".feather"):
        df = pd.read_feather(args.infile)
    else:
        df = pd.read_table(args.infile)

    if args.column not in df.columns:
        cerr(f"[ERR - no column named {args.column}]")

    if args.use_id:
        tokens = df[args.use_id].str.split(":", expand=True)
        df["CHROM"] = tokens.iloc[:, 0]
        df["POS"] = tokens.iloc[:, 1].astype(int)

    # translate chromosom
    if chrom_translation:
        df["CHROM_NO"] = [chrom_translation[ctg] for ctg in df["CHROM"]]
    else:
        df["CHROM_NO"] = df["CHROM"]

    if args.neg_log10:
        column = f"-log10({args.column})"
        df[column] = np.log10(df[args.column])
    elif args.neg_log:
        column = f"-log({args.column})"
        df[column] = np.log(df[args.column])
    else:
        column = args.column

    mark_significant(df, args.column, args.quantile, args.two_ways)
    filter_signals(df, args.minvariants, args.window_size)

    df.to_csv(args.outfile, sep="\t", index=False)


def mark_significant(data, column, quantile, two_ways):

    data["SIGNIFICANT"] = 0
    if quantile > 0:
        if two_ways:
            quantile = [1.0 - quantile, quantile]
            threshold = data.loc[:, column].quantile(quantile)
            data.loc[
                (data.loc[:, column] <= threshold.iloc[0])
                | (data.loc[:, column] >= threshold.iloc[1]),
                "SIGNIFICANT",
            ] = 1
            data["LOWER_THRESHOLD"] = threshold.iloc[0]
            data["UPPER_THRESHOLD"] = threshold.iloc[1]
        else:
            threshold = data.loc[:, column].quantile(quantile)
            data.loc[data.loc[:, column] >= threshold, "SIGNIFICANT"] = 1
            data["UPPER_THRESHOLD"] = threshold


def filter_signals(data, minvariants=5, window_size=10000):

    data["SIGNAL"] = 0
    for chrom in data.CHROM.unique():
        walk_along_chromosome(data, chrom, minvariants, window_size)


def walk_along_chromosome(data, chrom, minvariants, distance):

    for idx, row in data.loc[
        (data.CHROM == chrom) & (data.SIGNIFICANT == 1), :
    ].iterrows():

        significants = data.loc[
            (data.CHROM == chrom)
            & (data.POS > row.POS - distance)
            & (data.POS < row.POS + distance)
            & (data.SIGNIFICANT),
            :,
        ]
        if len(significants) >= minvariants:
            data.loc[significants.index, "SIGNAL"] = 1


def main(args):
    tab_annotate_signals(args)


# EOF
