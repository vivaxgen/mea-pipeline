#!/usr/bin/env spcli

from mea_pipeline import cout, cerr, cexit, arg_parser


def init_argparser():
    p = arg_parser("Filter signals")
    p.add_argument("-o", "--outfile")
    p.add_argument("--lower-quantile", default=0.05, type=float)
    p.add_argument("--upper-quantile", default=0.95, type=float)
    p.add_argument("--span-distance", default=50000, type=int)
    p.add_argument("--min-variant", default=2, type=int)
    p.add_argument(
        "--use-id",
        default=None,
        help="use the provided column name for getting id " "and split into CHROM:POS",
    )
    p.add_argument("--column", required=True, help="column name to get the values")
    p.add_argument("--neg-log10", default=False, action="store_true")
    p.add_argument("--neg-log", default=False, action="store_true")
    p.add_argument(
        "--translation-file", default=None, help="translation file for chromosome"
    )
    p.add_argument("infile")
    return p


def walk_along_chromosome(dataframe, chrom, distance, minvariants):

    for idx, row in dataframe.loc[
        (dataframe.CHROM == chrom) & (dataframe.in_quantile == 1), :
    ].iterrows():

        in_quantiles = dataframe.loc[
            (dataframe.CHROM == chrom)
            & (dataframe.POS > row.POS - distance)
            & (dataframe.POS < row.POS + distance)
            & (dataframe.in_quantile),
            :,
        ]
        if len(in_quantiles) >= minvariants:
            dataframe.loc[idx, "as_signal"] = 1


def filter_signals(args):

    import pandas as pd
    import numpy as np

    chrom_translation = None
    if args.translation_file:
        chrom_translation = {}
        cerr(f"Reading chromosome translation file")
        with open(args.translation_file) as f:
            for line in f:
                tokens = line.split()
                chrom_translation[tokens[0]] = tokens[1]

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

    if args.neg_log10:
        column = f"-log10({args.column})"
        df[column] = np.log10(df[args.column])
    elif args.neg_log:
        column = f"-log({args.column})"
        df[column] = np.log(df[args.column])
    else:
        column = args.column

    df["in_quantile"] = 0
    quantiles = df[column].quantile([args.lower_quantile, args.upper_quantile])
    df.loc[df[column] <= quantiles[0], "in_quantile"] = 1
    df.loc[df[column] >= quantiles[1], "in_quantile"] = 1

    df["as_signal"] = 0

    df.sort_values(by=["CHROM", "POS"], inplace=True).reset_index(
        drop=True, inplace=True
    )

    for chrom in df.CHROM.unique_values():
        walk_along_chromosome(
            df,
            chrom,
            args.span_distance,
        )

    if args.outfile.endswith(".feather"):
        df.to_feather(args.outfile)
    else:
        df.to_csv(args.outfile, index=False)


def main(args):
    filter_signals(args)


# EOF
