#!/usr/bin/env spcli

from mea_pipeline import cout, cerr, arg_parser


def init_argparser():
    p = arg_parser(
        "concat.py -- concatenate several dataframe files into a single dataframe file"
    )
    p.add_argument("-o", "--outfile", default="outfile.txt")
    p.add_argument(
        "-c",
        "--column",
        default=None,
        help="comma-separated column names to be concatenated, default is to concatenate "
        "all columns",
    )
    p.add_argument(
        "--add-header",
        default=None,
        help="comma-separated column names to be added to final file",
    )
    p.add_argument(
        "--add-pvalue",
        default=None,
        help="add pvalue column by normalizing the score to z-score"
        " and convert to pvalue",
    )
    p.add_argument("infiles", nargs="+")
    return p


def concat_tables(args):
    import pandas as pd

    columns = args.column.split(",") if args.column else None

    dfs = []
    c = 0
    for infile in args.infiles:
        df = pd.read_table(
            infile,
            sep="\t",
            names=args.add_header.split(",") if args.add_header else None,
        )
        if columns:
            dfs.append(df.loc[:, columns])
        else:
            dfs.append(df)
        c += 1

    # combine dfs
    cerr("[Combining %d tables]" % c)
    new_df = pd.concat(dfs, ignore_index=True)

    if args.add_pvalue:

        from mea_pipeline.utils import statsutils

        new_df["pvalue"] = statsutils.raw_to_pvalue(new_df[args.add_pvalue])

    new_df.to_csv(args.outfile, sep="\t", index=False)
    cerr("[Writing to %s]" % args.outfile)


def main(args):
    concat_tables(args)


# EOF
