#!/usr/bin/env mea-pl
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


from mea_pipeline import arg_parser, cerr


def init_argparser():

    p = arg_parser("select items with point of inflection as the boundary")

    p.add_argument(
        "--slope",
        type=float,
        default="1.01",
        help="slope to be used for base of inflection [1.01]",
    )
    p.add_argument("--outplot", default=None, help="output plot filename [None]")
    p.add_argument("--outyaml", default=None, help="output yaml for parameters [None]")
    p.add_argument("-o", "--outfile", help="output filename")
    p.add_argument("infile", help="qc file")

    return p


def select_POI(args):

    import pandas as pd
    from mea_pipeline.utils import statsutils

    columns = None
    if ":" in args.infile:
        infile, columns = args.infile.split(":")
        columns = columns.split(",")
    else:
        infile = args.infile

    df = pd.read_table(infile)
    if columns is None:
        columns = list(df.columns[:2])

    indexes_df = statsutils.select_POI(
        df.loc[:, columns], args.slope, outplot_file=args.outplot
    )

    if args.outfile:
        indexes_df = indexes_df.sort_values(by=list(indexes_df.columns))
        indexes_df.to_csv(args.outfile, index=False, header=False, sep="\t")
        cerr(f"[selected items written to {args.outfile}]")

    if args.outyaml:
        import yaml

        d = dict(
            max_values=df[columns].max(),
            slope=args.slope,
            selected_N=len(indexes_df),
            total_N=len(df),
        )
        yaml.dump(d, open(args.outyaml, "w"))
        cerr(f"[parameters written to {args.outyaml}]")


def main(args):
    select_POI(args)


# EOF
