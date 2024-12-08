#!/usr/bin/env mea-pl

from mea_pipeline import arg_parser, cerr


def init_argparser():
    p = arg_parser("gds2fws - calculate Fws from GDS file")

    p.add_argument("-o", "--outfile", required=True, help="output filename")
    p.add_argument("infile", help="GDS input filename")

    return p


def gds_calculate_fws(args):

    import rpy2
    from rpy2.robjects.packages import importr
    import pandas as pd

    cerr("[Importing R SeqArray..]")
    sa = importr("SeqArray")

    cerr("[Importing R moimix..]")
    mx = importr("moimix")

    cerr(f"[Reading {args.infile}]")
    isolates = sa.seqOpen(args.infile)

    cerr("[Calculating Fws...]")
    fws_all = mx.getFws(isolates)

    fws_df = pd.DataFrame({"SAMPLE": list(fws_all.names), "FWS": list(fws_all)})
    fws_df.to_csv(args.outfile, sep="\t", index=False)
    cerr(f"[FWS written to {args.outfile}]")


def main(args):
    gds_calculate_fws(args)


# EOF
