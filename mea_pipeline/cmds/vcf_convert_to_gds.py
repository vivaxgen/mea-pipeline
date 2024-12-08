#!/usr/bin/env mea-pl

from mea_pipeline import arg_parser, cerr


def init_argparser():
    p = arg_parser("vcf2gds - convert VCF to GDS format using R SeqArray")

    p.add_argument("-o", "--outfile", required=True, help="output filename")
    p.add_argument("infile", help="VCF input filename")

    return p


def vcf_convert_to_gds(args):

    import rpy2
    from rpy2.robjects.packages import importr

    cerr("[Importing R SeqArray..]")
    sa = importr("SeqArray")

    cerr(f"[Converting {args.infile} to {args.outfile}]")
    sa.seqVCF2GDS(args.infile, args.outfile)


def main(args):
    vcf_convert_to_gds(args)


# EOF
