#!/usr/bin/env mea-pl
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


import sys
from mea_pipeline import arg_parser, cerr


def init_argparser():

    p = arg_parser("set GT format of a VCF file")

    p.add_argument(
        "--min-minor-depth",
        default=-1,
        type=int,
        help="minimum number of either allele reads to be called hets, "
        "MalariaGEN value is 2 [-1]",
    )
    p.add_argument(
        "--min-minor-ratio",
        default=-1,
        type=float,
        help="minimum ratio of either alleles reads to be called hets, "
        "MalariaGEN value is 0.10 [-1]",
    )
    p.add_argument(
        "--minimum-depth",
        default=-1,
        type=int,
    )
    p.add_argument(
        "--set-het-to-ref",
        default=False,
        action="store_true",
        help="set GT to reference alleles for all het alleles",
    )
    p.add_argument(
        "--set-het-to-alt",
        default=False,
        action="store_true",
        help="set GT to alternate alleles for all het alleles, only correct for "
        "biallelic variants",
    )
    p.add_argument(
        "--set-het-to-missing",
        default=False,
        action="store_true",
        help="set GT to missing for all het alleles",
    )
    p.add_argument(
        "--set-missing-to-ref",
        default=False,
        action="store_true",
        help="set GT to reference alleles for all missing alleles",
    )
    p.add_argument(
        "--set-missing-to-alt",
        default=False,
        action="store_true",
        help="set GT to alternate alleles for all missing alleles, only correct for "
        "biallelic variants",
    )
    p.add_argument(
        "--set-missing-to-het",
        default=False,
        action="store_true",
        help="set GT to het alleles for all missing alleles, only correct for "
        "biallelic variants",
    )
    p.add_argument(
        "--set-id",
        default=False,
        action="store_true",
        help="set the ID of VCF to CHROM:POS",
    )
    p.add_argument("-o", "--outfile", default="-", help="file output [stdout]")
    p.add_argument("--outlog", default=None, help="output log [None]")
    p.add_argument(
        "infile",
        nargs="?",
        help="input VCF file (skipping this arg would default to stdin)",
    )

    return p


def vcf_set_GT(args):

    from mea_pipeline.utils import vcfutils

    if not args.infile:
        args.infile = "-"

    headers = ["##mea-pl_vcf-set-GTCmdLine= " + " ".join(sys.argv[1:])]

    cerr(f"[Reading file: {args.infile}]")
    vcfutils.set_GT(
        infile=args.infile,
        outfile=args.outfile,
        logfile=args.outlog,
        minimum_minor_depth=args.min_minor_depth,
        minimum_minor_ratio=args.min_minor_ratio,
        minimum_depth=args.minimum_depth,
        set_het_to_ref=args.set_het_to_ref,
        set_het_to_alt=args.set_het_to_alt,
        set_het_to_missing=args.set_het_to_missing,
        set_missing_to_ref=args.set_missing_to_ref,
        set_missing_to_alt=args.set_missing_to_alt,
        set_missing_to_het=args.set_missing_to_het,
        set_id=args.set_id,
        headers=headers,
    )
    cerr(f"[Output written to {args.outfile}]")


def main(args):
    vcf_set_GT(args)


# EOF
