#!/usr/bin/env mea-pl
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


import sys
from mea_pipeline import arg_parser, cerr


def init_argparser():

    p = arg_parser("generate QC metric from a VCF file based on GT values")

    p.add_argument(
        "--key", choices=["AC", "AN"], default="AC", help="key to reorder [AC]"
    )
    p.add_argument(
        "--action",
        choices=["deduplicate", "reorder"],
        default="deduplicate",
        help="action to perform on duplicated variants [deduplicate]",
    )
    p.add_argument("-o", "--outfile", default="-")

    p.add_argument(
        "infile",
        nargs="?",
        help="input VCF file (skipping this arg would default to stdin)",
    )

    return p


def vcf_process_duplicate(args):

    from mea_pipeline.utils import vcfutils

    if not args.infile:
        args.infile = "-"

    headers = ["##mea-pl_vcf-process-duplicateCmdLine= " + " ".join(sys.argv[1:])]

    count = vcfutils.process_duplicate_variant(
        args.infile, args.outfile, key=args.key, action=args.action, headers=headers
    )

    cerr(f"[Processed {count} duplicate variants]")
    cerr(f"[Output written to {args.outfile}]")


def main(args):
    vcf_process_duplicate(args)


# EOF
