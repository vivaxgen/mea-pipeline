#!/usr/bin/env mea-pl
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


import sys
from mea_pipeline import cerr, cexit, snakeutils
from mea_pipeline.cmds.run_snakefile import run_snakefile


def init_argparser():

    p = snakeutils.init_argparser("Filter and process VCF filter")
    p.add_argument(
        "-q",
        "--query",
        default=[],
        action="append",
        help="query in specification notation",
    )
    p.add_argument(
        "-o",
        "--outfile",
        default=[],
        action="append",
        help="output file(s) as target",
    )

    return p


def run_vcf_query(args):

    # check arguments
    if not any(args.query):
        cexit("Please provide --query argument")

    config = dict(queries=args.query, targets=args.target, outfiles=args.outfile)

    print(config)
    if args.snakefile is None:
        args.snakefile = "vcf_query.smk"
    run_snakefile(args, config=config)


def main(args):
    run_vcf_query(args)


# EOF
