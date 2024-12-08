#!/usr/bin/env mea-pl
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

import sys
from mea_pipeline import arg_parser, cerr


def init_argparser():

    p = arg_parser("generate QC metric from a VCF file based on GT values")

    p.add_argument("--fraction", default=False, action="store_true")
    p.add_argument("--outsample", default=None)
    p.add_argument("--outvariant", default=None)

    p.add_argument("infile")

    return p


def vcf_generate_QC_metric(args):

    import datetime
    from mea_pipeline.utils import vcfutils

    cerr(f"[Reading VCF file {args.infile}]")

    start_time = datetime.datetime.now()
    variant_metrics, sample_metrics = vcfutils.generate_QC_metrics(
        args.infile, add_fraction=args.fraction
    )
    finish_time = datetime.datetime.now()
    cerr(f"[QC process finished in {finish_time - start_time}]")

    if args.outsample:
        sample_metrics.to_csv(args.outsample, index=False, sep="\t")
        cerr(f"[Sample QC metrics written to {args.outsample}]")

    if args.outvariant:
        variant_metrics.to_csv(args.outvariant, index=False, sep="\t")
        cerr(f"[Variant QC mterics written to {args.outvariant}]")


def main(args):
    vcf_generate_QC_metric(args)


# EOF
