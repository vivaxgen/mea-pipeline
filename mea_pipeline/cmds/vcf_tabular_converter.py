#!/usr/bin/env mea-pl
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"

import sys
from mea_pipeline import arg_parser, cerr


def init_argparser():
    p = arg_parser("vcf-tabular-converter - convert VCF to any tabular format")
    p.add_argument("--format", default="index", choices=["index", "allele"])
    p.add_argument("-o", "--outfile", required=True, help="output filename")
    p.add_argument("infile", help="VCF input filename")

    return p


def vcf_tabular_converter(args):

    from mea_pipeline.utils import vcfutils

    converter = vcfutils.VCFConverter(func=vcfutils.GT_to_index)

    samples, genotypes = converter.convert(args.infile)


def main(args):
    vcf_convert_to_gds(args)


# pass
