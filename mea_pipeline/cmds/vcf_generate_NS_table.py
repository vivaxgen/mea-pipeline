#!/usr/bin/env mea-pl
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


from mea_pipeline import arg_parser, cerr


def init_argparser():

    p = arg_parser("generate non-synonymous table from a VCF file based on GT values")

    p.add_argument("-o", "--outfile")
    p.add_argument("infile")

    return p


def vcf_generate_NS_table(args):

    from mea_pipeline.utils import vcfutils
    import yaml

    table = vcfutils.generate_NS_table(args.infile)
    yaml.dump(table, open(args.outfile, "w"))


def main(args):
    vcf_generate_NS_table(args)


# EOF
