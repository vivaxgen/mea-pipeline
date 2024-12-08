#!/usr/bin/env mea-pl
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


from mea_pipeline import arg_parser, cerr


def init_argparser():
    p = arg_parser("run-fastlmm - run FastLMM")
    p.add_argument(
        "--phenotype-file",
        required=True,
        help="phenotype file in plink phenotype format",
    )
    p.add_argument(
        "--covars-file", default=None, help="covariates file in plink phenotype format"
    )
    p.add_argument(
        "--similarity-file",
        default=None,
        help="SNPs file (in PED/binary PED) format for " "similarity matrix",
    )
    p.add_argument(
        "-o",
        "--outfile",
        required=True,
        help="output file in feather format if extension is "
        ".feather, otherwise a tab-delimited text",
    )
    p.add_argument(
        "infile",
        help="test SNPs in PED/binary PED format (must have "
        "extension of .ped or .bed)",
    )
    return p


def load_snp_file(infile):

    from pysnptools.snpreader import Ped, Bed

    fn, ext = infile[:-4], infile[-4:]
    match ext:
        case ".ped":
            variant_file = Ped(fn)
        case ".bed":
            variant_file = Bed(fn)
        case _:
            raise ValueError(f"unrecognized extension: {ext}")

    return variant_file


def run_fastlmm(args):

    import numpy as np
    from fastlmm.association import single_snp, single_snp_linreg
    from matplotlib import pyplot as plt
    import fastlmm.util.util as flutil

    variant_file = load_snp_file(args.infile)
    K0 = None
    if args.similarity_file:
        K0 = load_snp_file(args.similarity_file)

    variant_data = variant_file.read()
    if args.similarity_file is None:
        results_df = single_snp_linreg(
            variant_data, args.phenotype_file, covar=args.covars_file, count_A1=False
        )
    else:
        results_df = single_snp(
            variant_data,
            args.phenotype_file,
            K0=K0,
            covar=args.covars_file,
            count_A1=False,
        )

    if args.outfile.endswith(".feather"):
        results_df.to_feather(args.outfile)
    else:
        results_df.to_csv(args.outfile, index=False, sep="\t")


def main(args):
    run_fastlmm(args)


# EOF
