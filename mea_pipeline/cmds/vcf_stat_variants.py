#!/usr/bin/env mea-pl
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


from mea_pipeline import arg_parser, cerr, cout


def init_argparser():
    p = arg_parser("vcf-stat-variants")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument(
        "-t",
        "--targets",
        default=None,
        help="comma-separated chromosom position, in the format of CHROM:POS or CHROM:START-END",
    )
    p.add_argument("-T", "--target-file", default=None)
    p.add_argument("infile", help="VCF file")
    return p


def vcf_stat_variants(args):

    from cyvcf2 import VCF

    vcf = VCF(args.infile, threads=args.threads)
    N = len(vcf.samples)
    targets = args.targets.split(",")

    cout(f"Samples: {N}\n")

    for target in targets:
        for v in vcf(target):
            cout(
                f"{v.CHROM} {v.POS} {v.REF} {v.ALT} aaf: {v.aaf} miss: {1 - v.call_rate}\n"
            )


def main(args):
    vcf_stat_variants(args)
