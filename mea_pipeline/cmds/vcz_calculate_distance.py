#!/usr/bin/env mea-pl
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


import sys
from mea_pipeline import arg_parser, cerr, cexit

# this utility can do multi-processing automatically


def init_argparser():
    p = arg_parser("calculate genetic distance from zarr format (vcz)")

    p.add_argument(
        "--genotype",
        default="GT",
        choices=["GT", "major"],
        help="genotype call to use [GT]",
    )
    p.add_argument("--mindepth", type=int, default=5)
    p.add_argument("--ploidy", type=int, default=2)
    p.add_argument(
        "--method",
        default="proportional",
        choices=["proportional", "ibs"],
        help="method to calculate distance [proportional]",
    )
    p.add_argument(
        "--encode",
        default=False,
        action="store_true",
        help="encode distance with Szudzik's pairing function",
    )
    p.add_argument(
        "-o", "--outfile", default="outdist.tsv", help="output filename [outdist.tsv]"
    )
    p.add_argument("infile", help="input filename in zarr-format (vcz)")

    return p


def vcz_calculate_distance(args):

    import pandas as pd
    import sgkit
    from mea_pipeline.utils import distutils

    cerr(f"[Reading genotype data from {args.infile}]")
    ds = sgkit.load_dataset(args.infile)

    cerr(f"[Generating index allele matrix...]")
    match args.genotype:
        case "GT":
            alleles = distutils.GT_to_alleles(ds)
        case "major":
            alleles = distutils.AD_to_alleles(
                ds.call_AD.values, min_depth=args.mindepth, allele_number=args.ploidy
            )
        case _:
            cexit(f"ERR: --genotype option {args.genotype} is unknown")

    cerr(f"[Calculating pairwise distances...]")
    distm = distutils.pairwise_distances(alleles)

    if not args.encode:
        # split distm to proportional distance
        cerr("[Decoding distance matrix...]")
        d, N = distutils.decode_distmatrix(distm)
        distm = d / N

    distm_df = pd.DataFrame(data=distm, columns=ds.sample_id)
    distm_df.to_csv(args.outfile, index=False, sep="\t")
    # tabutils.write_file(args.outfile, distm_df)
    cerr(f"[Distance written to {args.outfile}]")


def main(args):
    vcz_calculate_distance(args)


# EOF
