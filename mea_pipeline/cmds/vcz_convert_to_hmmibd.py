#!/usr/bin/env spcli

from mea_pipeline import cerr, arg_parser
from mea_pipeline.utils import posutils


def init_argparser(p=None):
    p = p if p else arg_parser("Converting VCZ to hmmIBD genotype format")
    p = posutils.init_argparser(p)

    p.add_argument(
        "-s",
        "--samplefile",
        default="",
        help="A headerless text file containing sample code per single line",
    )
    p.add_argument(
        "-t",
        "--translationfile",
        default="",
        help="TSV-format chromosome translation table",
    )

    # required if we read VCF file
    p.add_argument("-d", "--mindepth", type=int, default=5)
    p.add_argument(
        "--ploidy",
        type=int,
        default=2,
        help="Ploidy of samples (in VCF file), default = 2",
    )
    p.add_argument(
        "--max_alt_alleles",
        type=int,
        default=8,
        help="Maximum number of alternate alleles (in VCF file), default = 8",
    )

    p.add_argument("-o", "--outfile", required=True, help="Output filename")

    p.add_argument("infile")

    return p


def read_chrom_translation(infile):

    d = {}
    f = open(infile)
    for line in f:
        tokens = line.strip().split()
        d[tokens[0]] = int(tokens[1])
    return d


def vcz_convert_to_hmmibd(args):

    from mea_pipeline.utils import sgio, sgutils
    import numpy as np
    import pandas as pd

    ds, _ = sgio.prepare_dataset(
        args.infile, posfile=args.posfile, samplefile=args.samplefile
    )

    cerr("[Converting alleles to hmmIBD format...]")
    variants = sgutils.get_alleles(
        sgutils._allele_for_hmmibd,
        ds,
        hetratio=-1,
        mindepth=args.mindepth,
        minaltdepth=-1,
        useGT=False,
        threads=-1,
    )

    # convert chrom name to integer representation
    trans_d = read_chrom_translation(args.translationfile)
    chrom = [trans_d[x] for x in np.array(ds.contig_id)[ds.variant_contig]]
    coordinates = np.column_stack((chrom, ds.variant_position))

    # construct dataframe to hold genotype data
    columns = ["chrom", "pos"] + list(ds.sample_id.values)
    data = np.hstack((coordinates, variants))

    df = pd.DataFrame(data, columns=columns)
    df.sort_values(by=["chrom", "pos"], inplace=True)
    df.to_csv(args.outfile, sep="\t", index=False)
    cerr(f"Input file for hmmIBD written at: {args.outfile}")


def main(args):
    vcz_convert_to_hmmibd(args)


# EOF
