#!/usr/bin/env mea-pl
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


import sys
from mea_pipeline import arg_parser, cerr


def init_argparser():
    p = arg_parser(
        "generate a table run of homozygosity (RoH) from VCF file based on GT values"
    )

    p.add_argument("--threshold", type=float, default=0.005)
    p.add_argument("--winsize", type=int, default=100, help="window size in kb")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("-o", "--outfile", required=True)
    p.add_argument("infile")

    return p


def calculate_RoH_window(roh_maps, current_het_counts, winsize, threshold):
    for sample in current_het_counts:
        mean_het = current_het_counts[sample] / winsize
        if mean_het <= threshold:
            roh_maps[sample].append(1)
        else:
            roh_maps[sample].append(0)

        current_het_counts[sample] = 0


def vcf_calculate_RoH(args):

    import pandas as pd
    from cyvcf2 import VCF

    # read file
    cerr(f"Reading {args.infile}")
    vcf = VCF(args.infile, gts012=True, threads=args.threads)

    samples = vcf.samples
    winsize = args.winsize * 1000
    roh_maps = {sample: list() for sample in samples}
    current_het_counts = {sample: 0 for sample in samples}
    current_var_count = 0

    start_pos = -1
    current_chrom = None

    counter = 1
    for v in vcf:

        if start_pos < 0:
            cerr(f"change CHROM to {v.CHROM}")
            current_chrom = v.CHROM
            start_pos = v.POS
            current_var_count = 0

        elif current_chrom != v.CHROM:
            cerr(f"change CHROM to {v.CHROM}")
            calculate_RoH_window(roh_maps, current_het_counts, winsize, args.threshold)
            start_pos = v.POS
            current_chrom = v.CHROM
            current_var_count = 0

        elif v.POS - start_pos >= winsize:
            cerr(
                f"move with start_pos {start_pos} to block {v.POS} at counter {counter}"
            )
            calculate_RoH_window(roh_maps, current_het_counts, winsize, args.threshold)
            start_pos = v.POS
            current_var_count = 0

        counter += 1
        current_var_count += 1

        for sample, gt in zip(samples, v.gt_types):
            if gt == 1:
                current_het_counts[sample] += 1

    calculate_RoH_window(roh_maps, current_het_counts, winsize, args.threshold)

    # calculate total RoH
    SAMPLE = []
    COUNT = []
    HOMOZYG = []
    ROH = []

    for sample in samples:
        roh_map = roh_maps[sample]
        SAMPLE.append(sample)
        COUNT.append(len(roh_map))
        homozyg = roh_map.count(1)
        HOMOZYG.append(homozyg)
        ROH.append(homozyg / len(roh_map))

    df = pd.DataFrame(dict(SAMPLE=SAMPLE, COUNT=COUNT, HOMOZYG=HOMOZYG, ROH=ROH))
    df.to_csv(args.outfile, sep="\t", index=False)


def main(args):
    vcf_calculate_RoH(args)


# EOF
