#!/usr/bin/env mea-pl
# [https://github.com/vivaxgen/mea-pipeline]

__author__ = "Hidayat Trimarsanto"
__copyright__ = "(c) 2025, Hidayat Trimarsanto"
__email__ = "trimarsanto@gmail.com,hidayat.trimarsanto@menzies.edu.au"
__license__ = "MIT"

import sys
from mea_pipeline import arg_parser, cerr

"""
BIMBAM format
=============

*.geno.txt

5 <- number of samples
4 <- number of SNPs
IND, id1, id2, id3, id4, id5
rs1, AT, TT, ??, AT, AA
rs2, GG, CC, GG, CC, CG
rs3, CC, ??, ??, CG, GG
rs4, AC, CC, AA, AC, AA

*.pos.txt

rs1, 1200, 99
rs2, 4000, 99
rs3, 3320, 99
rs4, 100001750, 99

[note: add 100 000 000 per chromosome change]

*.pheno.txt

-0.1906689
0.6430579
1.0646248
0.9002399
-0.3561080

*.covar.txt

1, -0.6493242, 1.4506556
0, 0.4436516, 0.6158216
0, -1.6376914, 0.7719995
1, -1.4223895, -0.5368858
1, -1.2909029, 1.9088563

"""


def init_argparser():
    p = arg_parser("convert VCF file to BIMBAM format")
    p.add_argument(
        "--phenotype-file",
        default=None,
        help="a TSV file containing phenotype data, "
        "use filename:sample_column,phenotype_column, "
        "default will be using 1st and 2nd columns",
    )
    p.add_argument(
        "--covariant-file",
        default=None,
        help="a TSV file containing covariant data"
        "use filename:sample_column,covar1_column,covar2_column, "
        "default will be using 1st and 2nd columns",
    )
    p.add_argument(
        "-r", "--region", default=[], action="append", help="region to extract"
    )
    p.add_argument("--translation-file", required=True, help="translation file")
    p.add_argument("-o", "--outprefix", default="")
    p.add_argument("--outdir", default=".", help="output directory")
    p.add_argument("infile", help="VCF file")
    return p


def write_genotype(outdir, outprefix, n_samples, sample_df, variants, positions):
    # this write to genotype and position files
    print(f"outdir: {outdir} outprefix: {outprefix}")

    with open(outdir / (outprefix + ".geno.txt"), "w") as f:
        f.write(f"{n_samples}\n")
        f.write(f"{len(variants)}\n")
        f.write(",".join(["IND"] + sample_df["SAMPLE"].tolist()) + "\n")
        for v in variants:
            f.write(",".join(v) + "\n")
        cerr(f"write to {outprefix}.geno.txt")

    # write to position file
    with open(outdir / (outprefix + ".pos.txt"), "w") as f:
        for p in positions:
            f.write(",".join(map(str, p)) + "\n")
        cerr(f"write to {outprefix}.pos.txt")


def vcf_convert_to_bimbam(args):

    import pathlib
    import cyvcf2
    from mea_pipeline.utils import tabutils

    outdir = pathlib.Path(args.outdir)

    chrom_translation = {}
    with open(args.translation_file) as f:
        for line in f:
            chrom, num = line.strip().split()
            chrom_translation[chrom] = num

    vcf = cyvcf2.VCF(args.infile)
    samples = vcf.samples

    # merge with phenotype file and covariant file
    if args.phenotype_file:
        sample_df, _ = tabutils.join_metafile(samples, args.phenotype_file)

        if args.covariant_file:
            covars_df, _ = tabutils.join_metafile(samples, args.covariant_file)
            sample_df = sample_df.merge(covars_df)

        sample_df = sample_df.dropna()

    # set samples to be called from vcf
    vcf.set_samples(sample_df["SAMPLE"].tolist())
    n_samples = len(sample_df["SAMPLE"])

    current_chrom = None
    pos_offset = 0

    positions = []
    variants = []

    for v in vcf:

        if current_chrom is None:
            current_chrom = v.CHROM

        if current_chrom != v.CHROM:
            if any(args.region):
                if current_chrom in args.region:
                    write_genotype(
                        outdir,
                        current_chrom,
                        n_samples,
                        sample_df,
                        variants,
                        positions,
                    )
                variants = []
                positions = []
            else:
                pos_offset += 100_000_000
            current_chrom = v.CHROM

        ID = f"{v.CHROM}:{v.POS}" if v.ID is None else v.ID
        POS = v.POS + pos_offset

        positions.append((ID, POS, chrom_translation[v.CHROM]))

        # for each genotype, we extract the actual bases and combine them into a string
        # based on alphabetical order

        genotypes = []
        for gt in v.gt_bases:
            if gt == "./.":
                genotypes.append("??")
            else:
                genotypes.append("".join(sorted(gt.split("/"))))

        variants.append([ID] + genotypes)

    if any(args.region) and current_chrom in args.region:
        write_genotype(
            outdir,
            current_chrom,
            n_samples,
            sample_df,
            variants,
            positions,
        )
    else:
        write_genotype(
            outdir, args.outprefix, n_samples, sample_df, variants, positions
        )

    # write to file
    # with open(args.outprefix + ".geno.txt", "w") as f:
    #    f.write(f"{n_samples}\n")
    #    f.write(f"{len(variants)}\n")
    #    f.write(",".join(["IND"] + sample_df["SAMPLE"].tolist()) + "\n")
    #    for v in variants:
    #        f.write(",".join(v) + "\n")
    #    cerr(f"write to {args.outprefix}.geno.txt")

    # write to position file
    # with open(args.outprefix + ".pos.txt", "w") as f:
    #    for p in positions:
    #        f.write(",".join(map(str, p)) + "\n")
    #    cerr(f"write to {args.outprefix}.pos.txt")

    # write to phenotype file
    if args.phenotype_file:
        out_fn = outdir / "test.pheno.txt"
        with open(out_fn, "w") as f:
            for p in sample_df.iloc[:, 1]:
                f.write(f"{p}\n")
        cerr(f"write to {out_fn}")

    # write to covariant file
    if args.covariant_file:
        out_fn = outdir / "test.covar.txt"
        with open(out_fn, "w") as f:
            for c in sample_df.iloc[:, 2:].itertuples(index=False):
                f.write(",".join(map(str, c)) + "\n")
        cerr(f"write to {out_fn}")


def main(args):
    vcf_convert_to_bimbam(args)


# EOF
