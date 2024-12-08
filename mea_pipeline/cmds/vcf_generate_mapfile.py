#!/usr/bin/env mea-pl

from mea_pipeline import arg_parser, cerr


def init_argparser():
    p = arg_parser("vcf-generate-mapfile - generate a mapfile from a VCF file")

    p.add_argument("--region", default=None)
    p.add_argument("--translation-file", default=None)
    p.add_argument("-o", "--outfile", required=True, help="output filename")
    p.add_argument("infile", help="VCF input filename")

    return p


def vcf_generate_mapfile(args):

    from cyvcf2 import VCF

    vcf = VCF(args.infile)
    out_f = open(args.outfile, "w")
    if args.translation_file:
        translation = read_chrom_translation(args.translation_file)
    else:
        translation = None

    for v in vcf:

        if args.region and v.CHROM != args.region:
            continue

        chrom = translation[v.CHROM] if translation else v.CHROM

        out_f.write(f"{chrom}\t{v.ID}\t{int(v.POS)/10000}\t{v.POS}\n")

    out_f.close()


def read_chrom_translation(infile):

    d = {}
    f = open(infile)
    for line in f:
        tokens = line.strip().split()
        d[tokens[0]] = int(tokens[1])
    return d


def main(args):
    vcf_generate_mapfile(args)


# EOF
