# vcfutils.py
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


import sys
import pathlib
import numpy as np
from cyvcf2 import VCF, Writer
from mea_pipeline import cout, cerr, cexit


def set_genotypes(variant, indexes, value):
    for idx in indexes:
        variant.genotypes[indexes] = value


# after setting GT with this function, it is advised to pipe the result to:
#   bcftools +fill-tags -o {output} - -- -t AC,AC_Hom,AC_Het,AF,AN,AF_MISSING,MAF,NS,TYPE
# to recalculate parameters in INFO field (since bcftools might compute
# faster than python-based like this function)


def set_GT(
    infile: str | pathlib.Path,
    outfile: str | pathlib.Path,
    *,
    logfile: str | pathlib.Path | None = None,
    minimum_depth: int = -1,
    minimum_minor_depth: int = -1,
    minimum_minor_ratio: float = -1,
    set_het_to_ref: bool = False,
    set_het_to_alt: bool = False,
    set_het_to_missing: bool = False,
    set_missing_to_het: bool = False,
    set_missing_to_ref: bool = False,
    set_missing_to_alt: bool = False,
    set_id: bool = False,
    headers: list = [],
    threads: int = 1,
):

    # we rely that cyvcf2.Variant().gt_types is (gts012=False):
    #   0 for 0/0
    #   1 for any hets such as 0/1, 0/2, 1/2, etc
    #   2 for ./.
    #   3 for 1/1, 2/2, etc

    gt_HOM_REF = 0
    gt_HET = 1
    gt_MISSING = 2
    gt_HOM_ALT = 3

    genotype_HOM_REF = [0, 0, False]
    genotype_HET = [0, 1, False]
    genotype_MISSING = [-1, -1, False]
    genotype_HOM_ALT = [1, 1, False]

    vcf = VCF(infile, threads=threads)
    w = Writer(outfile, vcf)

    for header in headers:
        w.add_to_header(header)

    logs = []

    # prepare some flags
    USE_GT = True
    if minimum_minor_depth >= 0 or minimum_minor_ratio >= 0:
        if minimum_depth < 0:
            raise ValueError(
                "when setting for minimum_minor_depth or minimum_minor_ratio, "
                "minimum_depth also needs to be set"
            )
        USE_GT = False

    for v in vcf:

        AD = v.format("AD")

        if not USE_GT:

            if len(v.ALT) == 1:
                # we have biallelic alleles

                minor_depths = AD.min(axis=1)
                minor_args = np.argpartition(AD, -2, axis=1)
                major_alleles = minor_args[:, 1]
                minor_alleles = minor_args[:, 0]

            else:
                # we have multiple alleles, need only to take all minor depths and
                # check if 1st & 2nd minor depths is the same

                minor_part = np.partition(AD, -2, axis=1)
                minor_args = np.argsort(-AD, axis=1)

                # get the major alleles & temp minor alleles
                major_alleles = minor_args[:, 0]
                minor_alleles = minor_args[:, 1]

                # get cumulative minor depths (eg: all minor depths except the last one)
                minor_depths = minor_part[:, :-1].sum(axis=1)

                # check if 1st minor == 2nd minor, and the value is not zero
                equal_minor_depths = (minor_part[:, -2] == minor_part[:, -3]) & (
                    minor_depths > minimum_minor_depth
                )

                if any(equal_minor_depths):

                    # in each sample with equivalent minor depths, reassigned
                    # minor allele with minor_args[1] which is guarantee to be
                    # smaller allele index than minor_args[2] when the depths of
                    # those minor alleles are the same (by np.argsort(-AD))
                    minor_alleles[equal_minor_depths] = minor_args[
                        equal_minor_depths, 1
                    ]

                    # get two of the most common allele index
                    allele_indexes = (AD < minimum_minor_depth).sum(axis=0).argsort()
                    # cerr(f'Allele indexes: {allele_indexes}')

                    # sanity check if reference is major (allele_indexes[0]) or minor
                    # (allele_indexes[1])
                    if allele_indexes[0] != 0 and allele_indexes[1] != 0:
                        logs.append(
                            f"WARNING: reference allele is not part of major/minor alleles "
                            f"at variant {repr(v)}"
                        )

            # up to this point, there should be:
            # major_alleles
            # minor_alleles (after adjustment)
            # minor_depths (after adjustment)

            minor_ratios = minor_depths / v.gt_depths

            non_hets = (minor_depths < minimum_minor_depth) | (
                minor_ratios < minimum_minor_ratio
            )

            # set minor_alleles to major_alleles if minor_depths == 0;
            null_minor_indexes = minor_depths == 0
            minor_alleles[null_minor_indexes] = major_alleles[null_minor_indexes]

            # for all that are non-hets, set both alleles to major allele
            for idx in non_hets.nonzero()[0]:
                allele = major_alleles[idx]
                v.genotypes[idx] = [allele, allele, False]

            # hets are those that are not non-hets
            hets = (~non_hets).nonzero()[0]

        else:
            # get non-hets from genotypes / GT fields

            hets = v.gt_types != gt_HET
            non_hets = (v.gt_types == gt_HOM_REF) | (v.gt_types == gt_HOM_ALT)

        # -- handling hets -- #

        # for all hets, set their alleles
        if set_het_to_ref:
            set_genotypes(v, hets, genotype_HOM_REF)
        elif set_het_to_alt:
            set_genotypes(v, hets, genotype_HOM_ALT)
        elif set_het_to_missing:
            set_genotypes(v, hets, genotype_MISSING)
        elif minimum_minor_depth >= 0 or minimum_minor_ratio >= 0:
            for idx in hets:
                # sort alleles to give 0/1, 0/2 or 1/2 instead of 1/0, 2/0 or 2/1
                if (maj_allele := major_alleles[idx]) > (
                    min_allele := minor_alleles[idx]
                ):
                    alleles = [min_allele, maj_allele, False]
                else:
                    alleles = [maj_allele, min_allele, False]
                v.genotypes[idx] = alleles

        # -- handling missing -- #

        if minimum_depth >= 0:
            missings = (v.gt_depths < minimum_depth).nonzero()[0]
        else:
            missings = v.gt_types == gt_MISSING

        if set_missing_to_ref:
            set_genotypes(v, missings, genotype_HOM_REF)
        elif set_missing_to_alt:
            # WARNING: this is correct only for biallelic variants
            set_genotypes(v, missings, genotype_HOM_ALT)
        elif set_missing_to_het:
            # WARNING: this is correct only for biallelic variants
            set_genotypes(v, missings, genotype_HET)
        elif minimum_depth >= 0:
            set_genotypes(v, missings, genotype_MISSING)

        # reset genotypes
        v.genotypes = v.genotypes

        if set_id:
            v.ID = f"{v.CHROM}:{v.POS}"

        w.write_record(v)

    vcf.close()
    w.close()

    cerr(f"Processsed VCF file written to {outfile}")

    if logfile:
        with open(logfile, "w") as outlog_f:
            outlog_f.write("\n".join(logs))
        cerr(f"Log file written to {logfile}")


# EOF
