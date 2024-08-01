# vcfutils.py
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


import sys
import pathlib
import numpy as np
from cyvcf2 import VCF, Writer
from mea_pipeline import cout, cerr, cexit


def set_GT(
    infile: str | pathlib.Path,
    outfile: str | pathlib.Path,
    *,
    logfile: str | pathlib.Path | None = None,
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
):

    # after setting GT with this function, it is advised to run:
    # bcftools +setGT -- -t q -n . -i "FORMAT/DP<{mindepth}"
    # to set min depth and recalculate parameters in INFO field

    vcf = VCF(infile)
    w = Writer(outfile, vcf)

    for header in headers:
        w.add_to_header(header)

    logs = []

    for v in vcf:

        AD = v.format("AD")

        if minimum_minor_depth >= 0 or minimum_minor_ratio >= 0:

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
                    minor_alleles[equal_minor_depths] = minor_args[equal_minor_depths, 1]

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

        else:
            # get non-hets from genotypes / GT fields
            # XXX continue here
            non_hets = []

        # for all that are not non-hets, eg. the hets, set their alleles
        for idx in (~non_hets).nonzero()[0]:
            if set_het_to_ref:
                v.genotypes[idx] = [0, 0, False]
                continue
            if set_het_to_alt:
                # WARNING: for now this only correct for biallelic variants
                v.genotypes[idx] = [1, 1, False]
                continue
            if set_het_to_missing:
                v.genotypes[idx] = [-1, -1, False]
                continue
            # sort alleles
            if (maj_allele := major_alleles[idx]) > (min_allele := minor_alleles[idx]):
                alleles = [min_allele, maj_allele, False]
            else:
                alleles = [maj_allele, min_allele, False]
            v.genotypes[idx] = alleles

        for idx in (v.gt_depths == 0).nonzero()[0]:
            if set_missing_to_ref:
                v.genotypes[idx] = [0, 0, False]
                continue
            if set_missing_to_alt:
                # WARNING: for now this only correct for biallelic variants
                v.genotypes[idx] = [1, 1, False]
                continue
            if set_missing_to_het:
                # WARNING: for now this only correct for biallelic variants
                v.genotypes[idx] = [0, 1, False]
            v.genotypes[idx] = [-1, -1, False]

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
