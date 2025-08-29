# vcfutils.py
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"
__version__ = "2024.08.16.01"


import sys
import pathlib
import collections
import typing
from enum import Enum
import numpy as np
import numpy.typing as npt
import pandas as pd
from numba import njit, jit
from cyvcf2 import VCF, Writer
from mea_pipeline import cout, cerr, cexit


def set_genotypes(variant, indexes, value):
    for idx in indexes:
        variant.genotypes[idx] = value


# -- the following functions only work for diploid settings


def set_alt2_gt(variant, allele):
    # find positions that have allele index > 1
    # for idx in range(len(variant.genotypes)):
    #    if variant.genotypes[idx][0]
    for idx in range(len(variant.genotypes)):
        if variant.genotypes[idx][0] > 1:
            variant.genotypes[idx] = [
                allele,
                variant.genotypes[idx][1],
                variant.genotypes[idx][2],
            ]
        if variant.genotypes[idx][1] > 1:
            variant.genotypes[idx] = [
                variant.genotypes[idx][0],
                allele,
                variant.genotypes[idx][2],
            ]


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
    set_alt2_to_ref: bool = False,
    set_alt2_to_alt1: bool = False,
    set_alt2_to_missing: bool = False,
    set_id: bool = False,
    headers: list = [],
    threads: int = 1,
) -> None:

    # we rely that cyvcf2.Variant().gt_types is (gts012=False):
    #   0 for 0/0
    #   1 for any hets such as 0/1, 0/2, 1/2, etc
    #   2 for ./.
    #   3 for 1/1, 2/2, etc

    # TODO:
    # - implement set_alt2_*

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

            if len(v.ALT) == 0:
                # in case no alternate base
                major_alleles = minor_alleles = minor_depths = np.zeros(AD.shape)

            elif len(v.ALT) == 1:
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
                v.genotypes[idx] = [allele, allele, v.genotypes[idx][-1]]

            # hets are those that are not non-hets
            hets = (~non_hets).nonzero()[0]
            non_hets = non_hets.nonzero()[0]

        else:
            # get non-hets from genotypes / GT fields

            hets = (v.gt_types != gt_HET).nonzero()[0]
            non_hets = (
                (v.gt_types == gt_HOM_REF) | (v.gt_types == gt_HOM_ALT)
            ).nonzero()[0]

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
            missings = (v.gt_types == gt_MISSING).nonzero()[0]

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

        # handling other alternate alleles

        if set_alt2_to_ref:
            set_alt2_gt(v, allele=0)
        elif set_alt2_to_alt1:
            set_alt2_gt(v, allele=1)
        elif set_alt2_to_missing:
            set_alt2_gt(v, allele=-1)
            # XXX: need to check if alleles = ./? or ?/., then set to ./.

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


# QC-ing VCF based on GT values


@njit(fastmath=True)
def observe_genotypes(
    gt_types: npt.NDArray, sample_metrics: npt.NDArray, gt_list: npt.NDArray
) -> None:

    gt_list[0] = gt_list[1] = gt_list[2] = gt_list[3] = 0

    for idx in range(gt_types.shape[0]):
        gt = gt_types[idx]
        sample_metrics[idx, gt] += 1
        gt_list[gt] += 1


def generate_QC_metrics(
    infile: str | pathlib.Path, threads: int = 1, add_fraction: bool = False
) -> tuple[pd.DataFrame, pd.DataFrame]:

    VariantMetric = collections.namedtuple(
        "VariantMetric",
        ["CHROM", "POS", "N_SAMPLES", "HOM_REF", "HETS", "HOM_ALT", "MISSING"],
    )

    vcf = VCF(infile, gts012=True, threads=threads)

    samples = vcf.samples
    N_samples = len(samples)

    variant_metrics = []  # this will be list of namedtuple

    # for sample metrics, use gt_types with 0=HOM_REF, 1=HET, 2=HOM_ALT, 3=MISSING (UNKNOWN)
    # refer to: https://brentp.github.io/cyvcf2/docstrings.html
    # sample_metrics = pd.DataFrame(dict(SAMPLE=vcf.samples, HOM_REF=0, HETS=0, HOM_ALT=0, MISSING=0))
    sample_metric_values = np.zeros((N_samples, 4), dtype=np.int32)
    gt_list = np.zeros(4, dtype=np.int32)

    for v in vcf:
        # iterate over variants
        # cerr(f'{v.CHROM}:{v.POS}')

        gt_types = v.gt_types

        observe_genotypes(gt_types, sample_metric_values, gt_list)

        hom_ref, hets, hom_alt, unknown = gt_list

        if (total := hom_ref + hets + hom_alt + unknown) != N_samples:
            raise ValueError(
                f"number of allele types {total} does not match number of samples {N_samples}"
            )

        variant_metrics.append(
            VariantMetric(v.CHROM, v.POS, N_samples, hom_ref, hets, hom_alt, unknown)
        )

    sample_metrics = pd.DataFrame(
        dict(
            SAMPLE=vcf.samples,
            HOM_REF=sample_metric_values[:, 0],
            HETS=sample_metric_values[:, 1],
            HOM_ALT=sample_metric_values[:, 2],
            MISSING=sample_metric_values[:, 3],
        )
    )

    # convert to dataframe
    if any(variant_metrics):
        variant_metrics = pd.DataFrame(variant_metrics)
    else:
        variant_metrics = pd.DataFrame(
            dict(
                CHROM=[],
                POS=[],
                N_SAMPLES=[],
                HOM_REF=[],
                HETS=[],
                HOM_ALT=[],
                MISSING=[],
            )
        )

    sample_metrics["N_VARIANTS"] = (
        sample_metrics.HOM_REF
        + sample_metrics.HETS
        + sample_metrics.HOM_ALT
        + sample_metrics.MISSING
    )

    if add_fraction:

        sample_metrics["F_HOM_REF"] = sample_metrics.HOM_REF / sample_metrics.N_VARIANTS
        sample_metrics["F_HETS"] = sample_metrics.HETS / sample_metrics.N_VARIANTS
        sample_metrics["F_HOM_ALT"] = sample_metrics.HOM_ALT / sample_metrics.N_VARIANTS
        sample_metrics["F_MISSING"] = sample_metrics.MISSING / sample_metrics.N_VARIANTS

        variant_metrics["F_HOM_REF"] = (
            variant_metrics.HOM_REF / variant_metrics.N_SAMPLES
        )
        variant_metrics["F_HETS"] = variant_metrics.HETS / variant_metrics.N_SAMPLES
        variant_metrics["F_HOM_ALT"] = (
            variant_metrics.HOM_ALT / variant_metrics.N_SAMPLES
        )
        variant_metrics["F_MISSING"] = (
            variant_metrics.MISSING / variant_metrics.N_SAMPLES
        )

    return variant_metrics, sample_metrics


_break_flag = False


def count_alleles(alleles):
    d = {}
    for allele in alleles:
        for a in allele.split("/"):
            if a == ".":
                continue
            try:
                d[a] += 1
            except KeyError:
                d[a] = 1
    return d


def calculate_MAF(
    infile: str | pathlib.Path,
    threads: int = 1,
    regions: list[str] = [],
):

    vcf = VCF(infile, gts012=True, threads=threads)
    results = []

    for reg in regions:

        # sanity checks
        if "-" not in reg:
            tokens = reg.split(":")
            reg = f"{tokens[0]}:{int(tokens[1])-1}-{tokens[1]}"

        for v in vcf(region=reg):
            if _break_flag:
                break

            allele_counts = count_alleles(v.gt_bases)
            if len(allele_counts) == 1:
                # only one allele
                results.append((v.CHROM, v.POS, 0.0))
            else:
                counts = sorted(allele_counts.values())
                minor_count = counts[-2]
                results.append((v.CHROM, v.POS, minor_count / sum(counts)))

    vcf.close()
    return results


def _break():
    global _break_flag
    _break_flag = True


# deduplicate variant positions


class SimpleCounter(object):

    def __init__(self, start=0):
        self.counter = start

    def add(self, value=1):
        self.counter += value

    def __str__(self):
        return str(self.counter)


def keep_max_AC(duplicated_variants):
    variant_kept = duplicated_variants[0]
    for variant in duplicated_variants[1:]:
        if variant_kept.INFO["AC"] < variant.INFO["AC"]:
            variant_kept = variant
    return variant_kept


def keep_max_AN(duplicated_variants):
    variant_kept = duplicated_variants[0]
    for variant in duplicated_variants[1:]:
        if variant_kept.INFO["AN"] < variant.INFO["AN"]:
            variant_kept = variant
    return variant_kept


def deduplicate_variants(
    infile: str | pathlib.Path,
    outfile: str | pathlib.Path,
    strategy: str = "keep-max-AC",
    headers: list[str] = [],
):

    vcf = VCF(infile)
    w = Writer(outfile, vcf)

    for header in headers:
        w.add_to_header(header)

    dedup_count = SimpleCounter()

    match (strategy):

        case "keep-max-AC":
            func = keep_max_AC

        case "keep-max-AN":
            func = keep_max_AN

        case _:
            raise ValueError("strategy = keep-max-AC | keep-max-AN")

    curr_chrom = None
    prev_position = -1
    duplicated_variants = []

    def _check_duplicated_variant():
        if len(duplicated_variants) == 1:
            w.write_record(duplicated_variants[0])
            duplicated_variants.clear()

        else:
            w.write_record(func(duplicated_variants))
            duplicated_variants.clear()
            dedup_count.add()

        for v in vcf:

            if curr_chrom != v.CHROM:
                curr_chrom = v.CHROM
                prev_position = -1

            if prev_position > v.start:
                cexit("[Fatal Error - VCF file is not sorted by position]")

            if prev_position == v.start:
                duplicated_variants.append(v)
                continue

            # after this line, we processes the previous line
            if any(duplicated_variants):
                _check_duplicated_variant()
            duplicated_variants.append(v)
            prev_position = v.start

        if any(duplicated_variants):
            _check_duplicated_variant()

    cerr(f"[Deduplicated {dedup_count} positions]")
    vcf.close()
    w.close()


def get_INFO_AC(variant):
    return variant.INFO["AC"]


def get_INFO_AN(variant):
    return variant.INFO["AN"]


def process_duplicate_variant(
    infile: str | pathlib.Path,
    outfile: str | pathlib.Path,
    key: str = "AC",
    action: str = "deduplicate",
    headers: list[str] = [],
):

    vcf = VCF(infile)
    w = Writer(outfile, vcf)

    for header in headers:
        w.add_to_header(header)

    dedup_count = SimpleCounter()

    match (key):

        case "AC":
            key_func = lambda v: v.INFO["AC"]

        case "AN":
            key_func = lambda v: v.INFO["AN"]

        case _:
            raise ValueError("key = AC | AN")

    match (action):
        case "deduplicate":
            action_func = lambda variants, writer: writer.write_record(variants[0])

        case "reorder":
            action_func = lambda variants, writer: [
                writer.write_record(v) for v in variants
            ]

        case _:
            raise ValueError("action = deduplicate | reorder")

    curr_chrom = None
    prev_position = -1
    duplicated_variants = []

    def _process_duplicated_variant():
        if len(duplicated_variants) == 1:
            w.write_record(duplicated_variants[0])
            duplicated_variants.clear()

        else:
            duplicated_variants.sort(key=key_func, reverse=True)
            cerr(f"Process: {duplicated_variants[0].POS}")
            action_func(duplicated_variants, w)
            duplicated_variants.clear()
            dedup_count.add()

    for v in vcf:

        if curr_chrom != v.CHROM:
            curr_chrom = v.CHROM
            prev_position = -1

        if prev_position > v.start:
            cexit("[Fatal Error - VCF file is not sorted by position]")

        if prev_position == v.start:
            duplicated_variants.append(v)
            continue

        # after this line, we processes the previous line
        if any(duplicated_variants):
            _process_duplicated_variant()
        duplicated_variants.append(v)
        prev_position = v.start

    if any(duplicated_variants):
        _process_duplicated_variant()

    vcf.close()
    w.close()

    return dedup_count.counter


def GT_to_index(variant):

    # use GT to convert to alleles, assume all to be homozygotes
    return variant.gt_types[:, 0]


class VCFConverter(object):

    def __init__(self, func: typing.Callable):

        self.func = func

    def convert(
        self,
        infile: str | pathlib.Path,
        threads: int = 1,
    ):

        vcf = VCF(infile, threads=threads)
        samples = vcf.samples
        genotypes = []

        logs = []

        for variant in vcf:

            genotypes.append(self.func(variant))

        return (samples, genotypes)


# EOF
