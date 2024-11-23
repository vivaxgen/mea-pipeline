# vcfutils.py
# [https://github.com/vivaxgen/mea-pipeline]

__copyright__ = "(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


from mea_pipeline import cerr

# the following code are taken and extended from seqpy distance.py

from math import isqrt
import numpy as np
import sklearn.metrics

# this distance module provides functions to handle proportional genetic distance
# alleles can be coded to up 9999 types

MISSING_CODE = -1


def calculate_diff_N(X, Y):
    # return pair of (diff, N) as Szudzik's pairing function

    L = X.shape[0]
    X_miss = X == MISSING_CODE
    Y_miss = Y == MISSING_CODE
    any_missing = np.count_nonzero(X_miss | Y_miss)
    one_missing = np.count_nonzero(X_miss ^ Y_miss)
    total_diff = np.count_nonzero(X != Y)

    diff = total_diff - one_missing
    N = L - any_missing

    # return as Szudzik's pairing function
    return N**2 + diff


def decode_distmatrix(distm):
    import numpy as np

    isqrt_vec = np.vectorize(isqrt)
    N = isqrt_vec(distm)
    d = distm - N**2
    return d, N


def pairwise_distances(alleles):
    return sklearn.metrics.pairwise_distances(
        alleles, metric=calculate_diff_N, n_jobs=-1, force_all_finite=True
    ).astype(int)


def AD_to_alleles(allele_depths, min_depth=5, allele_number=-1):
    """return allele array from xarray allele_depths"""

    if allele_number > 0:
        allele_depths = allele_depths[:, :, :allele_number]

    allele_depths[allele_depths < 0] = 0
    alleles = allele_depths.argmax(axis=2)

    # calculate total depths and get missing mask (True == missing)
    variant_depths = allele_depths.sum(axis=2)
    missing_mask = variant_depths < min_depth
    alleles[missing_mask] = MISSING_CODE

    return alleles.T


def any_to_alleles(genotypes, missing=np.nan):
    """return allele array from any object (char/string/numeric) genotypes"""

    translation_dict = {}
    cerr("[Finding unique allele values...]")
    items = np.unique(genotypes)
    cerr("[Building allele translation table...]")
    for idx, item in enumerate(items, 1):
        # we pass missing because the new allele matrix will originally be filled with missing codes
        if item is missing:
            continue
        translation_dict[item] = idx
    alleles = np.full(genotypes.shape, MISSING_CODE, np.uint8)
    for k, i in translation_dict.items():
        cerr(f"[Converting allele {k} to {i}...]")
        alleles[genotypes == k] = i

    return alleles


def GT_to_alleles(datastore):

    import sgkit

    alleles = sgkit.convert_call_to_index(datastore, merge=False)
    return alleles.call_genotype_index.values.T


# EOF
