import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import linregress, zscore, norm


def slope_f(data):
    return linregress(np.arange(data.shape[0]), data).slope


def sigmoid_f(x, L, x0, k, b):
    y = L / (1 + np.exp(-k * (x - x0))) + b
    return y


def exponential_f(x, a, b, c):
    return a * np.exp(b * x) + c


def select_POI(data, slope=1.01, outplot_file=None):

    import seaborn as sns
    import matplotlib.pyplot as plt

    value_column = data.columns[-1]

    sorted_df = data.sort_values(by=value_column).reset_index(drop=True)

    xdata = range(len(sorted_df))
    sorted_values = sorted_df.loc[:, value_column]
    ydata_norm = sorted_values / sorted_values.max()
    xdata_norm = np.array(xdata) / len(sorted_df)

    slopes = ydata_norm.rolling(max(21, int(len(ydata_norm) / 100))).apply(
        slope_f, raw=True
    ) * len(sorted_df)

    boundary_idx = -1
    for i in range(len(slopes) - 1, 0, -1):
        if slopes[i] > slope:
            continue
        boundary_idx = i + 1
        break

    if outplot_file:
        ax = sns.scatterplot(x=xdata_norm, y=ydata_norm, s=3, c="b")
        ax.set_box_aspect(1)
        plt.axvline(boundary_idx / len(slopes), c="r")
        plt.plot([0, 1], [0, 1 * slope])
        plt.savefig(outplot_file)
        plt.close()

    return sorted_df.iloc[:boundary_idx, range(len(data.columns) - 1)]


# general functions


def raw_to_pvalue(values, tail: int = 2):

    z_scores = zscore(values)

    match tail:

        case 2:
            # Two-tailed test: p-value = 2 * (1 - CDF(|z|))
            p_value = 2 * norm.sf(abs(z_scores))

        case 1:
            # One-tailed test: p-value = 1 - CDF(z) greater than or right-tailed
            p_value = norm.sf(z_scores)

        case -1:
            # one-tailed test: less than or left-tailed
            p_value = norm.sf(abs(z_scores))

        case _:
            raise ValueError("tail can only be 2, 1 or -1")

    return p_value


# EOF
