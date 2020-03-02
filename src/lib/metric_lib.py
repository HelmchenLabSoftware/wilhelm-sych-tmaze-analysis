import numpy as np

from src.lib.fc_lib import cumulative_orderability_3D
from src.lib.matrix_lib import get_off_diag


def get_metric_by_name(metricName):
    apply_all = lambda metric, data: list(map(metric, data))

    metricDict = {
        "mean"      : lambda data : apply_all(np.nanmean, data),
        "std"       : lambda data : apply_all(np.nanstd, data),
        "stdFreq"   : lambda data : apply_all(num_non_zero_std, data),
        "avgCorr"   : lambda data : apply_all(avg_corr, data),
        "avgOrd"    : lambda data : [avg_cumulative_orderability(data)]
    }
    return metricDict[metricName]


# Compute average off-diagonal correlation coefficient
def num_non_zero_std(data2D):
    std = np.nanstd(data2D, axis=1)
    return np.mean(std > 1.0e-10)


def avg_corr(data2D):
    absCorrAll = np.abs(np.corrcoef(data2D))

    if np.any(np.isnan(absCorrAll)):
        print("Warning, some values had 0 variance")

    return np.mean(get_off_diag(absCorrAll))


def avg_cumulative_orderability(data2DLst):
    mat = cumulative_orderability_3D(data2DLst)
    return np.mean(get_off_diag(mat))
