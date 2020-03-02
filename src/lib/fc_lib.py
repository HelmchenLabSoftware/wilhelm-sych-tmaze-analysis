import numpy as np

from src.lib.stat_lib import discrete_CDF
from src.lib.baseline_lib import convert_pmf

# For each pair of channels, decide which one is earlier by comparing their time courses, converted to CDF
# Result is bool [nCell, nCell] array
def cumulative_orderability_2D(data2D, baselineByCell):
    nNode, nTime = data2D.shape

    # Convert traces to probability distributions
    pvec = [convert_pmf(dataCell, baseline) for dataCell, baseline in zip(data2D, baselineByCell)]
    cdfVec = [discrete_CDF(p) for p in pvec]

    rez = np.zeros((nNode, nNode), dtype=bool)
    for i in range(nNode):
        for j in range(i+1, nNode):
            rez[i, j] = np.sum(cdfVec[i] - cdfVec[j]) > 0
            rez[j, i] = rez[i, j]

    return rez


# For each pair of channels, compute orderability metric
#   * Metric in [-1, 1]
#   * {0: order is undecided, 1: first cell later, -1: second cell later}
# DataShape [nTrial, nCell, nTime]. nCell is fixed, but nTime may vary between trials, hence data passed as list of 2D arrays
# Result is float [nCell, nCell] array
def cumulative_orderability_3D(data2DLst):
    data2DFilter = [dataTrial for dataTrial in data2DLst if np.prod(dataTrial.shape) != 0]

    # Determine baseline by cell
    baselineByCell = np.array([np.min(dataTrial, axis=1) for dataTrial in data2DFilter])

    print(baselineByCell.shape)

    baselineByCell = np.min(baselineByCell, axis=0)

    ordByTrial = np.array([cumulative_orderability_2D(dataTrial, baselineByCell) for dataTrial in data2DFilter])
    phat2D = np.mean(ordByTrial, axis=0)  # Frequency with which first cell later than second
    return np.abs(2 * phat2D - 1)
