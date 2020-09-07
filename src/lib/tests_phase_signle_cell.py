import numpy as np
from scipy.stats import mannwhitneyu

from src.lib.metric_wrapper import metric_by_selector_all


# For each interval, perform rank-sum test against all other intervals
# data2D shape [nInterval, nTime]
def _test_inverse(data2D, i, alternative="greater"):
    dataThis = data2D[i]
    dataOther = np.hstack([data2D[:i].flatten(), data2D[i+1:].flatten()])
    return mannwhitneyu(dataThis, dataOther, alternative=alternative)[1]


def _test_pair(dataA, dataB, alternative):
    return np.array([mannwhitneyu(dataA[i], dataB[i], alternative=alternative)[1] for i in range(len(dataA))])


# For each cell and each interval determine a pvalue of a rank-sum test against all other intervals
# Return array of shape [nChannel, nInterval]
def test_inverse_all_selectors(dataDB, queryDict, phaseType, metricName='mean', alternative="greater", ranges=None, haveWaiting=True, settings=None):
    settings = {} if settings is None else settings

    rez3D = metric_by_selector_all(dataDB, queryDict, phaseType, metricName, 'pr', settings, ranges=ranges, haveWaiting=haveWaiting)

    nInterval, nChannel, nTrial = rez3D.shape
    pVals2D = []
    for iChannel in range(nChannel):
        pVals2D += [[_test_inverse(rez3D[:, iChannel], iInterv, alternative=alternative) for iInterv in range(nInterval)]]

    return np.array(pVals2D)


# # For each cell determine a pvalue of a rank-sum test between two given intervals or phases
# # Returns array of shape [nChannel, 1]
# def test_difference_selectors(dataDB, queryDict, selector1, selector2, metricName='mean', alternative="greater", settings=None):
#     settings = {} if settings is None else settings
#
#     rez1 = metric_by_selector(dataDB, queryDict, metricName, 'pr', selector1, settings, channelFilter=None)
#     rez2 = metric_by_selector(dataDB, queryDict, metricName, 'pr', selector2, settings, channelFilter=None)
#
#     testAB = _test_pair(rez1, rez2, alternative)
#     testBA = _test_pair(rez2, rez1, alternative)
#
#     return np.array([testAB, testBA]).T


def pvalues_2_significant_cells(pVal2D, pTHR=0.01):
    sign2D = pVal2D < pTHR
    return [set(np.where(sign2D[:, idx])[0]) for idx in range(pVal2D.shape[1])]


# Calculate confusion matrix
def significance_confusion_matrix(signCellsByInterv):
    '''
    :param signCellsByInterv: A list of sets of indices of significant cells by interval
    :return: An array of shape [nInterval, nInterval] of the number of cells significant simultaneously in a pair of intervals
    '''

    nPhase = len(signCellsByInterv)
    confMat = np.zeros((nPhase, nPhase))
    for i in range(nPhase):
        for j in range(i, nPhase):
            confMat[i][j] = len(signCellsByInterv[i] & signCellsByInterv[j])
            confMat[j][i] = confMat[i][j]
    return confMat


def find_exclusive_sets(signCellsByInterv, setIdxs):
    allIdxs = np.hstack(setIdxs)
    if len(allIdxs) != len(set(allIdxs)):
        raise ValueError("Requested sets have overlapping indices")

    # 1. Find elements of each set
    setList = []
    for thisSetIdxs in setIdxs:
        thisSets = [signCellsByInterv[idx] for idx in thisSetIdxs]
        setList += [set().union(*thisSets)]

    # 2. For each set, exclude all other sets from it
    setListUnique = []
    for iSet in range(len(setIdxs)):
        otherSets = setList[:iSet] + setList[iSet+1:]
        otherSetsUnion = set().union(*otherSets)
        setListUnique += [setList[iSet] - otherSetsUnion]

    return setListUnique