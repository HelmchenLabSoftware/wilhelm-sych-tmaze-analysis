import numpy as np

from scipy.cluster.hierarchy import linkage, fcluster
from scipy.stats import mannwhitneyu, wilcoxon, binom
rstest_twosided = lambda x, y : mannwhitneyu(x, y, alternative='two-sided')
from mesostat.stat.permtests import difference_test


# Convert discrete PDF into CDF
def discrete_CDF(p):
    nCDF = len(p) + 1
    x = np.zeros(nCDF)
    for i in range(nCDF-1):
        x[i + 1] = x[i] + p[i]
    return x


def tolerance_interval(data, p):
    return np.percentile(data, p), np.percentile(data, 1-p)

# def metric(t1, t2):
#     ord12 = (t1 < t2).astype(float)
#     return np.abs(2 * np.mean(ord12) - 1)


# Cycle data backwards by n steps
def cycle_data(x, n):
    return np.hstack((x[n:], x[:n]))


# Compute clustering given distance matrix and distance threshold
def cluster_dist_matrix(M, t):
    distTril = np.tril(M, 1)
    # linkageMatrix = linkage(distTril, method='weighted', metric='euclidean')
    # return fcluster(linkageMatrix, nCluster, criterion='maxclust')# - 1  # Original numbering starts at 1 for some reason
    linkageMatrix = linkage(distTril, method='centroid', metric='euclidean')
    return fcluster(linkageMatrix, t, criterion='distance')# - 1  # Original numbering starts at 1 for some reason


def test_quantity(dataA, dataB, pval, proxyFunc=None, nResample=1000):
    nObject, nSamplesA = dataA.shape
    nObject, nSamplesB = dataB.shape
    pValByTest = np.zeros(nObject)

    # We are testing the values directly via wilcoxon or mann-whitney u test
    if proxyFunc is None:
        test_func = wilcoxon if nSamplesA == nSamplesB else rstest_twosided
        for iObject in range(nObject):
            pValByTest[iObject] = test_func(dataA[iObject], dataB[iObject])[1]

    # We are testing a function of the values, for which true distribution is unknown
    # Instead use permutation testing
    else:
        for iObject in range(nObject):
            pless, pmore = difference_test(proxyFunc, dataA[iObject], dataB[iObject], nResample, sampleFunction="permutation")
            pValByTest[iObject] = np.min([pless, pmore])

    nObjectSignificant = np.sum(pValByTest < pval)
    negLogPValPop = binom_ccdf(nObject, nObjectSignificant, pval)
    
    return pValByTest, nObjectSignificant, negLogPValPop


def binom_ccdf(nObject, nObjectSignificant, pval):
    # Compute probability of seeing at least nObjectSignificant positive outcomes from nObject by chance
    # Given that all tests are probability pval of being true by chance
    binomPMF = binom.pmf(np.arange(0, nObject), nObject, pval)
    pValPop = np.sum(binomPMF[nObjectSignificant:])
    return -np.log10(pValPop)
    