import numpy as np

from scipy.cluster.hierarchy import linkage, fcluster
from scipy.stats import mannwhitneyu, wilcoxon, binom
rstest_twosided = lambda x, y : mannwhitneyu(x, y, alternative='two-sided')


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


def test_quantity(dataA, dataB, pval):
    nTest, nSamplesA = dataA.shape
    nTest, nSamplesB = dataB.shape
    
    test_func = wilcoxon if nSamplesA == nSamplesB else rstest_twosided
    
    pValByTest = [test_func(dataA[iTest], dataB[iTest])[1] for iTest in range(nTest)]
    nTestSignificant, negLogPValPop = test_quantity_from_pval(pValByTest, pval)
    
    return pValByTest, nTestSignificant, negLogPValPop


def test_quantity_from_pval(pValByTest, pval):
    # Compute number of datasets with significant differences
    nTest = len(pValByTest)
    nTestSignificant = np.sum(np.array(pValByTest) < pval)
    
    # Compute probability of seeing at least that many by chance
    binomPMF = binom.pmf(np.arange(0, nTest), nTest, pval)
    pValPop = np.sum(binomPMF[nTestSignificant:])
    negLogPValPop = np.round(-np.log10(pValPop), 1)
    return nTestSignificant, negLogPValPop
    