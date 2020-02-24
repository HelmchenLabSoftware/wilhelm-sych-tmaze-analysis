import numpy as np

from scipy.cluster.hierarchy import linkage, fcluster

from src.lib.matrix_lib import set_diag_zero

def CDF(p):
    nCDF = len(p) + 1
    y = np.linspace(0, 1, nCDF)
    x = np.zeros(nCDF)
    for i in range(nCDF-1):
        x[i + 1] = x[i] + p[i]
    return x, y

def tolerance_interval(data, p):
    return np.percentile(data, p), np.percentile(data, 1-p)

def metric(t1, t2):
    ord12 = (t1 < t2).astype(float)
    return np.abs(2 * np.mean(ord12) - 1)

# Cycle data backwards by n steps
def cycle_data(x, n):
    return np.hstack((x[n:], x[:n]))

def order_func(pvec):
    nNode = pvec.shape[0]
    cdfVec = [CDF(p) for p in pvec]

    ord1D = np.zeros((nNode, nNode))
    for i in range(nNode):
        for j in range(nNode):
            ord1D[i, j] += [np.sum(cdfVec[i] - cdfVec[j]) > 0]
    return ord1D


# Compute metric from comparisons matrix
# Comparisons matrix dimensions [nTrial, nCell, nCell]
def order_metric(comparisonsByTrial):
    phat = np.mean(comparisonsByTrial, axis=0)
    metricByConn = np.abs(2 * phat - 1)
    metricByConn = set_diag_zero(metricByConn)
    return metricByConn


# Compute clustering given distance matrix and distance threshold
def cluster_dist_matrix(M, nCluster):
    distTril = np.tril(M, 1)
    linkageMatrix = linkage(distTril, method='weighted', metric='euclidean')
    return fcluster(linkageMatrix, nCluster, criterion='maxclust') - 1  # Original numbering starts at 1 for some reason
