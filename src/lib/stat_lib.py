import numpy as np

from scipy.cluster.hierarchy import linkage, fcluster


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