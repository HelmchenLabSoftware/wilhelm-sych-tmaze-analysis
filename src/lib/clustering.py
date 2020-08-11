import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster

# Compute clustering given distance matrix and distance threshold
def cluster_dist_matrix(M, t):
    distTril = np.tril(M, 1)
    # linkageMatrix = linkage(distTril, method='weighted', metric='euclidean')
    # return fcluster(linkageMatrix, nCluster, criterion='maxclust')# - 1  # Original numbering starts at 1 for some reason
    linkageMatrix = linkage(distTril, method='centroid', metric='euclidean')
    return fcluster(linkageMatrix, t, criterion='distance')# - 1  # Original numbering starts at 1 for some reason