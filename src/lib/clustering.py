import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.cluster import AffinityPropagation, SpectralClustering, OPTICS

# Compute clustering given distance matrix and distance threshold
def cluster_dist_matrix(M, t, method='hierarchic'):
    if method=='hierarchic':
        distTril = np.tril(M, 1)
        # linkageMatrix = linkage(distTril, method='weighted', metric='euclidean')
        # return fcluster(linkageMatrix, nCluster, criterion='maxclust')# - 1  # Original numbering starts at 1 for some reason
        linkageMatrix = linkage(distTril, method='centroid', metric='euclidean')
        rez = fcluster(linkageMatrix, t, criterion='distance')
    elif method == 'affinity':
        clustering = AffinityPropagation(affinity='precomputed', damping=t).fit(M)
        rez =  clustering.labels_
    elif method == 'spectral':
        clustering = SpectralClustering(affinity='precomputed', assign_labels="discretize", n_init=100).fit(M)
        rez =  clustering.labels_
    elif method == 'optics':
        clustering = OPTICS(metric='precomputed', min_samples=t).fit(M)
        rez =  clustering.labels_
    else:
        raise ValueError("Unknown method", method)

    # Original numbering may start at something other than 0 for some methods
    rez = np.array(rez, dtype=int)
    return rez - np.min(rez).astype(int)

