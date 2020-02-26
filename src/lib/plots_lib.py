import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from IPython.display import display
from copy import deepcopy

from src.lib.pandas_lib import outer_product_df
from src.lib.classifier_lib import binary_classifier



def plot_baseline_percentile(data):
    nCell, nTrial, nTime = data.shape

    cmap = plt.get_cmap("tab10")
    perc = [2, 5, 10, 20]
    fig, ax = plt.subplots(ncols=nCell, figsize=(4 * nCell, 4))
    for i in range(nCell):
        dataThis = data[i].flatten()
        ax[i].hist(dataThis, bins='auto')

        for ip, p in enumerate(perc):
            cut = np.percentile(dataThis, p)
            ax[i].axvline(x=cut, color=cmap(ip+1), linestyle='--', label=str(p))
        ax[i].set_yscale('log')
        ax[i].legend()
    plt.show()


def orderability_plots(metricByConn, thrBi, dataLabel):
    nCell = metricByConn.shape[0]
    metricByConnOffDiag = metricByConn[~np.eye(nCell, dtype=bool)]
    meanMetricByCell = np.sum(metricByConn, axis=0) / (nCell - 1)
    meanMetricByCellSorted = np.sort(meanMetricByCell)
    meanMetricByCellSortIdx = np.argsort(meanMetricByCell)

    print("Mean order param", np.mean(metricByConnOffDiag))
    print("Number above 1%", np.mean(metricByConnOffDiag > thrBi))

    fig, ax = plt.subplots(ncols=4, figsize=(16, 4))
    fig.suptitle(dataLabel)
    ax[0].set_title("Orderability by connection")
    ax[0].hist(metricByConnOffDiag, bins='auto')
    ax[0].axvline(x=thrBi, linestyle='--', color='y', label='p=1%')
    ax[0].legend()
    ax[0].set_xlim(0, 1)
    ax[1].set_title("Orderability matrix")
    ax[1].imshow(metricByConn, vmin=0, vmax=1)
    ax[2].set_title("Mean orderability by cell, sorted")
    ax[2].plot(meanMetricByCellSorted)
    ax[3].set_title("Orderability matrix, sorted")
    ax[3].imshow(metricByConn[meanMetricByCellSortIdx][:, meanMetricByCellSortIdx], vmin=0, vmax=1)


def clustering_plots(matRescaled, metricByConn, clustering):
    nTime = matRescaled.shape[2]
    clusterSortIdxs = np.argsort(clustering)  # Node indices so that clusters appear consecutive
    nCluster = np.max(clustering)
    nodePerClusterCumul = [np.sum(clustering <= j + 1) for j in range(nCluster)]

    # Plot orderability metric matrix, and cluster separators with red lines
    fig2, ax2 = plt.subplots(ncols=2, figsize=(8, 4))
    ax2[0].imshow(metricByConn[clusterSortIdxs][:, clusterSortIdxs])
    for clusterSep in nodePerClusterCumul:
        ax2[0].axvline(x=clusterSep - 0.5, color='r', alpha=0.5)
        ax2[0].axhline(y=clusterSep - 0.5, color='r', alpha=0.5)

    # Plot average rescaled activity for each cluster
    for iCluster in range(nCluster):
        mu = np.mean(matRescaled[clustering == iCluster], axis=(0, 1)) * 50 + iCluster
        std = np.std(matRescaled[clustering == iCluster], axis=(0, 1)) * 50
        ax2[1].fill_between(np.arange(nTime), mu - std, mu + std, alpha=0.3)
        ax2[1].plot(mu, label=str(iCluster))


def table_test_metric_phase_vs_all(dataDB, phase, metricName, metricFunc):
    sweepDF = outer_product_df({
        # "mousename"   : ["m060"],
        "performance": ["Correct", "Mistake", "All"],
        "direction": ["L", "R", "All"],
        "datatype": ["raw", "high"]
    })

    keyPhase = metricName + "(" + phase + ")"
    keyAll = metricName + "(WholeTrial)"

    rezDF = sweepDF.copy()
    rezDF[keyPhase] = 0.0
    rezDF[keyAll] = 0.0
    rezDF["pVal_log10"] = 0.0
    rezDF["nTrial"] = 0

    dict_drop_value = lambda d, val: {k: v for k, v in d.items() if v != val}

    for idx, row in sweepDF.iterrows():
        queryDict = dict_drop_value(row.to_dict(), "All")

        dataAll = dataDB.get_data_from_phase(None, queryDict)
        dataPhase = dataDB.get_data_from_phase(phase, queryDict)

        metricAll = [metricFunc(data) for data in dataAll]
        metricPhase = [metricFunc(data) for data in dataPhase]

        rezDF.at[idx, keyPhase] = np.mean(metricPhase)
        rezDF.at[idx, keyAll] = np.mean(metricAll)
        rezDF.at[idx, "pVal_log10"] = np.log10(mannwhitneyu(metricAll, metricPhase)[1])
        rezDF.at[idx, "nTrial"] = len(dataAll)

    rezDF.style.set_caption("Difference in " + metricName + " for " + phase + " phase vs all trial")
    display(rezDF)


def table_binary_classification(dataDB, phase, queryDict, binaryDimension):
    binaryValSet = set(dataDB.metaDataFrames['behaviorStates'][binaryDimension])
    assert(len(binaryValSet) == 2)

    queryDictCopy = deepcopy(queryDict)

    x = []
    y = []
    for i, colName in enumerate(binaryValSet):
        queryDictCopy[binaryDimension] = colName
        data = dataDB.get_data_from_phase(phase, queryDictCopy)
        dataMetric = np.array([np.mean(d, axis=0) for d in data])

        print(queryDictCopy, dataMetric.shape)

        x += [dataMetric]
        y += [i] * len(data)

    return binary_classifier(np.vstack(x), y, nShuffle=100, pTest=0.1)