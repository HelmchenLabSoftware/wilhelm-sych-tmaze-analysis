import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from mesostat.utils.matrix import tril_1D, triu_1D, pairwise_differences
import mesostat.metric.sequence as sequence
from mesostat.stat.permtests import perm_test, difference_test
from mesostat.stat.testing.quantity_test import test_quantity
from mesostat.stat.resampling import sample
from mesostat.visualization.mpl_matrix import imshow

from src.lib.clustering import cluster_dist_matrix
from src.lib.metric_wrapper import metric_by_selector
import src.lib.plots_lib as plots_lib


######################
# Helper Functions
######################

# Reorder rows and columns of a matrix based on a permutation
def mat_order_by_idx(M, idx):
    return np.copy(M[idx][:, idx])


def _abbo_diff_func(x, y):
    return sequence.avg_bivariate_binary_orderability_from_temporal_mean(x) -\
           sequence.avg_bivariate_binary_orderability_from_temporal_mean(y)


######################
# Methods
######################

def count_significant_cells(dataDB, datatype, selector, condition, metricName, nResample=1000, pval=0.01, proxyFunc=None):
    nMice = len(dataDB.mice)
    condValues = set(dataDB.metaDataFrames['behaviorStates'][condition])
    mouseResults = []

    settings = {"dropShortTrialsTHR": 5}

    for iMouse, mousename in enumerate(sorted(dataDB.mice)):
        means = []
        for condVal in condValues:
            queryDictCond = {"datatype": datatype, "mousename": mousename, condition: condVal}
            means += [metric_by_selector(dataDB, queryDictCond, metricName, "pr", selector, settings)]

        nCells = dataDB.get_nchannel(mousename, datatype)
        pValByCell, nCellSignificant, negLogPValPop = test_quantity(means[0], means[1], pval, nResample=nResample,
                                                                    proxyFunc=proxyFunc)

        mouseResults += [[mousename, nCellSignificant, nCells, np.round(negLogPValPop, 1)]]

    rezDF = pd.DataFrame(mouseResults, columns=["mouse", "nCellSignificant", "nCellTot", "-log10(pval)"])
    return rezDF


def count_significant_cells_legendre(dataDB, datatype, selector, condition, basisIdx, nResample=1000, pval=0.01, proxyFunc=None):
    nMice = len(dataDB.mice)
    condValues = set(dataDB.metaDataFrames['behaviorStates'][condition])
    mouseResults = []

    settings = {"dropShortTrialsTHR": 5}

    for iMouse, mousename in enumerate(sorted(dataDB.mice)):
        means = []
        for condVal in condValues:
            queryDictCond = {"datatype": datatype, "mousename": mousename, condition: condVal}
            means += [
                metric_by_selector(dataDB, queryDictCond, 'temporal_basis', "pr", selector, settings)[:, :, basisIdx]]

        nCells = dataDB.get_nchannel(mousename, datatype)
        pValByCell, nCellSignificant, negLogPValPop = test_quantity(means[0], means[1], pval, nResample=nResample,
                                                                    proxyFunc=proxyFunc)

        mouseResults += [[mousename, nCellSignificant, nCells, np.round(negLogPValPop, 1)]]

    rezDF = pd.DataFrame(mouseResults, columns=["mouse", "nCellSignificant", "nCellTot", "-log10(pval)"])
    #     display(rezDF)
    return rezDF


def test_orderability(dataDB, datatype, selector, condition, nResample=1000, pval=0.01, proxyFunc=None, signCellsMouseDict=None):
    nMice = len(dataDB.mice)
    condValues = list(set(dataDB.metaDataFrames['behaviorStates'][condition]))
    mouseResults = []

    settings = {"dropShortTrialsTHR": 5}
    settingsTest = {"haveEffectSize": True, "haveMeans": True}

    if signCellsMouseDict is not None:
        channelFilter = signCellsMouseDict
    else:
        channelFilter = None

    rezLst = []
    for iMouse, mousename in enumerate(sorted(dataDB.mice)):
        means = []
        for condVal in condValues:
            queryDictCond = {"datatype": datatype, "mousename": mousename, condition: condVal}
            means += [metric_by_selector(dataDB, queryDictCond, "ord_mean", "rp", selector, settings,
                                         channelFilter=channelFilter)]

        pLess, pMore, eff, mu1, mu2 = difference_test(_abbo_diff_func, means[0], means[1], nSample=nResample,
                                                      sampleFunction="permutation", settings=settingsTest)

        rezLst += [[mousename, pLess, pMore, eff, mu1, mu2]]

    columnNames = ['mousename', 'pValLess', 'pValMore', 'EffSize', 'diffTrue', 'diffRand']

    rezDF = pd.DataFrame(rezLst, columns=columnNames)
    #     display(rezDF)
    return rezDF


def plot_stability_temporal_mean(dataDB):
    nMice = len(dataDB.mice)
    fig, ax = plt.subplots(ncols=nMice, figsize=(4 * nMice, 4), tight_layout=True)

    for iMouse, mousename in enumerate(sorted(dataDB.mice)):
        queryDictThis = {'datatype': 'deconv', 'mousename': mousename}

        # Calculate orderability
        dataTempMu = metric_by_selector(dataDB,
                                        queryDictThis,
                                        "ord_mean",
                                        "pr",
                                        {'interval': 8},
                                        {},
                                        channelFilter=None)

        muTempMu = np.mean(dataTempMu, axis=1)
        stdTempMu = np.std(dataTempMu, axis=1)

        ax[iMouse].plot(np.abs(muTempMu - 0.5), stdTempMu, '.')
        ax[iMouse].set_ylabel("|MU temporal mean - 0.5|")
        ax[iMouse].set_ylabel("STD temporal mean")

    plt.show()


def plot_undirected_orderability(dataDB, queryDict, selector, nResample=100, signCellsMouseDict=None):
    nMice = len(dataDB.mice)

    figureProp = {
        "ord": (4, 4, "Orderability matrices"),
        "ordSort": (4, 4, "Orderability matrices, sorted by max orderability"),
        "ordSort1D": (4, 4, "Orderability averaged for each cell"),
        "ordCorr": (4, 4, "Absolute Correlation between sorted orderability rows"),
        "ordCorrClust": (4, 4, "Absolute Correlation, sorted by clustering"),
        "ordClust": (4, 4, "Orderability matrices, clustered by highest correlation"),
        # "tempMeanHist"   : (4, 4, "Temporal mean histogram"),
        "avgResampledTemp": (10, 10, "Average resampled activity, colored by temporal mean"),
        "avgResampledClust": (10, 10, "Average resampled activity, colored by cluster"),
        "avgResampledClustAvg": (5, 5, "Average resampled activity, colored by cluster, cell-average")
    }

    figs = {}
    for key, (xScale, yScale, title) in figureProp.items():
        fig, ax = plt.subplots(ncols=nMice, figsize=(xScale * nMice, yScale))
        fig.suptitle(title)
        figs[key] = (fig, ax)

    for iMouse, mousename in enumerate(sorted(dataDB.mice)):
        print("Doing mouse", mousename)
        queryDictThis = {**queryDict, **{'mousename': mousename}}

        ########################
        # Calculations
        ########################

        # Calculate orderability and temporal means
        dataBinOrd = metric_by_selector(dataDB,
                                        queryDictThis,
                                        "ord_binary",
                                        "",
                                        selector,
                                        {},
                                        channelFilter=signCellsMouseDict)

        dataTempMu = metric_by_selector(dataDB,
                                        queryDictThis,
                                        "ord_mean",
                                        "p",
                                        selector,
                                        {},
                                        channelFilter=signCellsMouseDict)

        # Find indices of sorted temporal order
        idxSortTempMu = np.argsort(dataTempMu)

        # Sort orderability matrix in descending order of average orderability per cell
        meanBinOrd1D = np.nanmean(dataBinOrd, axis=0)
        idxSortBinOrd = np.argsort(meanBinOrd1D)[::-1]
        dataBinOrdSorted = mat_order_by_idx(dataBinOrd, idxSortBinOrd)

        # Correlation matrix of cells
        dataBinOrdSortedZero = np.copy(dataBinOrdSorted)
        np.fill_diagonal(dataBinOrdSortedZero, 0)
        dataBinOrdSortedCorr = np.abs(np.corrcoef(dataBinOrdSortedZero))

        # Clustering of correlation matrix
        clustering = cluster_dist_matrix(dataBinOrdSortedCorr, 0.5, method='affinity')
        idxClustering = np.argsort(clustering)
        clusteringSorted = clustering[idxClustering]

        # Orderability and correlation sorted by clustering
        dataBinOrdSortedCorrClust = mat_order_by_idx(dataBinOrdSortedCorr, idxClustering)
        dataBinOrdClust = mat_order_by_idx(dataBinOrdSorted, idxClustering)

        # Resampled cell activities sorted by clustering
        cellResampledActivity = metric_by_selector(dataDB,
                                                   queryDictThis,
                                                   "resample_fixed",
                                                   "p",
                                                   selector,
                                                   {"metricSettings": {"nResamplePoint": nResample}},
                                                   channelFilter=signCellsMouseDict)

        cellResampledActivitySortCluster = cellResampledActivity[idxClustering]
        cellResampledActivitySortTemporal = cellResampledActivity[idxSortTempMu]

        ########################
        # Plotting
        ########################
        getfigs = lambda figname: (figs[figname][0], figs[figname][1][iMouse])

        imshow(*getfigs("ord"), dataBinOrd, title=mousename, haveColorBar=True, limits=[0, 1])

        imshow(*getfigs("ordSort"), dataBinOrdSorted, title=mousename, haveColorBar=True, limits=[0, 1])

        figs["ordSort1D"][1][iMouse].plot(meanBinOrd1D[idxSortBinOrd])

        imshow(*getfigs("ordCorr"), dataBinOrdSortedCorr, title=mousename, haveColorBar=True, limits=[0, 1])

        imshow(*getfigs("ordCorrClust"), dataBinOrdSortedCorrClust, title=mousename, haveColorBar=True, limits=[0, 1])
        plots_lib.clustering_imshow_overplot_lines(figs["ordCorrClust"][1][iMouse], clustering)

        imshow(*getfigs("ordClust"), dataBinOrdClust, title=mousename, haveColorBar=True, limits=[0, 1])
        plots_lib.clustering_imshow_overplot_lines(figs["ordClust"][1][iMouse], clustering)

        imshow(*getfigs("avgResampledTemp"), cellResampledActivitySortTemporal, title=mousename, haveColorBar=True)

        imshow(*getfigs("avgResampledClust"), cellResampledActivitySortCluster, title=mousename, haveColorBar=True)
        plots_lib.clustering_imshow_overplot_lines(figs["avgResampledClust"][1][iMouse], clustering, haveV=False)

        xResample = np.arange(nResample)
        for cluster in set(clustering):
            thisClusterData = cellResampledActivity[np.array(clustering) == cluster]
            muTrace = np.mean(thisClusterData, axis=0)
            stdTrace = np.std(thisClusterData, axis=0)
            figs["avgResampledClustAvg"][1][iMouse].plot(xResample, muTrace, label=cluster)
            figs["avgResampledClustAvg"][1][iMouse].fill_between(xResample, muTrace - stdTrace, muTrace + stdTrace,
                                                                 alpha=0.3)

        figs["avgResampledClustAvg"][1][iMouse].legend()

    #         nCluster = len(set(clustering))
    #         colors = random_colors(nCluster)
    #         for iCell, cellAct in enumerate(cellResampledActivity):
    #             thisColor = colors[clusteringSorted[iCell]]
    #             figs["avgResampledAct"][1][iMouse].plot(cellAct + iCell, color=thisColor)

    fnameSuffix = '_'.join([
        queryDict['datatype'],
        list(selector.keys())[0],
        list(selector.values())[0],
    ])

    # Save all figures
    for key, (fig, ax) in figs.items():
        plt.figure(fig.number)
        plt.savefig(key + '_' + fnameSuffix + '.pdf', dpi=600)


def plot_directed_orderability(dataDB, datatype, selector, signCellsSelector=None):
    # Determine cell filtering
    if signCellsSelector == None:
        signCellsSelector = {'None': None}

    signCellsName, signCellsMouseDict = list(signCellsSelector.items())[0]
    selectorName, selectorVal = list(selector.items())[0]

    nMice = len(dataDB.mice)
    figureProp = {
        "ord": (4, 4, "Directed Orderability matrices"),
        "ordSort": (4, 4, "Directed Orderability matrices, sorted by max directed orderability"),
        "ordSortShuffleNeuron": (4, 4, "Directed Orderability matrices, sorted by max directed orderability, shuffled by neurons"),
        "ordSortShuffleTrial": (4, 4, "Directed Orderability matrices, sorted by max directed orderability, shuffled by trials"),
        "ordMuScatter": (4, 4, "Relationship of average temporal mean and average directed orderability"),
        "ordShiftInv": (4, 4, "Relationship of difference of temporal means and directed orderability"),
    }

    figs = {}
    for key, (xScale, yScale, title) in figureProp.items():
        fig, ax = plt.subplots(ncols=nMice, figsize=(xScale * nMice, yScale))
        fig.suptitle(title)
        figs[key] = (fig, ax)

    for iMouse, mousename in enumerate(sorted(dataDB.mice)):
        print("Doing mouse", mousename)
        queryDictThis = {'datatype': datatype, 'mousename': mousename}

        ########################
        # Calculations
        ########################

        # Calculate orderability and temporal means
        settings = {"metricSettings": {"directed": True}, "dropShortTrialsTHR": 4}
        dataBinOrd = metric_by_selector(dataDB,
                                        queryDictThis,
                                        "ord_binary",
                                        "",
                                        selector,
                                        settings,
                                        channelFilter=signCellsMouseDict)

        dataTempMu = metric_by_selector(dataDB,
                                        queryDictThis,
                                        "ord_mean",
                                        "rp",
                                        selector,
                                        settings,
                                        channelFilter=signCellsMouseDict)

        # Sort orderability matrix in descending order of average orderability per cell
        meanBinOrd1D = np.nanmean(dataBinOrd, axis=0)
        idxSortBinOrd = np.argsort(meanBinOrd1D)[::-1]
        dataBinOrdSorted = mat_order_by_idx(dataBinOrd, idxSortBinOrd)

        #         print(idxSortBinOrd)

        # Calculate orderability permuted by neurons
        dataTempMuShuffleNeurons = sample(dataTempMu, 'permutation', permAxis=1, iterAxis=0)
        dataBinOrdShuffleNeurons = sequence.bivariate_orderability_from_temporal_mean(dataTempMuShuffleNeurons, {"directed": True})
        meanBinOrd1DShuffleNeurons = np.nanmean(dataBinOrdShuffleNeurons, axis=0)
        idxSortBinOrdShuffleNeurons = np.argsort(meanBinOrd1DShuffleNeurons)[::-1]
        dataBinOrdSortedShuffleNeurons = mat_order_by_idx(dataBinOrdShuffleNeurons, idxSortBinOrdShuffleNeurons)

        # Calculate orderability permuted by trials
        dataTempMuShuffleTrials = sample(dataTempMu, 'permutation', permAxis=0, iterAxis=1)
        dataBinOrdShuffleTrials = sequence.bivariate_orderability_from_temporal_mean(dataTempMuShuffleTrials, {"directed": True})
        meanBinOrd1DShuffleTrials = np.nanmean(dataBinOrdShuffleTrials, axis=0)
        idxSortBinOrdShuffleTrials = np.argsort(meanBinOrd1DShuffleTrials)[::-1]
        dataBinOrdSortedShuffleTrials = mat_order_by_idx(dataBinOrdShuffleTrials, idxSortBinOrdShuffleTrials)

        # Calculate difference matrix of temporal means
        dataTempMuAvg = np.mean(dataTempMu, axis=0)
        diffTempMu = pairwise_differences(dataTempMuAvg)
        diffTempMu1D = triu_1D(diffTempMu)
        dataBinOrd1D = triu_1D(dataBinOrd)

        ########################
        # Plotting
        ########################
        getfigs = lambda figname: (figs[figname][0], figs[figname][1][iMouse])

        imshow(*getfigs("ord"), dataBinOrd, title=mousename, aspect=None, haveTicks=True, haveColorBar=True,
               limits=[-1, 1], cmap='jet')

        imshow(*getfigs("ordSort"), dataBinOrdSorted, title=mousename, aspect=None, haveTicks=True, haveColorBar=True,
               limits=[-1, 1], cmap='jet')

        imshow(*getfigs("ordSortShuffleNeuron"), dataBinOrdSortedShuffleNeurons, title=mousename, aspect=None,
               haveTicks=True, haveColorBar=True, limits=[-1, 1], cmap='jet')

        imshow(*getfigs("ordSortShuffleTrial"), dataBinOrdSortedShuffleTrials, title=mousename, aspect=None,
               haveTicks=True, haveColorBar=True, limits=[-1, 1], cmap='jet')

        getfigs("ordMuScatter")[1].plot(dataTempMuAvg, meanBinOrd1D, '.')
        getfigs("ordMuScatter")[1].set_xlim([0, 1])
        getfigs("ordMuScatter")[1].set_ylim([-1, 1])
        getfigs("ordMuScatter")[1].axhline(y=0, linestyle='--', color='r')
        getfigs("ordMuScatter")[1].axvline(x=0.5, linestyle='--', color='r')

        getfigs("ordShiftInv")[1].plot(diffTempMu1D, dataBinOrd1D, '.')
        getfigs("ordShiftInv")[1].axhline(y=0, linestyle='--', color='r')
        getfigs("ordShiftInv")[1].axvline(x=0, linestyle='--', color='r')

    suffix = '_'.join([datatype, selectorName, str(selectorVal), signCellsName])

    # Save all figures
    for key, (fig, ax) in figs.items():
        plt.figure(fig.number)
        plt.savefig("directed_" + key + '_' + suffix + '.pdf', dpi=600)
        plt.close()


def plots_directed_orderability_subphase_comparison(dataDB, datatype, phase):
    nMice = len(dataDB.mice)
    semiPhases = dataDB.get_phasetype_keys_from_phase(phase, 'semiphase', 'Correct')
    selectorPairs = [('phase', phase)] + [('semiphase', semiphase) for semiphase in semiPhases]
    nRows = len(selectorPairs) + 3

    fig, ax = plt.subplots(nrows=nRows, ncols=nMice, figsize=(4 * nMice, 4 * nRows))
    fig.suptitle("Comparison of directed orderability for phases and intervals")

    for iMouse, mousename in enumerate(sorted(dataDB.mice)):
        print("Doing mouse", mousename)
        queryDictThis = {'datatype': datatype, 'mousename': mousename}
        ax[0, iMouse].set_title(mousename)

        ord1Ddict = {}
        for iSelector, selectorPair in enumerate(selectorPairs):
            ########################
            # Calculations
            ########################

            # Calculate orderability and temporal means
            settings = {"metricSettings": {"directed": True}, "dropShortTrialsTHR": 4}
            dataBinOrd = metric_by_selector(dataDB,
                                            queryDictThis,
                                            "ord_binary",
                                            "",
                                            {selectorPair[0] : selectorPair[1]},
                                            settings)

            if iSelector == 0:
                temporalMeanAvg = metric_by_selector(dataDB,
                                                    queryDictThis,
                                                    "ord_mean",
                                                    "rp",
                                                    {selectorPair[0] : selectorPair[1]},
                                                    settings)
                temporalMeanAvg = np.mean(temporalMeanAvg, axis=0)

            ord1Ddict[selectorPair[1]] = np.abs(tril_1D(dataBinOrd))

            # Sort orderability matrix in descending order of average orderability per cell
            # All data will be sorted based on phase orderability.
            if iSelector == 0:
                meanBinOrd1D = np.nanmean(dataBinOrd, axis=0)
                idxSortBinOrd = np.argsort(meanBinOrd1D)[::-1]

            dataBinOrdSorted = mat_order_by_idx(dataBinOrd, idxSortBinOrd)
            meanBinOrd1DSorted = np.nanmean(dataBinOrdSorted, axis=0)


            ########################
            # Plotting
            ########################
            imshow(fig, ax[iSelector, iMouse], dataBinOrdSorted,
                   aspect=None, haveTicks=True, haveColorBar=True, limits=[-1, 1], cmap='jet')

            if iMouse == 0:
                ax[iSelector, 0].set_ylabel(selectorPair[1])

            if iSelector == 0:
                ax[-3, iMouse].plot(temporalMeanAvg[idxSortBinOrd])

            ax[-1, iMouse].plot(meanBinOrd1DSorted, label=selectorPair[1])


        ord1Ddf = pd.DataFrame(ord1Ddict)
        sns.violinplot(ax=ax[-2, iMouse], data=ord1Ddf, cut=0)
        ax[-1, iMouse].legend()


    plt.savefig("directed_comparison_" + datatype + '_' + phase + '.pdf', dpi=600)
    plt.close()


def write_cell_avg_orderability(dataDB, outfname):
    rezDict = {}
    for interval in [6, 7, 8]:
        for datatype in ['raw', 'deconv']:
            for mousename in sorted(dataDB.mice):
                key = (datatype, mousename, interval)
                print(key)

                queryDictThis = {'datatype': datatype, 'mousename': mousename}

                # Calculate orderability
                dataBinOrd = metric_by_selector(dataDB,
                                                queryDictThis,
                                                "ord_binary",
                                                "",
                                                {'interval': interval},
                                                {},
                                                channelFilter=None)

                rezDict[str(key)] = [list(row) for row in dataBinOrd]

    with open(outfname, 'w') as fp:
        json.dump(rezDict, fp)


def plot_two_cell_temporal_means_bytrial(dataDB, datatype, mousename, selector, cellIdx1, cellIdx2):
    queryDict = {"datatype": datatype, "mousename": mousename}

    # get data for all cells
    data = dataDB.get_data_from_selector(selector, queryDict)

    # extract data for cells 1 and 2
    data1 = [d[cellIdx1] for d in data]
    data2 = [d[cellIdx2] for d in data]

    settings = {"metricSettings": {"directed": True}, "dropShortTrialsTHR": 4}
    dataTempMu = metric_by_selector(dataDB,
                                    queryDict,
                                    "ord_mean",
                                    "pr",
                                    selector,
                                    settings,
                                    channelFilter=None)

    muCell1 = dataTempMu[cellIdx1]
    muCell2 = dataTempMu[cellIdx2]

    # plot series
    shiftPref = 1.5
    plt.figure(figsize=(7, 40))
    for shIdx, (d1, d2, mu1, mu2) in enumerate(zip(data1, data2, muCell1, muCell2)):
        plt.plot(d1 + shiftPref * shIdx, color='blue')
        plt.plot(d2 + shiftPref * shIdx, color='green')

        plt.plot([mu1 * len(d1), mu1 * len(d1)], [shiftPref * shIdx, shiftPref * shIdx + 1], color='darkblue')
        plt.plot([mu2 * len(d2), mu2 * len(d2)], [shiftPref * shIdx, shiftPref * shIdx + 1], color='darkgreen')

    # plot temporal means for both
    plt.figure()
    plt.plot(muCell1, muCell2, '.')
    plt.xlabel("mu_" + str(cellIdx1))
    plt.ylabel("mu_" + str(cellIdx2))
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.axhline(y=0.5, linestyle='--', color='r')
    plt.axvline(x=0.5, linestyle='--', color='r')
    plt.title('Temporal means')
    plt.show()


def test_average_orderability(dataDB, datatype, phasetype, performance,
                              signCellsSelector=None, haveWaiting=False, withSEM=False, permStrategy='neuron'):

    # Determine cell filtering
    if signCellsSelector == None:
        signCellsSelector = {'None': None}

    # Determine permutation axes
    if permStrategy == 'neuron':
        permAxis = 1
        iterAxis = 0
    elif permStrategy == 'trial':
        permAxis = 0
        iterAxis = 1
    else:
        raise ValueError('Unexpected permutation strategy', permStrategy)

    signCellsName, signCellsMouseDict = list(signCellsSelector.items())[0]

    nMice = len(dataDB.mice)

    figTot, axTot = plt.subplots(figsize=(4, 4))
    axTot.set_ylabel("AvgBinOrd")

    fig, ax = plt.subplots(nrows=2, ncols=nMice, figsize=(nMice * 4, 2 * 4))
    ax[0, 0].set_ylabel("AvgBinOrd")
    ax[1, 0].set_ylabel("Effect tize")

    phaseKeys = dataDB.get_phasetype_keys(phasetype, performance, haveWaiting=haveWaiting)

    print(phaseKeys)

    xDummy = np.arange(1, len(phaseKeys) + 1)

    fTrueAllLst = []
    fRandAllLst = []

    for iMouse, mousename in enumerate(sorted(dataDB.mice)):
        print("Doing mouse", mousename)
        queryDictThis = {'datatype': datatype, 'mousename': mousename, 'performance': performance}

        fTrueLst = []
        fRandLst = []
        effSizeLst = []

        for phaseKey in phaseKeys:
            # Calculate orderability and temporal means
            settings = {"dropShortTrialsTHR": 4}
            dataTempMu = metric_by_selector(
                dataDB,
                queryDictThis,
                "ord_mean",
                "rp",
                {phasetype: phaseKey},
                settings,
                channelFilter=signCellsMouseDict)

            func_ABO = sequence.avg_bivariate_binary_orderability_from_temporal_mean

            permSettings = {"haveEffectSize": True, "haveMeans": True}
            pValL, pValR, effSize, fTrue, fRand = perm_test(
                func_ABO, dataTempMu, permAxis=permAxis, iterAxis=iterAxis, nSample=200, settings=permSettings)

            fTrueLst += [fTrue]
            fRandLst += [fRand]
            effSizeLst += [effSize]

        fTrueAllLst += [fTrueLst]
        fRandAllLst += [fRandLst]

        ax[0, iMouse].set_title(mousename)
        ax[0, iMouse].plot(xDummy, fTrueLst)
        ax[0, iMouse].plot(xDummy, fRandLst, '--r')
        ax[1, iMouse].plot(xDummy, effSizeLst)

        #         ax[0, iMouse].axvline(x=5.5, linestyle='--', color='y')
        #         ax[0, iMouse].axvline(x=8.5, linestyle='--', color='y')
        #         ax[1, iMouse].axvline(x=5.5, linestyle='--', color='y')
        #         ax[1, iMouse].axvline(x=8.5, linestyle='--', color='y')

        ax[0, iMouse].set_ylim([0, None])
        ax[1, iMouse].set_ylim([0, None])
        ax[0, iMouse].set_xticks(xDummy)
        ax[1, iMouse].set_xticks(xDummy)
        ax[0, iMouse].set_xticklabels(phaseKeys)
        ax[1, iMouse].set_xticklabels(phaseKeys)

    muTrue = np.mean(fTrueAllLst, axis=0)
    muRand = np.mean(fRandAllLst, axis=0)
    stdTrue = np.std(fTrueAllLst, axis=0)
    stdRand = np.std(fRandAllLst, axis=0)

    if withSEM:
        stdTrue /= np.sqrt(nMice)
        stdRand /= np.sqrt(nMice)

    axTot.set_title('Allmice')
    axTot.plot(xDummy, muTrue, color='blue')
    axTot.fill_between(xDummy, muTrue - stdTrue, muTrue + stdTrue, color='blue', alpha=0.3)
    axTot.fill_between(xDummy, muRand - stdRand, muRand + stdRand, color='blue', alpha=0.3)
    axTot.plot(xDummy, muRand, '--r')

    #     axTot.axvline(x=5.5, linestyle='--', color='y')
    #     axTot.axvline(x=8.5, linestyle='--', color='y')
    axTot.set_ylim([0, None])
    axTot.set_xticks(xDummy)
    axTot.set_xticklabels(phaseKeys)

    plotSuffix = '_'.join([datatype, signCellsName, 'permstrat', permStrategy])

    plt.figure(figTot.number)
    plt.savefig("orderability_test_by_intervals_allmice_" + plotSuffix + ".pdf")

    plt.figure(fig.number)
    plt.savefig("orderability_test_by_intervals_" + plotSuffix + ".pdf")
    plt.show()
