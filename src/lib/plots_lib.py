import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import interpolate
from scipy.stats import mannwhitneyu, wilcoxon
from IPython.display import display
from copy import deepcopy
from statannot import add_stat_annotation

from mesostat.utils.pandas_helper import outer_product_df
from mesostat.stat.machinelearning import binary_classifier

from src.lib.metric_wrapper import metric_by_sweep, metric_by_phase


def _test_rank_sum_nan_aware(data1, data2):
    data1nonan = data1[~np.isnan(data1)]
    data2nonan = data2[~np.isnan(data2)]

    if (len(data1nonan) == 0) or (len(data2nonan) == 0) or (len(set(data1nonan) | set(data2nonan)) < 2):
        print("NANI",
              len(data1nonan),
              len(data2nonan),
              len(set(data1nonan)),
              len(set(data2nonan)),
              data1nonan[0],
              data2nonan[0],
              )

        logPval = 0
    else:
        logPval = np.log10(mannwhitneyu(data1nonan, data2nonan)[1])
    return logPval, [len(data1nonan), len(data2nonan)]


def _test_signed_rank_nan_aware(data1, data2):
    data1flat = data1.flatten()
    data2flat = data2.flatten()

    nanIdx = np.isnan(data1flat) | np.isnan(data2flat)
    data1nonan = data1flat[~nanIdx]
    data2nonan = data2flat[~nanIdx]

    if (len(data1nonan) == 0) or (len(data2nonan) == 0) or (len(set(data1nonan) | set(data2nonan)) < 2):
        print("NANI",
              len(data1nonan),
              len(data2nonan),
              len(set(data1nonan)),
              len(set(data2nonan)),
              data1nonan[0],
              data2nonan[0],
              )

        logPval = 0
    else:
        logPval = np.log10(wilcoxon(data1nonan, data2nonan)[1])
    nNoNan = np.sum(~nanIdx)
    return logPval, nNoNan


def plot_labeled_violins(ax, dataLst, dataLabels, paramName, metricName, joinMeans=True, haveLog=False, sigTestPairs=None, printLogP=False):
    data = np.concatenate(dataLst)
    nDataLst = [len(d) for d in dataLst]

    labels = np.repeat(dataLabels, nDataLst)
    df = pd.DataFrame({paramName: labels, metricName: data})

    if haveLog:
        ax.set_yscale("log")

    sns.violinplot(ax=ax, data=df, x=paramName, y=metricName, cut=0)
    if sigTestPairs is not None:
        labelPairs = [(dataLabels[i], dataLabels[j]) for i,j in sigTestPairs]
        dataSizes = [(len(df[df[paramName] == l1]), len(df[df[paramName] == l2])) for l1, l2 in labelPairs]
        dataTuples = [(df[df[paramName] == l1][metricName], df[df[paramName] == l2][metricName]) for l1, l2 in labelPairs]

        print("For", labelPairs[0], "of data size", dataSizes[0], "rank-sum-test is", mannwhitneyu(*dataTuples[0], alternative='two-sided')[1])

        add_stat_annotation(ax, data=df, x=paramName, y=metricName, box_pairs=labelPairs, test='Mann-Whitney', loc='inside', verbose=0)  # , text_format='full'

    tickCoords = np.arange(len(dataLabels))
    # plt.setp(ax, xticks=tickCoords, xticklabels=dataLabels)
    if joinMeans:
        ax.plot(tickCoords, [np.nanmean(d) for d in dataLst], '*--', color='blue')

    if printLogP:
        nData = len(dataLst)
        logPValArr = np.zeros((nData, nData))
        for iData in range(nData):
            for jData in range(iData + 1, nData):
                logPValArr[iData][jData] = np.round(_test_signed_rank_nan_aware(dataLst[iData], dataLst[jData])[0], 2)

        logPValArr = logPValArr + logPValArr.T
        display(logPValArr)


def plot_labeled_bars(ax, data, dataLabels, plotLabel=None, alpha=None):
    dataEff = data if data.ndim == 2 else data[:, None]
    nDim, nData = dataEff.shape
    x = np.arange(nDim)
    y = np.nanmean(dataEff, axis=1)

    ax.bar(x, y, width=1, alpha=alpha, label=plotLabel)
    plt.setp(ax, xticks=x, xticklabels=dataLabels)


def plot_bar_bounds(ax, bounds):
    for vline in bounds:
        ax.axvline(x=vline-0.5, color='r', linestyle='--')


def plot_stretched_intervals(ax, dataDB, queryDict, idxStart, idxEnd, nInterp=100):

    rezLst = []  # [nInterv, nDataPoint]
    timesLst = []
    for iInterv in range(idxStart, idxEnd):
        dataLst = dataDB.get_data_from_interval(iInterv, iInterv + 1, queryDict)
        rez = []
        times = []
        for data in dataLst:
            rez += [np.nanmean(data, axis=0)]
            times += [np.linspace(iInterv, iInterv+1, len(rez[-1]))]

        rezLst += [rez]
        timesLst += [times]


    nInterv = idxEnd-idxStart
    nTrial = len(timesLst[0])

    timesArr = [np.hstack([timesLst[iInterv][iTrial] for iInterv in range(nInterv)]) for iTrial in range(nTrial)]
    rezArr = [np.hstack([rezLst[iInterv][iTrial] for iInterv in range(nInterv)]) for iTrial in range(nTrial)]

    xRez = np.linspace(idxStart, idxEnd, nInterp)
    yLst = [interpolate.interp1d(x, y, kind="linear")(xRez) for x,y in zip(timesArr, rezArr)]
    muY = np.mean(yLst, axis=0)
    stdY = np.std(yLst, axis=0)

    ax.fill_between(xRez, muY-stdY, muY+stdY)
    ax.plot(xRez, muY, color='blue')


def table_discriminate_behavior(dataDB, selector, condition, sweepDict, metricName, trgDimOrder="r", settings=None, multiplexKey=None):
    if settings == None:
        settings = dict()
    sweepDF = outer_product_df(sweepDict)
    condValues = list(set(dataDB.metaDataFrames['behaviorStates'][condition]))

    metricValuesDict = {}
    for i, condVal in enumerate(condValues):
        sweepDFThis = sweepDF.copy()
        sweepDFThis[condition] = condVal
        rez = metric_by_sweep(dataDB, sweepDFThis, metricName, trgDimOrder, selector, settings, multiplexKey=multiplexKey)

        # If metric is not scalar, just average over all values
        ndim = rez[0].ndim
        if ndim > 1:
            print("Warning: Using nonscalar metric of", ndim-1, "dimensions")
            axis = tuple(range(1, ndim))
            rez = [np.mean(r, axis=axis) for r in rez]

        metricValuesDict[condVal] = rez

    rezDF = sweepDF.copy()
    for condVal in condValues:
        rezDF[condVal] = 0.0

    nCond = len(condValues)
    for i in range(nCond):
        for j in range(i+1, nCond):
            rezDF["-logp" + str((i, j))] = 0.0
    rezDF["nTrial"] = ":)"

    for idx, row in sweepDF.iterrows():
        for condVal in condValues:
            rezDF.at[idx, condVal] = np.nanmean(metricValuesDict[condVal][idx])

        nNoNanLst = []
        for i in range(nCond):
            for j in range(i + 1, nCond):
                logp, nNoNan = _test_rank_sum_nan_aware(metricValuesDict[condValues[i]][idx], metricValuesDict[condValues[j]][idx])
                rezDF.at[idx, "-logp" + str((i, j))] = np.round(-logp, 2)
                nNoNanLst += [nNoNan]

        rezDF.at[idx, "nTrial"] = str(nNoNanLst)

    display(rezDF)


def table_discriminate_phases(dataDB, sweepDict, phases, metricName, trgDimOrder="r", settings=None, multiplexKey=None):
    if settings == None:
        settings = dict()
    sweepDF = outer_product_df(sweepDict)

    metricValuesDict = {}
    for phase in phases:
        rez = metric_by_sweep(dataDB, sweepDF, metricName, trgDimOrder, {"phase" : phase}, settings, multiplexKey=multiplexKey)

        # If metric is not scalar, just average over all values
        ndim = rez[0].ndim
        if ndim > 1:
            print("Warning: Using nonscalar metric of", ndim-1, "dimensions")
            axis = tuple(range(1, ndim))
            rez = [np.mean(r, axis=axis) for r in rez]

        metricValuesDict[phase] = rez

        # if ndim > len(trgDimOrder) + int(multiplexKey is not None):
        #     raise ValueError("Cannot test non-scalar functions")

    rezDF = sweepDF.copy()
    for phase in phases:
        rezDF[phase] = 0.0

    nPhase = len(phases)
    for i in range(nPhase):
        for j in range(i+1, nPhase):
            rezDF["-logp" + str((i, j))] = 0.0
    rezDF["nTrial"] = ":)"

    for idx, row in sweepDF.iterrows():
        for phase in phases:
            rezDF.at[idx, phase] = np.nanmean(metricValuesDict[phase][idx])

        nNoNanLst = []
        for i in range(nPhase):
            for j in range(i + 1, nPhase):
                logp, nNoNan = _test_signed_rank_nan_aware(metricValuesDict[phases[i]][idx], metricValuesDict[phases[j]][idx])
                rezDF.at[idx, "-logp" + str((i, j))] = np.round(-logp, 2)
                nNoNanLst += [nNoNan]

        rezDF.at[idx, "nTrial"] = str(nNoNanLst)

    display(rezDF)


def table_binary_classification(dataDB, phase, binaryDimension, metric, dimOrdTarget, queryDict, settings, multiplexKey=None):
    paramValues = set(dataDB.metaDataFrames['behaviorStates'][binaryDimension])

    nMetricLabels = []
    metricLabels = []
    metricValsFlat = []
    for paramVal in paramValues:
        queryDict[binaryDimension] = paramVal
        metricVals = metric_by_phase(dataDB, queryDict, metric, dimOrdTarget, settings, phases=[phase], multiplexKey=multiplexKey)[0]
        if metricVals[0].ndim > 1:
            metricValsFlat += [val.flatten() for val in metricVals]
        else:
            metricValsFlat += list(metricVals)
        metricLabels += [paramVal] * len(metricVals)
        nMetricLabels += [len(metricVals)]

    bc = binary_classifier(np.array(metricValsFlat), np.array(metricLabels), method="looc", balancing=True)
    bc['nLabels'] = tuple(nMetricLabels)

    return bc


def table_binary_classification_bymouse(dataDB, phase, binaryDimension, metric, dimOrdTarget, queryDict, settings):
    rezDF = pd.DataFrame(columns=["mousename", "nLabels", "acc_train", "acc_test", "acc_naive", "p-value"])

    for mousename in dataDB.mice:
        print("Doing mouse", mousename)
        queryDict["mousename"] = mousename

        bc = table_binary_classification(dataDB, phase, binaryDimension, metric, dimOrdTarget, queryDict, settings)
        bc['mousename'] = mousename
        rezDF = rezDF.append(bc, ignore_index=True)

    display(rezDF)


# def clustering_plots(matRescaled, metricByConn, clustering):
#     nTime = matRescaled.shape[2]
#     clusterSortIdxs = np.argsort(clustering)  # Node indices so that clusters appear consecutive
#     nCluster = np.max(clustering)
#     nodePerClusterCumul = [np.sum(clustering <= j + 1) for j in range(nCluster)]
#
#     # Plot orderability metric matrix, and cluster separators with red lines
#     fig2, ax2 = plt.subplots(ncols=2, figsize=(8, 4))
#     ax2[0].imshow(metricByConn[clusterSortIdxs][:, clusterSortIdxs])
#     for clusterSep in nodePerClusterCumul:
#         ax2[0].axvline(x=clusterSep - 0.5, color='r', alpha=0.5)
#         ax2[0].axhline(y=clusterSep - 0.5, color='r', alpha=0.5)
#
#     # Plot average rescaled activity for each cluster
#     for iCluster in range(nCluster):
#         mu = np.mean(matRescaled[clustering == iCluster], axis=(0, 1)) * 50 + iCluster
#         std = np.std(matRescaled[clustering == iCluster], axis=(0, 1)) * 50
#         ax2[1].fill_between(np.arange(nTime), mu - std, mu + std, alpha=0.3)
#         ax2[1].plot(mu, label=str(iCluster))


def clustering_plots_chill(ax, metricByConn, clustering):
    clusterSortIdxs = np.argsort(clustering)  # Node indices so that clusters appear consecutive
    nCluster = np.max(clustering)
    nodePerClusterCumul = [np.sum(clustering <= j + 1) for j in range(nCluster)]

    # Plot orderability metric matrix, and cluster separators with red lines
    ax.imshow(metricByConn[clusterSortIdxs][:, clusterSortIdxs])
    for clusterSep in nodePerClusterCumul:
        ax.axvline(x=clusterSep - 0.5, color='r', alpha=0.5)
        ax.axhline(y=clusterSep - 0.5, color='r', alpha=0.5)