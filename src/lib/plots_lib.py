import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import interpolate
from scipy.stats import mannwhitneyu, wilcoxon
from IPython.display import display
from statannot import add_stat_annotation

from mesostat.utils.plotting import imshowAddFakeColorBar
from src.lib.stat_lib import test_rank_sum_nan_aware, test_signed_rank_nan_aware


# Plot y(x) with color given by z(x)
def plot_coloured_1D(ax, x, y, z, cmap='jet', vmin=None, vmax=None, haveColorBar=False):
    nP = len(x)

    vminEff = np.min(z) if vmin is None else vmin
    vmaxEff = np.max(z) if vmax is None else vmax

    zNorm = (z - vminEff) / (vmaxEff - vminEff)
    zNorm = np.clip(zNorm, 0, 1)

    cmapFunc = plt.get_cmap(cmap)
    for i in range(1, nP):
        # ax.plot(x[i - 1:i + 1], y[i - 1:i + 1], color=(zNorm[i], 0, 0))
        ax.plot(x[i - 1:i + 1], y[i - 1:i + 1], color=cmapFunc(zNorm[i]))

    if haveColorBar:
        imshowAddFakeColorBar(ax.figure, ax, cmap=cmap, vmin=vmin, vmax=vmax)


def plot_labeled_violins(ax, dataLst, dataLabels, paramName, metricName, joinMeans=True, haveLog=False, sigTestPairs=None, printLogP=False):
    data = np.concatenate(dataLst)
    nDataLst = [len(d) for d in dataLst]

    labels = np.repeat(dataLabels, nDataLst)
    df = pd.DataFrame({paramName: labels, metricName: data})

    if haveLog:
        ax.set_yscale("log")

    sns.violinplot(ax=ax, data=df, x=paramName, y=metricName, cut=0)
    if sigTestPairs is not None:
        labelPairs = [(dataLabels[i], dataLabels[j]) for i, j in sigTestPairs]
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
                logPValArr[iData][jData] = np.round(test_signed_rank_nan_aware(dataLst[iData], dataLst[jData])[0], 2)

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


def plot_stretched_intervals(ax, dataDB, queryDict, intervalStart, intervalEnd, nInterp=100):

    rezLst = []  # [nInterv, nDataPoint]
    timesLst = []
    intervals = list(range(intervalStart, intervalEnd+1))
    for iInterv in intervals:
        dataLst = dataDB.get_data_from_interval(iInterv, iInterv + 1, queryDict)
        rez = []
        times = []
        for data in dataLst:
            rez += [np.nanmean(data, axis=0)]
            times += [np.linspace(iInterv, iInterv+1, len(rez[-1]))]

        rezLst += [rez]
        timesLst += [times]


    nInterv = len(intervals)
    nTrial = len(timesLst[0])

    timesArr = [np.hstack([timesLst[iInterv][iTrial] for iInterv in range(nInterv)]) for iTrial in range(nTrial)]
    rezArr = [np.hstack([rezLst[iInterv][iTrial] for iInterv in range(nInterv)]) for iTrial in range(nTrial)]

    xRez = np.linspace(intervalStart, intervalEnd, nInterp)
    yLst = [interpolate.interp1d(x, y, kind="linear")(xRez) for x,y in zip(timesArr, rezArr)]
    muY = np.mean(yLst, axis=0)
    stdY = np.std(yLst, axis=0)

    ax.fill_between(xRez, muY-stdY, muY+stdY)
    ax.plot(xRez, muY, color='blue')


def clustering_imshow_overplot_lines(ax, clustering, haveH=True, haveV=True):
    nCluster = np.max(clustering)
    nodePerClusterCumul = [np.sum(clustering <= j + 1) for j in range(nCluster)]

    # Plot orderability metric matrix, and cluster separators with red lines
    for clusterSep in nodePerClusterCumul:
        if haveV:
            ax.axvline(x=clusterSep - 0.5, color='r', alpha=0.5)
        if haveH:
            ax.axhline(y=clusterSep - 0.5, color='r', alpha=0.5)