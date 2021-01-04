import numpy as np
from scipy import interpolate


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