import numpy as np
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA

'''
    TODO:
    [+] Generic plot by interval
    [ ] Generic plot by session
    [ ] Trial-average plot (concat to minimum)
    [ ] Arbitrary PCA combination
'''


def plot_pca(ax, dataDB, selector, queryDict, condition, plot_func):
    if "phase" in selector.keys():
        data_func = lambda queryDictThis: dataDB.get_data_from_phase(selector["phase"], queryDictThis)
    elif "interval" in selector.keys():
        data_func = lambda queryDictThis: dataDB.get_data_from_interval(selector["interval"], selector["interval"]+1, queryDictThis)
    else:
        raise ValueError("Unexpected selector", selector)

    data3DAll = data_func(queryDict)
    data2DAll = np.hstack(data3DAll)

    pca = PCA(n_components=2)
    pca.fit(data2DAll.T)

    cmap = plt.get_cmap("tab10")
    for i, condVal in enumerate(set(dataDB.metaDataFrames['behaviorStates'][condition])):
        color = cmap(i)
        dataMouseCond = data_func({**queryDict, **{condition: condVal}})
        plot_func(ax, pca, dataMouseCond, condVal, color)

    ax.legend()
    ax.set_xlabel("PCA1")
    ax.set_ylabel("PCA2")


def _plot_time_avg_scatter(ax, pca, dataMouseCond, condVal, color):
    data2DCond = np.array([np.mean(d, axis=1) for d in dataMouseCond]).T
    x, y = pca.transform(data2DCond.T).T
    ax.plot(x, y, '.', label=condVal, alpha=0.5, color=color)


def _plot_pca_time_trial_avg(ax, pca, dataMouseCond, condVal, color):
    nTimeMin = np.min([data.shape[1] for data in dataMouseCond])
    dataCropped3D = np.array([data[:, :nTimeMin] for data in dataMouseCond])
    dataCropped2D = np.mean(dataCropped3D, axis=0)
    x, y = pca.transform(dataCropped2D.T).T
    ax.plot(x, y, label=condVal, alpha=0.5, color=color)


def _plot_pca_concat(ax, pca, dataMouseCond, condVal, color):
    data2DCond = np.hstack(dataMouseCond)
    x, y = pca.transform(data2DCond.T).T
    ax.plot(x, y, label=condVal, alpha=0.5, color=color)


def _plot_pca_concat_cumul(ax, pca, dataMouseCond, condVal, color):
    data2DCond = np.hstack(dataMouseCond)

    nChannel, nTime = data2DCond.shape
    for iTime in range(1, nTime):
        data2DCond[:, iTime] += data2DCond[:, iTime - 1]

    x, y = pca.transform(data2DCond.T).T
    ax.plot(x, y, label=condVal, alpha=0.5, color=color)


def _plot_pca_bytrial(ax, pca, dataMouseCond, condVal, color):
    haveLabel = False
    for dataTrial in dataMouseCond:
        label = None if haveLabel else condVal
        haveLabel = True
        x, y = pca.transform(dataTrial.T).T
        ax.plot(x, y, color=color, label=label, alpha=0.5)


def _plot_pca_bytrial_cumul(ax, pca, dataMouseCond, condVal, color):
    haveLabel = False
    for dataTrial in dataMouseCond:
        nChannel, nTime = dataTrial.shape
        for iTime in range(1, nTime):
            dataTrial[:, iTime] += dataTrial[:, iTime - 1]

        label = None if haveLabel else condVal
        haveLabel = True
        x, y = pca.transform(dataTrial.T).T
        ax.plot(x, y, color=color, label=label, alpha=0.5)
