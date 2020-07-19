import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

from sklearn.decomposition import PCA
from scipy.ndimage import gaussian_filter1d

'''
    TODO:
    [+] Generic plot by interval
    [ ] Generic plot by session
    [ ] Trial-average plot (concat to minimum)
    [ ] Arbitrary PCA combination
'''

class PCAPlots:
    def __init__(self, dataDB, selector, queryDict):
        self.dataDB = dataDB
        self.selector = selector
        self.queryDict = queryDict

        data3DAll = self._get_data(queryDict)
        data2DAll = np.hstack(data3DAll)

        self.pca = PCA(n_components=2)
        self.pca.fit(data2DAll.T)


    def _get_data(self, queryDict):
        if "phase" in self.selector.keys():
            return self.dataDB.get_data_from_phase(self.selector["phase"], queryDict)
        elif "interval" in self.selector.keys():
            interval = self.selector["interval"]
            return self.dataDB.get_data_from_interval(interval, interval+1, queryDict)
        else:
            raise ValueError("Unexpected selector", self.selector)

    def _process_trials(self, data3D, strat):
        if strat == 'concat':
            return np.array([np.hstack(data3D)])
        elif strat == 'avg':
            return np.array([np.mean(data3D, axis=0)])
        else:
            return data3D

    def _accumulate(self, data2D, strat):
        if strat == 'cumulative':
            return data2D.cumsum(axis=1)
        elif strat == 'gaussfilter':
            return gaussian_filter1d(data2D, 2, axis=1, mode="reflect")
            # convFunc = lambda m: np.convolve(m, filt, mode='full')
            # return np.apply_along_axis(convFunc, axis=0, arr=data2D)
        else:
            return data2D

    def _crop(self, data3D, strat):
        if strat == "cropmin":
            nTimeMin = np.min([data.shape[1] for data in data3D])
            return np.array([data[:, :nTimeMin] for data in data3D])
        else:
            return data3D

    def _transform(self, data2D):
        return self.pca.transform(data2D.T).T

    def _plot_pca(self, ax, data, param):
        dataEff = deepcopy(data)

        # Apply cropping if requested
        dataEff = self._crop(dataEff, param['cropStrategy'])

        # Apply concatenation if requested
        dataEff = self._process_trials(dataEff, param['trialStrategy'])

        haveLabel = False
        for dataTrial in dataEff:
            dataTrialEff = self._accumulate(dataTrial, param['accStrategy'])

            labelEff = None if haveLabel else param['label']
            haveLabel = True
            x, y = self._transform(dataTrialEff)
            ax.plot(x, y, color=param['color'], label=labelEff, alpha=0.5)

    def plot_pca(self, ax, param):
        iColor = 0
        cmap = plt.get_cmap("tab10")
        for performance in ['Correct', 'Mistake']:
            for direction in ['L', 'R']:
                queryDictThis = {**self.queryDict, **{'performance' : performance, 'direction' : direction}}
                data = self._get_data(queryDictThis)

                param['color'] = cmap(iColor)
                param['label'] = str([direction, performance])
                self._plot_pca(ax, data, param)

                iColor += 1

        ax.legend()
        ax.set_xlabel("PCA1")
        ax.set_ylabel("PCA2")

    def plot_time_avg_scatter(self, ax):
        iColor = 0
        cmap = plt.get_cmap("tab10")
        for performance in ['Correct', 'Mistake']:
            for direction in ['L', 'R']:
                queryDictThis = {**self.queryDict, **{'performance' : performance, 'direction' : direction}}
                data = self._get_data(queryDictThis)

                data2DCond = np.array([np.mean(d, axis=1) for d in data]).T
                x, y = self._transform(data2DCond)

                label = str([direction, performance])
                ax.plot(x, y, '.', label=label, color=cmap(iColor))

                iColor += 1

        ax.legend()
        ax.set_xlabel("PCA1")
        ax.set_ylabel("PCA2")

    def plot_pca_vs_time(self, ax, iPCA):
        iColor = 0
        cmap = plt.get_cmap("tab10")
        for performance in ['Correct', 'Mistake']:
            for direction in ['L', 'R']:
                queryDictThis = {**self.queryDict, **{'performance' : performance, 'direction' : direction}}
                data = self._get_data(queryDictThis)

                haveLabel = False
                label = str([direction, performance])
                for data2D in data:
                    x = self._transform(data2D)[iPCA]

                    labelEff = None if haveLabel else label
                    haveLabel = True
                    ax.plot(x, label=labelEff, alpha=0.5, color=cmap(iColor))

                iColor += 1

        ax.legend()
        ax.set_xlabel("timestep")
        ax.set_ylabel("PCA" + str(iPCA))
