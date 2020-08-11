import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

from sklearn.decomposition import PCA
from scipy.ndimage import gaussian_filter1d

from mesostat.utils.signals import resample_stretch
from mesostat.stat.permtests import difference_test

from src.lib.plots_lib import plot_coloured_1D

'''
    TODO:
    [+] Generic plot by interval
    [ ] Generic plot by session
    [ ] Trial-average plot (concat to minimum)
    [ ] Arbitrary PCA combination
'''

class PCAPlots:
    def __init__(self, dataDB, selector, queryDict, channelFilterMouse=None):
        self.dataDB = dataDB
        self.selector = selector
        self.queryDict = queryDict

        data3DAll = self._get_data(queryDict)
        data2DAll = np.hstack(data3DAll)

        if channelFilterMouse is not None:
            data2DAll = data2DAll[channelFilterMouse]

        print('yolo', data2DAll.shape)

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
            return gaussian_filter1d(data2D, 3, axis=1, mode="reflect")
            # convFunc = lambda m: np.convolve(m, filt, mode='full')
            # return np.apply_along_axis(convFunc, axis=0, arr=data2D)
        else:
            return data2D

    def _crop(self, data3D, strat):
        if strat == "cropmin":
            nTimeMin = np.min([data.shape[1] for data in data3D])
            return np.array([data[:, :nTimeMin] for data in data3D])
        elif strat == "stretch":
            return np.array([[resample_stretch(d, self.nTimeStretch) for d in dataTrial] for dataTrial in data3D])
        else:
            return data3D

    def _transform(self, data2D):
        return self.pca.transform(data2D.T).T

    def _preprocess(self, data3D, param):
        dataEff = deepcopy(data3D)

        # Apply cropping if requested
        dataEff = self._crop(dataEff, param['cropStrategy'])

        # Apply concatenation if requested
        dataEff = self._process_trials(dataEff, param['trialStrategy'])

        # Apply accumulation if requested
        dataEff = [self._accumulate(dataTrial, param['accStrategy']) for dataTrial in dataEff]

        return dataEff

    # This is required to have comparable number of timesteps between different behavioral conditions for distance calculation
    def set_stretch_timesteps(self, n):
        self.nTimeStretch = n

    def plot_pca(self, ax, param):
        iColor = 0
        cmap = plt.get_cmap("tab10")
        for performance in ['Correct', 'Mistake']:
            for direction in ['L', 'R']:
                queryDictThis = {**self.queryDict, **{'performance' : performance, 'direction' : direction}}
                data = self._get_data(queryDictThis)

                dataEff = self._preprocess(data, param)

                param['color'] = cmap(iColor)
                param['label'] = str([direction, performance])

                haveLabel = False
                for dataTrial in dataEff:
                    labelEff = None if haveLabel else param['label']
                    haveLabel = True
                    x, y = self._transform(dataTrial)
                    ax.plot(x, y, color=param['color'], label=labelEff, alpha=0.5)

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

    def plot_distances(self, ax, param):
        '''
        Plan:
            * Align data (crop or stretch)
            * Consider 3 binary tests (Dir(Corr), Perf(L), Perf(R))
            * Calc all-cell trial-average euclidean distance as function if timestep, plot
            * Do permutation test over labels, get pval for each timestep, use a colormap for above
        '''

        # Avg over trials -> Norm Dist over neurons
        dist_func = lambda a, b: np.linalg.norm(np.mean(a, axis=0) - np.mean(b, axis=0))  #/ a.shape[1]

        tests = [
            {"condition": "direction", "performance" : "Correct"},
            {"condition": "performance", "direction": "L"},
            {"condition": "performance", "direction": "R"}
        ]

        for iTest, test in enumerate(tests):
            if test["condition"] == "direction":
                queryDictExtraA = {'performance': test["performance"], 'direction': 'L'}
                queryDictExtraB = {'performance': test["performance"], 'direction': 'R'}
            else:
                queryDictExtraA = {'direction': test["direction"], 'performance': 'Correct'}
                queryDictExtraB = {'direction': test["direction"], 'performance': 'Mistake'}

            dataA = self._get_data({**self.queryDict, **queryDictExtraA})
            dataB = self._get_data({**self.queryDict, **queryDictExtraB})

            dataAEff = np.array(self._preprocess(dataA, param)).transpose((2, 0, 1))
            dataBEff = np.array(self._preprocess(dataB, param)).transpose((2, 0, 1))

            distTrue = np.array([dist_func(a, b) for a, b in zip(dataAEff, dataBEff)])
            distPval = [difference_test(dist_func, a, b, 1000, sampleFunction="permutation")[1] for a, b in zip(dataAEff, dataBEff)]
            distNegLogPval = np.array([-np.log10(p) for p in distPval])

            label = str("-".join(list(test.values())))
            # ax.plot(distTrue + iTest, label=label)

            x = np.arange(len(distTrue))
            plot_coloured_1D(ax, x, distTrue + iTest, distNegLogPval, vmin=0, vmax=4, haveColorBar=(iTest==0))
            ax.text(len(distTrue) / 2, iTest, label, ha='left', wrap=True)



