import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

from IPython.display import display

from sklearn.decomposition import PCA
from scipy.ndimage import gaussian_filter1d

from mesostat.utils.signals import resample_stretch
from mesostat.stat.permtests import difference_test
from mesostat.utils.hdf5_io import DataStorage
from mesostat.utils.pandas_helper import get_rows_colvals


from src.lib.plots_lib import plot_coloured_1D

'''
    TODO:
    [+] Generic plot by interval
    [ ] Generic plot by session
    [ ] Trial-average plot (concat to minimum)
    [ ] Arbitrary PCA combination
'''

class PCAPlots:
    def __init__(self, dataDB, selector, queryDict, channelFilterMouse=None, h5name='pcaplotdata.h5'):
        self.dataDB = dataDB
        self.selector = selector
        self.queryDict = queryDict
        self.channelFilterMouse = channelFilterMouse
        self.h5name = h5name

        data3DAll = self.dataDB.get_data_from_selector(self.selector, queryDict)
        data2DAll = np.hstack(data3DAll)
        data2DAll = self._filter_channels(data2DAll)

        self.pca = PCA(n_components=2)
        self.pca.fit(data2DAll.T)

    def _filter_channels(self, data2D):
        if self.channelFilterMouse is not None:
            return data2D[self.channelFilterMouse]
        else:
            return data2D

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
                data = self.dataDB.get_data_from_selector(self.selector, queryDictThis)

                dataEff = self._preprocess(data, param)

                param['color'] = cmap(iColor)
                param['label'] = str([direction, performance])

                haveLabel = False
                for dataTrial in dataEff:
                    labelEff = None if haveLabel else param['label']
                    haveLabel = True
                    x, y = self._transform(self._filter_channels(dataTrial))
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
                data = self.dataDB.get_data_from_selector(self.selector, queryDictThis)

                data2DCond = np.array([np.mean(d, axis=1) for d in data]).T
                data2DCond = self._filter_channels(data2DCond)

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
                data = self.dataDB.get_data_from_selector(self.selector, queryDictThis)

                haveLabel = False
                label = str([direction, performance])
                for data2D in data:
                    x = self._transform(self._filter_channels(data2D))[iPCA]

                    labelEff = None if haveLabel else label
                    haveLabel = True
                    ax.plot(x, label=labelEff, alpha=0.5, color=cmap(iColor))

                iColor += 1

        ax.legend()
        ax.set_xlabel("timestep")
        ax.set_ylabel("PCA" + str(iPCA))

    def calc_distances(self, param):
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

        # cmaps = ['Reds_r', 'Greens_r', 'Blues_r']

        ds = DataStorage(self.h5name)

        for iTest, test in enumerate(tests):
            if test["condition"] == "direction":
                queryDictExtraA = {'performance': test["performance"], 'direction': 'L'}
                queryDictExtraB = {'performance': test["performance"], 'direction': 'R'}
            else:
                queryDictExtraA = {'direction': test["direction"], 'performance': 'Correct'}
                queryDictExtraB = {'direction': test["direction"], 'performance': 'Mistake'}

            dataA = self.dataDB.get_data_from_selector(self.selector, {**self.queryDict, **queryDictExtraA})
            dataB = self.dataDB.get_data_from_selector(self.selector, {**self.queryDict, **queryDictExtraB})

            dataA = [self._filter_channels(data2D) for data2D in dataA]
            dataB = [self._filter_channels(data2D) for data2D in dataB]

            dataAEff = np.array(self._preprocess(dataA, param)).transpose((2, 0, 1))
            dataBEff = np.array(self._preprocess(dataB, param)).transpose((2, 0, 1))

            distPval = [difference_test(dist_func, a, b, 1000, sampleFunction="permutation")[1] for a, b in zip(dataAEff, dataBEff)]
            distNegLogPval = np.array([-np.log10(p) for p in distPval])

            attrDict = {**self.selector, **self.queryDict, **test, **param}
            attrDict = {k : str(v) for k,v in attrDict.items()}
            ds.save_data('pca_dist', distNegLogPval, attrDict)

    # def plot_distances(self, ax, param):
    #     '''
    #     Plan:
    #         * Align data (crop or stretch)
    #         * Consider 3 binary tests (Dir(Corr), Perf(L), Perf(R))
    #         * Calc all-cell trial-average euclidean distance as function if timestep, plot
    #         * Do permutation test over labels, get pval for each timestep, use a colormap for above
    #     '''
    #
    #     # Avg over trials -> Norm Dist over neurons
    #     dist_func = lambda a, b: np.linalg.norm(np.mean(a, axis=0) - np.mean(b, axis=0))  #/ a.shape[1]
    #
    #     tests = [
    #         {"condition": "direction", "performance" : "Correct"},
    #         {"condition": "performance", "direction": "L"},
    #         {"condition": "performance", "direction": "R"}
    #     ]
    #
    #     cmaps = ['Reds_r', 'Greens_r', 'Blues_r']
    #
    #     for iTest, test in enumerate(tests):
    #         if test["condition"] == "direction":
    #             queryDictExtraA = {'performance': test["performance"], 'direction': 'L'}
    #             queryDictExtraB = {'performance': test["performance"], 'direction': 'R'}
    #         else:
    #             queryDictExtraA = {'direction': test["direction"], 'performance': 'Correct'}
    #             queryDictExtraB = {'direction': test["direction"], 'performance': 'Mistake'}
    #
    #         dataA = self.dataDB.get_data_from_selector(self.selector, {**self.queryDict, **queryDictExtraA})
    #         dataB = self.dataDB.get_data_from_selector(self.selector, {**self.queryDict, **queryDictExtraB})
    #
    #         dataA = [self._filter_channels(data2D) for data2D in dataA]
    #         dataB = [self._filter_channels(data2D) for data2D in dataB]
    #
    #         dataAEff = np.array(self._preprocess(dataA, param)).transpose((2, 0, 1))
    #         dataBEff = np.array(self._preprocess(dataB, param)).transpose((2, 0, 1))
    #
    #         # distTrue = np.array([dist_func(a, b) for a, b in zip(dataAEff, dataBEff)])
    #         distPval = [difference_test(dist_func, a, b, 1000, sampleFunction="permutation")[1] for a, b in zip(dataAEff, dataBEff)]
    #         distNegLogPval = np.array([-np.log10(p) for p in distPval])
    #
    #         x = np.arange(len(distNegLogPval))
    #         label = str("-".join(list(test.values())))
    #         ax.plot(x, distNegLogPval, label=label)
    #
    #         # x = np.arange(len(distTrue))
    #         # plot_coloured_1D(ax, x, distTrue + iTest, distNegLogPval, cmap=cmaps[iTest], vmin=0, vmax=4, haveColorBar=(iTest==0))
    #         # ax.text(len(distTrue) / 2, iTest, label, ha='left', wrap=True)
    #
    #     ax.set_ylim(0, 3.5)
    #     ax.axhline(y=2, linestyle='--', color='r')
    #     ax.legend()


def piecewise_scaling(lens, scales):
    # Stretch each piece according to its scale
    points = [np.linspace(0, s, l) for l,s in zip(lens, scales)]

    # Shift each piece, so that it goes after previous one
    for i in range(1, len(lens)):
        points[i] += 2*points[i-1][-1] - points[i-1][-2]

    xCat = np.hstack(points)
    xEnds = np.array([p[-1] for p in points])

    return xCat, xEnds


def plot_distances(h5name, dataDB):
    ds = DataStorage(h5name)
    dfMetaData = ds.list_dsets_pd()

    dfMetaData = dfMetaData[dfMetaData['name'] == 'pca_dist']
    dfMetaData = dfMetaData.drop(['name', 'shape', 'datetime'], axis=1)

    dfMetaDataIter = dfMetaData.copy()
    dfMetaDataIter = dfMetaDataIter.drop(['dset', 'mousename', 'semiphase', 'condition', 'direction', 'performance'], axis=1)
    dfMetaDataIter = dfMetaDataIter.drop_duplicates()

    for idx, row in dfMetaDataIter.iterrows():
        dfFiltered = get_rows_colvals(dfMetaData, dict(row), dropQuery=True)

        mice = sorted(set(dfFiltered['mousename']))
        conditions = sorted(set(dfFiltered['condition']))
        semiphases = sorted(set(dfFiltered['semiphase']))

        plt.rcParams.update({'font.size': 12})
        fig, ax = plt.subplots(ncols=len(mice), figsize=(6*len(mice), 6))

        for iMouse, mousename in enumerate(mice):
            # Calculate relative durations of each semiphase
            semiRelLen = {}
            for semiphase in semiphases:
                dataTMPLst = dataDB.get_data_from_semiphase(semiphase, {'mousename': mousename})
                semiRelLen[semiphase] = np.mean([d.shape[1] for d in dataTMPLst])

            for condition in conditions:
                dfFiltered2 = get_rows_colvals(dfFiltered, {'mousename' : mousename, 'condition' : condition}, dropQuery=True)

                keyCond = 'direction' if condition == 'performance' else 'performance'
                keyCondVals = set(dfFiltered2[keyCond])

                for condVal in keyCondVals:
                    dfFiltered3 = dfFiltered2[dfFiltered2[keyCond] == condVal]

                    dataLstThis = []
                    scalesThis = []

                    for idx2, row2 in dfFiltered3.sort_values(by='semiphase').iterrows():
                        dataLstThis += [ds.get_data(row2['dset'])]
                        scalesThis += [semiRelLen[row2['semiphase']]]

                    lensThis = [len(d) for d in dataLstThis]

                    xThis, xEnds = piecewise_scaling(lensThis, scalesThis)
                    dataArrThis = np.concatenate(dataLstThis)
                    ax[iMouse].plot(xThis, dataArrThis, label=condition+'_'+condVal)

                    for x in xEnds[:-1]:
                        ax[iMouse].axvline(x=x, linestyle='--', color='gray')


            ax[iMouse].set_ylim(0, 3.5)
            ax[iMouse].set_title(mousename)
            ax[iMouse].axhline(y=2, linestyle='--', color='r')
            ax[iMouse].legend()

        fnameSuffix = '_'.join([k+'_'+v for k,v in dict(row).items()])
        plt.savefig("manifold_dist_vs_time_" + fnameSuffix + ".pdf", dpi=600)
        plt.close()
