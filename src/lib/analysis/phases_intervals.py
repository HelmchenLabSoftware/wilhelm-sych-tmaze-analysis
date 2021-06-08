import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import mannwhitneyu, wilcoxon, binom_test, combine_pvalues

from mesostat.utils.pandas_helper import pd_query
from mesostat.visualization.mpl_matrix import imshow
from mesostat.visualization.mpl_violin import violins_labeled
from mesostat.visualization.mpl_barplot import barplot_labeled
from mesostat.visualization.mpl_stat_annot import stat_annot_patches

from src.lib.metric_wrapper import metric_by_selector_all, metric_by_selector
import src.lib.tests_phase_signle_cell as single_cell_tests


def plot_avg_firing_rate_by_neuron(dataDB, datatype, phaseType, haveWaiting=True, cmap=None):
    settings = {"zscoreChannel": False, "serial": True, "metricSettings": {}}

    for mousename in sorted(dataDB.mice):
        print(mousename)

        fig, ax = plt.subplots(ncols=4, figsize=(4 * 2, 6), tight_layout=True)
        fig.suptitle(mousename)
        iCol = 0

        for performance in dataDB.get_performances():
            if phaseType != 'phase':
                boundingLines = dataDB.get_phase_bounding_lines(phaseType, performance, haveWaiting=haveWaiting)

            for direction in dataDB.get_directions():
                queryDict = {'datatype': datatype, 'mousename': mousename, 'performance': performance,
                             'direction': direction}
                means = metric_by_selector_all(dataDB, queryDict, phaseType, 'mean', 'p', settings,
                                               haveWaiting=haveWaiting)
                means = means.T

                # Sort by max argument
                # But only based on the first one (CL) so that all 4 combinations have the same sorting
                if iCol == 0:
                    idxsMaxArgs = np.argsort(np.array([np.argmax(m) for m in means]))
                means = means[idxsMaxArgs]

                title = direction[0] + performance[0]
                imshow(fig, ax[iCol], means, limits=None, title=title, haveColorBar=True, cmap=cmap)
                if phaseType != 'phase':
                    for bline in boundingLines:
                        ax[iCol].axvline(x=bline - 1, linestyle='--', color='r', alpha=0.5)

                iCol += 1

        plt.savefig('avgrate_' + datatype + '_' + mousename + '_' + phaseType + '.pdf', dpi=600)
        plt.close()


def plot_count_active_trials_by_neuron(dataDB, datatype, phaseType, thrAct=0.2, thrFreq=0.1, haveWaiting=True):
    nMice = len(dataDB.mice)
    for performance in dataDB.get_performances():
        fig, ax = plt.subplots(nrows=nMice, ncols=3, figsize=(12, 4*nMice), tight_layout=True)

        for iMouse, mousename in enumerate(sorted(dataDB.mice)):
            print(performance, mousename)

            freqActive = []
            phaseTypeKeys = dataDB.get_phasetype_keys(phaseType, performance, haveWaiting=haveWaiting)
            for phaseTypeKey in phaseTypeKeys:
                queryDict = {'datatype': datatype, 'mousename': mousename, 'performance': performance}
                dataThis = dataDB.get_data_from_selector({phaseType : phaseTypeKey}, queryDict)

                # Consider trial as active if for any timestep activity exceeded certain threshold
                dataActive = np.array([np.any(d > thrAct, axis=1) for d in dataThis])

                # Compute frequency over trials with which cell is active
                freqActive += [np.mean(dataActive, axis=0)]

            # Load average activity sorting
            fname = '_'.join(['sorted_cell_idxs', mousename, datatype, phaseType, 'Correct', 'L']) + '.txt'
            with open(fname, 'r') as txtFile:
                sortIdxs = txtFile.read().replace('\n', ' ')[1:-1].split(' ')
                sortIdxs = np.array([s for s in sortIdxs if s != '']).astype(int)

            freqActive = np.array(freqActive).T
            freqActive = freqActive[sortIdxs]
            isActive = (freqActive > thrFreq).astype(int)

            # Plot frequency distribution
            imshow(fig, ax[iMouse][0], freqActive, limits=[0, 1], haveColorBar=True, xTicks=phaseTypeKeys)
            imshow(fig, ax[iMouse][1], isActive, limits=[0, 1], haveColorBar=True, xTicks=phaseTypeKeys)

            # Plot histogram of max activities
            maxFreqActiveCDF = np.hstack(([0], np.max(freqActive, axis=1)))
            maxFreqActiveCDF.sort()
            pValCDF = np.arange(len(maxFreqActiveCDF))
            ax[iMouse][2].plot(maxFreqActiveCDF, pValCDF)

            ax[iMouse][0].set_ylabel(mousename)
            ax[iMouse][2].set_xlabel('Frequency active')
            ax[iMouse][2].set_ylabel('Cell index (sorted)')

        ax[0][0].set_title('Frequency over trials when neuron active')
        ax[0][1].set_title('Active vs Passive Cells')
        ax[0][2].set_title('Best frequency over all ' + phaseType + 's')
        plt.savefig('freq_active_' + datatype + '_' + phaseType + '_' + performance + '.pdf', dpi=600)
        plt.close()


def plot_activity_vs_active_frequency(dataDB, datatype, phaseType, mousename, performance, thrAct=0.2, haveWaiting=True):
    FPSdata = 10
    phaseTypeKeys = dataDB.get_phasetype_keys(phaseType, performance, haveWaiting=haveWaiting)
    nPhaseType = len(phaseTypeKeys)
    nCell = dataDB.get_nchannel(mousename, datatype)
    freqActive = np.zeros((nCell, nPhaseType))
    avgActivityAll = np.zeros((nCell, nPhaseType))
    avgActivityActive = np.zeros((nCell, nPhaseType))
    phaseDuration = []

    for iPhaseType, phaseTypeKey in enumerate(phaseTypeKeys):
        queryDict = {'datatype': datatype, 'mousename': mousename, 'performance': performance}
        dataThis = dataDB.get_data_from_selector({phaseType : phaseTypeKey}, queryDict)

        # Consider trial as active if for any timestep activity exceeded certain threshold
        isActive2D = np.array([np.any(d > thrAct, axis=1) for d in dataThis])
        avgAct2D = np.array([np.mean(d, axis=1) for d in dataThis])

        # Compute frequency over trials with which cell is active
        freqActive[:, iPhaseType] = np.mean(isActive2D, axis=0)

        avgActivityAll[:, iPhaseType] = np.mean(avgAct2D, axis=0)

        # avgActivityActive[:, iPhaseType] = np.mean(isActive2D * avgAct2D, axis=0)
        avgActivityActive[:, iPhaseType] = [np.mean(avgAct2D[iCell][isActive2D[iCell]]) for iCell in range(nCell)]
        # avgActivityActive[:, iPhaseType] = np.mean([np.mean(d, axis=1) for d in dataThis], axis=0)

        phaseDuration += [[d.shape[1] for d in dataThis]]

    phaseAvgDuration = np.mean(phaseDuration, axis=1)
    expFiringRate = freqActive / phaseAvgDuration * FPSdata

    mpl.rcParams['legend.fontsize'] = 10
    fig, ax = plt.subplots(ncols=6, figsize=(20, 4), tight_layout=True)

    violins_labeled(ax[0], avgActivityAll.T, phaseTypeKeys, phaseType, 'Average activity all cells')
    violins_labeled(ax[1], freqActive.T, phaseTypeKeys, phaseType, 'Frequency active')
    violins_labeled(ax[2], avgActivityActive.T, phaseTypeKeys, phaseType, 'Average activity active cells')
    violins_labeled(ax[3], phaseDuration, phaseTypeKeys, phaseType, 'Phase Duration',  violinScale='width')
    violins_labeled(ax[4], expFiringRate.T, phaseTypeKeys, phaseType, 'Expected Activity Rate, $s^{-1}$')

    for iPhase, phase in enumerate(phaseTypeKeys):
        ax[5].plot(freqActive[:, iPhase], avgActivityActive[:, iPhase], '.', label=phase)

    ax[5].set_xlabel('Frequency active')
    ax[5].set_ylabel('Average Activity')
    ax[5].legend()

    plt.show()


def plot_significant_firing_rate_by_neuron(dataDB, datatype, phaseType, cmapConfusion, confThr=0.01, haveWaiting=True):
    settings = {"zscoreChannel": False, "serial": True, "metricSettings": {}}
    colorBoundingLines = 'white'

    performanceValues = dataDB.get_performances()
    directionValues = dataDB.get_directions()

    directionValues += ['All']
    if phaseType == 'phase':
        performanceValues += ['All']

    nCol = len(performanceValues) * len(directionValues)
    cumulConfDict = {}
    cumulMeansSignDict = {}

    for mousename in sorted(dataDB.mice):
        print('doing mouse', mousename)
        nChannel = dataDB.get_nchannel(mousename, datatype)

        fig1, ax1 = plt.subplots(ncols=nCol, figsize=(nCol * 2, 6), tight_layout=True)
        fig2, ax2 = plt.subplots(ncols=nCol, figsize=(nCol * 4, 4), tight_layout=True)

        fig1.suptitle(mousename)
        fig2.suptitle(mousename)

        if phaseType == 'phase':
            fig3, ax3 = plt.subplots(ncols=nCol, figsize=(nCol * 4, 4), tight_layout=True)
            fig3.suptitle(mousename)

        iCol = 0

        for performance in performanceValues:
            if phaseType != 'phase':
                boundingLines = dataDB.get_phase_bounding_lines(phaseType, performance, haveWaiting=haveWaiting)

            for direction in directionValues:
                queryDict = {'datatype': datatype, 'mousename': mousename}
                if performance != 'All':
                    queryDict = {**queryDict, **{'performance': performance}}
                if direction != 'All':
                    queryDict = {**queryDict, **{'direction': direction}}

                # Calculate p-values
                pVals2D = single_cell_tests.test_inverse_all_selectors(dataDB, queryDict, phaseType, metricName='mean',
                                                                       alternative="greater", haveWaiting=haveWaiting,
                                                                       settings=settings)
                negLogPVals2D = -np.log10(pVals2D)

                # Calculate significant cells and confusion matrix
                signCellsByPhase = single_cell_tests.pvalues_2_significant_cells(pVals2D, confThr)
                confMat = single_cell_tests.significance_confusion_matrix(signCellsByPhase) / nChannel

                # For each cell find phasetype where it has highest significance
                # Compute sorting index for cells by significant phasetype
                # But only based on the first one (CL) so that all 4 combinations have the same sorting
                if iCol == 0:
                    idxsMaxArgs = np.argsort(np.array([np.argmax(p) for p in negLogPVals2D]))

                    # Store cell incices
                    fname = '_'.join(['sorted_cell_idxs', mousename, datatype, phaseType, performance, direction]) + '.txt'
                    with open(fname, 'w') as txtFile:
                        txtFile.write(str(idxsMaxArgs))

                # Sort result by phasetype index
                negLogPVals2D = negLogPVals2D[idxsMaxArgs]

                ###########################
                # Plot significance for each cell and phaseType
                ###########################
                cmapSignificant = 'copper'  # 'YlOrBr' # 'viridis'
                title = direction[0] + performance[0]
                if phaseType != 'phase':
                    imshow(fig1, ax1[iCol], negLogPVals2D, limits=[0, 4], title=title, haveColorBar=True,
                           cmap=cmapSignificant)
                    for bline in boundingLines:
                        ax1[iCol].axvline(x=bline - 1, linestyle='--', color=colorBoundingLines, alpha=1)
                else:
                    imshow(fig1, ax1[iCol], negLogPVals2D, limits=[0, 4], title=title, haveTicks=True,
                           haveColorBar=False, cmap=cmapSignificant)

                ax1[iCol].set_yticks(np.arange(0, len(negLogPVals2D), 10))

                ###########################
                # Plot confusion matrix
                ###########################
                imshow(fig2, ax2[iCol], confMat, limits=[0, 1], title=title, haveColorBar=True, cmap=cmapConfusion)
                if phaseType != 'phase':
                    for bline in boundingLines:
                        ax2[iCol].axvline(x=bline - 1, linestyle='--', color=colorBoundingLines, alpha=1)
                        ax2[iCol].axhline(y=bline - 1, linestyle='--', color=colorBoundingLines, alpha=1)

                # Store confusion matrices
                cumulConfDict[(mousename, direction, performance)] = confMat

                #####################################
                # Compare rate/DFF for cells of different phases
                #####################################
                if phaseType == 'phase':
                    meansByPhase = metric_by_selector_all(dataDB, queryDict, phaseType, 'mean', 'p', settings,
                                                          haveWaiting=haveWaiting)
                    meansSignByPhase = [means[list(idxs)] for means, idxs in zip(meansByPhase, signCellsByPhase)]

                    phaseNames = dataDB.get_phasetype_keys(phaseType, 'Correct', haveWaiting=haveWaiting)

                    sigTestPairs = []
                    nSignCells = [len(s) for s in signCellsByPhase]
                    if (nSignCells[0] > 1) and (nSignCells[1] > 1):
                        sigTestPairs += [(0, 1)]
                    if (nSignCells[1] > 1) and (nSignCells[2] > 1):
                        sigTestPairs += [(1, 2)]
                    if len(sigTestPairs) == 0:
                        sigTestPairs = None


                    print(mousename, datatype, phaseType, direction, performance,
                          np.round([np.mean(x) for x in meansSignByPhase], 3),
                          np.round([np.std(x) for x in meansSignByPhase], 3),
                          np.round([np.std(x) / np.sqrt(len(x)) for x in meansSignByPhase], 3))

                    ax3[iCol].set_title(title)
                    violins_labeled(ax3[iCol], meansSignByPhase, phaseNames, phaseType, "mean",
                                    joinMeans=False, printLogP=False, sigTestPairs=sigTestPairs,
                                    violinInner='box')  # 'point')

                    # Combine means for all mice
                    cumulMeansSignDict[(mousename, direction, performance)] = meansSignByPhase

                iCol += 1

        fig1.savefig(mousename + '_significantrate_' + datatype + '_' + phaseType + '.pdf', dpi=600)
        plt.close()
        fig2.savefig(mousename + '_significantrate_' + datatype + '_' + phaseType + '_confusion.pdf', dpi=600)
        plt.close()

        if phaseType == 'phase':
            fig3.savefig(mousename + '_significantcell_' + datatype + '_' + phaseType + '_rate_comparison.pdf', dpi=600)
            plt.close()

    iCol = 0
    fig2All, ax2All = plt.subplots(ncols=nCol, figsize=(nCol * 4, 4), tight_layout=True)
    fig3All, ax3All = plt.subplots(ncols=nCol, figsize=(nCol * 4, 4), tight_layout=True)
    for performance in performanceValues:
        for direction in directionValues:
            title = direction[0] + performance[0]

            ###########################
            # Plot confusion matrix
            ###########################
            confThisCondition = [cumulConfDict[(mousename, direction, performance)] for mousename in dataDB.mice]
            confAvg = np.mean(confThisCondition, axis=0)

            imshow(fig2All, ax2All[iCol], confAvg, limits=[0, 1], title=title, haveColorBar=True, cmap=cmapConfusion)
            if phaseType != 'phase':
                boundingLines = dataDB.get_phase_bounding_lines(phaseType, performance, haveWaiting=haveWaiting)

                for bline in boundingLines:
                    ax2All[iCol].axvline(x=bline - 1, linestyle='--', color=colorBoundingLines, alpha=1)
                    ax2All[iCol].axhline(y=bline - 1, linestyle='--', color=colorBoundingLines, alpha=1)

            #####################################
            # Compare rate/DFF for cells of different phases
            #####################################
            if phaseType == 'phase':
                phaseNames = dataDB.get_phasetype_keys(phaseType, 'Correct', haveWaiting=haveWaiting)

                meansThisParam = [cumulMeansSignDict[(mousename, direction, performance)] for mousename in dataDB.mice]
                meansThisParam = [np.hstack([m[i] for m in meansThisParam]) for i in range(len(phaseNames))]

                print(datatype, phaseType, direction, performance,
                      np.round([np.mean(x) for x in meansThisParam], 3),
                      np.round([np.std(x) for x in meansThisParam], 3),
                      np.round([np.std(x) / np.sqrt(len(x)) for x in meansThisParam], 3))

                ax3All[iCol].set_title(title)
                violins_labeled(ax3All[iCol], meansThisParam, phaseNames, phaseType, "mean",
                                joinMeans=False, printLogP=False, sigTestPairs=sigTestPairs,
                                violinInner='box')  # 'point')

            iCol += 1

    fig2All.savefig('Allmice_significantrate_' + datatype + '_' + phaseType + '_confusion.pdf', dpi=600)
    plt.close()

    if phaseType == 'phase':
        fig3All.savefig('Allmice_significantcell_' + datatype + '_' + phaseType + '_rate_comparison.pdf', dpi=600)
        plt.close()


def plot_df_count_by_mice(ax, data1, data2, label1, label2, xLabels, nCellPerMouse):
    dataNorm1 = np.array(data1) / np.array(nCellPerMouse)
    dataNorm2 = np.array(data2) / np.array(nCellPerMouse)
    pVals = [binom_test(d1, d1 + d2) for d1, d2 in zip(data1, data2)]

    xInd = np.arange(len(xLabels))
    width = 0.35  # the width of the bars

    rects1 = ax.bar(xInd - width / 2, dataNorm1, width, label=label1)
    rects2 = ax.bar(xInd + width / 2, dataNorm2, width, label=label2)

    for patch1, patch2, pVal in zip(rects1.patches, rects2.patches, pVals):
        stat_annot_patches(ax, patch1, patch2, pVal, fontsize=20)

    ax.set_xticks(xInd)
    ax.set_xticklabels(xLabels)


def plot_df_count_combined(ax, data1, data2, label1, label2, nCellPerMouse):
    # Convert from counts to fractions
    dataNorm1 = np.array(data1) / np.array(nCellPerMouse)
    dataNorm2 = np.array(data2) / np.array(nCellPerMouse)
    dataMean1 = np.mean(dataNorm1)
    dataMean2 = np.mean(dataNorm2)
    dataStdMean1 = np.std(dataNorm1) / np.sqrt(len(dataNorm1))
    dataStdMean2 = np.std(dataNorm2) / np.sqrt(len(dataNorm2))

    # Compute p-values and combine
    pVals = [binom_test(d1, d1 + d2) for d1, d2 in zip(data1, data2)]
    pValMean = combine_pvalues(pVals)[1]

    # Plot Bars
    width = 0.30  # the width of the bars
    dist = 0.05
    rects1 = ax.bar([-(width + dist) / 2], [dataMean1], width, yerr=dataStdMean1, label=label1)
    rects2 = ax.bar([+(width + dist) / 2], [dataMean2], width, yerr=dataStdMean2, label=label2)

    # Plot individual lines
    for d1, d2 in zip(dataNorm1, dataNorm2):
        ax.plot([-(width + dist) / 2, (width + dist) / 2], [d1, d2], color='gray')

    # Annotate bars
    for patch1, patch2 in zip(rects1.patches, rects2.patches):
        stat_annot_patches(ax, patch1, patch2, pValMean, fontsize=20)

    ax.set_xticks([])
    ax.set_ylim([0, 1.1 * np.max(np.hstack([dataNorm1, dataNorm2]))])


def plot_save_significantly_firing_neurons(dataDB, datatype, phaseType, exclusiveIndices, exclusiveLabels, ranges=None, confThr=0.01, haveAll=False):
    settings = {"zscoreChannel": False, "serial": True, "metricSettings": {}}

    exclusiveLabelsKey = '_'.join(exclusiveLabels)

    significantCellsDict = {label: [] for label in exclusiveLabels}

    if not haveAll:
        performanceValues = ['Correct', 'Mistake']
        directionValues = ['L', 'R']
    else:
        performanceValues = ['Correct', 'Mistake', 'All']
        directionValues = ['L', 'R', 'All']

    mice = list(sorted(dataDB.mice))
    for mousename in mice:
        print('doing mouse', mousename)

        for performance in performanceValues:
            for direction in directionValues:
                queryDict = {'datatype': datatype, 'mousename': mousename}
                if performance != 'All':
                    queryDict = {**queryDict, **{'performance': performance}}
                if direction != 'All':
                    queryDict = {**queryDict, **{'direction': direction}}

                # Calculate p-values and significant cells
                pVals2D = single_cell_tests.test_inverse_all_selectors(dataDB, queryDict, phaseType, metricName='mean',
                                                                       alternative="greater", settings=settings,
                                                                       ranges=ranges)
                signCellsByPhase = single_cell_tests.pvalues_2_significant_cells(pVals2D, confThr)

                # Determine maintenance-significant cells for storage
                exclusiveSets = single_cell_tests.find_exclusive_sets(signCellsByPhase, exclusiveIndices)

                for thisSet, thisLabel in zip(exclusiveSets, exclusiveLabels):
                    significantCellsDict[thisLabel] += [[mousename, performance, direction, np.array(list(thisSet))]]

    # Store significant cells in HDF5
    significantCellsDFDict = {}
    for label, signCellsData in significantCellsDict.items():
        dfCells = pd.DataFrame(signCellsData, columns=['mousename', 'performance', 'direction', 'cells'])
        dfCells.to_hdf('significant_cells_' + datatype + '_' + label + '.h5', key='df', mode='w')
        significantCellsDFDict[label] = dfCells

    nCellPerMouse = np.array([dataDB.get_nchannel(mousename, datatype) for mousename in mice])

    # Plot selector-specific number of significant cells vs mice
    for performance in performanceValues:
        for direction in directionValues:
            queryDict = {"performance": performance, "direction": direction}

            nSignCellByMouse = []
            for label, signCellsDF in significantCellsDFDict.items():
                rows = pd_query(signCellsDF, queryDict, dropQuery=True)
                nSignCellByMouse += [[len(list(row.values)[1]) for idx, row in rows.iterrows()]]

            ######################
            #  Plot by mouse
            ######################
            print('plotting', exclusiveIndices, performance, direction)

            figCount, axCount = plt.subplots(figsize=(5, 5))
            plot_df_count_by_mice(axCount, *nSignCellByMouse, *exclusiveLabels, mice, nCellPerMouse)

            figKey = '_'.join([exclusiveLabelsKey, datatype, phaseType, performance, direction])

            axCount.legend()
            figCount.savefig('significantrate_by_mouse_' + figKey + '.pdf')
            plt.close()

            ######################
            #  Plot combined
            ######################

            figCountAll, axCountAll = plt.subplots(figsize=(5, 5))
            plot_df_count_combined(axCountAll, *nSignCellByMouse, *exclusiveLabels, nCellPerMouse)
            print(figKey, exclusiveLabels, nSignCellByMouse)

            figCountAll.savefig('significantrate_combined_' + figKey + '.pdf')
            plt.close()


def plot_ratio_enc_mt(dataDB, datatype):
    settings = {"zscoreChannel": False, "serial": True, "metricSettings": {}}

    nPlots = len(dataDB.mice)
    fig, ax = plt.subplots(ncols=nPlots, figsize=(4 * nPlots, 4), tight_layout=True)

    for iMouse, mousename in enumerate(sorted(dataDB.mice)):
        queryDict = {'datatype': datatype, 'mousename': mousename}

        rezEnc = metric_by_selector(dataDB, queryDict, 'mean', 'p', {'phase': 'Encoding'}, settings, channelFilter=None)
        rezMt = metric_by_selector(dataDB, queryDict, 'mean', 'p', {'phase': 'Maintenance'}, settings,
                                   channelFilter=None)

        rezRatio = (rezMt - rezEnc) / (rezMt + rezEnc)

        cdf = lambda x: (np.sort(x), np.linspace(0, 1, len(x)))

        cdfX, cdfY = cdf(rezRatio)
        idx0 = np.argmin(cdfX ** 2)
        yMin = cdfY[idx0]

        ax[iMouse].plot(cdfX, cdfY)
        ax[iMouse].set_title(mousename)
        ax[iMouse].axvline(x=0, linestyle='--', color='r')
        ax[iMouse].axhline(y=yMin, linestyle='--', color='r')

    plt.show()


def plot_violins_by_phase(dataDB, datatype, phaseType, metricName, settings, haveWaiting=True, signCellsSelector=None):
    if signCellsSelector == None:
        signCellsSelector = {'None': None}

    signCellsKey, signCellsMouseDict = list(signCellsSelector.items())[0]

    for performance in ["Correct", "Mistake"]:
        fig, ax = plt.subplots(ncols=3, figsize=(15, 5))
        fig.suptitle(performance)

        intervalOrigIndices = dataDB.get_phasetype_keys(phaseType, performance)
        boundingLines = dataDB.get_phase_bounding_lines(phaseType, performance, haveWaiting=haveWaiting) - 1
        phases = dataDB.get_phasetype_keys('phase', performance, haveWaiting=haveWaiting)

        for mousename in sorted(dataDB.mice)[:-1]:
            queryThis = {"datatype": datatype, "mousename": mousename, "performance": performance}  # "direction" : "R"
            rez2D = metric_by_selector_all(dataDB, queryThis, phaseType, metricName, 'p', settings,
                                           haveWaiting=haveWaiting, channelFilter=signCellsMouseDict)
            #         rez2D = rez2D[..., 0]
            barplot_labeled(ax[0], rez2D, intervalOrigIndices, plotLabel=mousename, alpha=0.2, vlines=boundingLines)

        rez2DIntervAll = metric_by_selector_all(dataDB, queryThis, phaseType, metricName, 'p', settings,
                                                haveWaiting=haveWaiting, channelFilter=signCellsMouseDict)
        rez2DPhaseAll = metric_by_selector_all(dataDB, queryThis, 'phase', metricName, 'p', settings,
                                               haveWaiting=haveWaiting, channelFilter=signCellsMouseDict)
        barplot_labeled(ax[1], rez2DIntervAll, intervalOrigIndices, vlines=boundingLines)
        violins_labeled(ax[2], rez2DPhaseAll, phases, "phase", metricName, joinMeans=True, printLogP=True,
                        violinInner='box', sigTestPairs=[[0, 1], [1, 2]])

        ax[0].legend()
        plt.savefig(phaseType + '_' + datatype + '_' + signCellsKey + '_avg_violins_' + performance + '.pdf')
        plt.show()
