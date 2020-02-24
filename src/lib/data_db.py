import json
import numpy as np
import pandas as pd

from os.path import basename, dirname, join, isfile, splitext

from src.lib.sys_lib import strlst2date
# from codes.lib.data_io.yaro.mouse_performance import mouse_performance_allsessions
from src.lib.pandas_lib import filter_rows_colval, filter_rows_colvals, get_one_row
from src.lib.matlab_lib import loadmat
from src.lib.os_lib import getfiles_walk

from IPython.display import display
from ipywidgets import IntProgress

class BehaviouralNeuronalDatabase :
    def __init__(self, param):

        # Find and parse Data filenames
        self.mice = set()
        self.metaDataFrames = {}

        ##################################
        # Read Phase Description
        ##################################
        with open('behaviour_phases.json', 'r') as f:
            tmp = json.load(f)
            self.metaDataFrames['interval_maps'] = {
                "Correct" : pd.DataFrame(tmp['Correct'], columns=['index', 'phase', 'label']),
                "Mistake" : pd.DataFrame(tmp['Mistake'], columns=['index', 'phase', 'label'])
            }

        ##################################
        # Find and parse data files
        ##################################
        self._find_parse_data_files(param["root_path_data"])


    def _phase_to_interval(self, phase, performance):
        thisMap = self.metaDataFrames['interval_maps'][performance]
        thisPhaseMap = thisMap[thisMap['phase'] == phase]
        return min(thisPhaseMap.index), max(thisPhaseMap.index) + 1   # Note upper bound non-inclusive


    def _find_parse_data_files(self, rootPath):
        dataWalk = getfiles_walk(rootPath, [".mat", "AcceptedCells"])
        behavWalk = getfiles_walk(rootPath, [".mat", "Behavior"])

        drop_non_digit = lambda s: ''.join([i for i in s if i.isdigit()])
        
        def digits2date(s):
            sd = drop_non_digit(s)
            return strlst2date([sd[:4], sd[4:6], sd[6:]])


        dataSplit = [splitext(name)[0].split('_') for path, name in dataWalk]
        behavSplit = [splitext(name)[0].split('_') for path, name in behavWalk]

        dataDict = {
            "mousename" : [sp[0].lower() for sp in dataSplit],
            "date"      : [digits2date(drop_non_digit(sp[-1])) for sp in dataSplit],
            "mousekey"  : ['_'.join([sp[0].lower(), drop_non_digit(sp[-1])]) for sp in dataSplit],
            "path"      : [join(path, name) for path, name in dataWalk]
        }

        behavDict = {
            "mousename" : [sp[1].lower() for sp in behavSplit],
            "date"      : [digits2date(drop_non_digit(sp[-1])) for sp in behavSplit],
            "mousekey"  : ['_'.join([sp[1].lower(), drop_non_digit(sp[-1])]) for sp in behavSplit],
            "path"      : [join(path, name) for path, name in behavWalk]
        }

        self.metaDataFrames['neuro'] = pd.DataFrame(dataDict)
        self.metaDataFrames['behavior'] = pd.DataFrame(behavDict)

        self.mice = set(dataDict['mousename'])
        self.mice.update(behavDict['mousename'])


    def read_neuro_files(self):
        if 'neuro' in self.metaDataFrames.keys():
            nNeuroFiles = self.metaDataFrames['neuro'].shape[0]

            self.dataNeuronal = []
            progBar = IntProgress(min=0, max=nNeuroFiles, description='Read Neuro Data:')
            display(progBar)  # display the bar
            for idx, filepath in enumerate(self.metaDataFrames['neuro']['path']):
                self.dataNeuronal += [list(loadmat(filepath, waitRetry=3).values())[0]]
                progBar.value += 1

        else:
            print("No Neuro files loaded, skipping reading part")


    def read_behavior_files(self):
        FPS_RAW = 20.0

        if 'behavior' in self.metaDataFrames.keys():
            self.behaviorStateFrames = []
            behaviorStateKeys = []

            nBehaviorFiles = self.metaDataFrames['behavior'].shape[0]

            progBar = IntProgress(min=0, max=nBehaviorFiles, description='Read Neuro Data:')
            display(progBar)  # display the bar
            for idx, row in self.metaDataFrames['behavior'].iterrows():
                print("Reading", row['path'])
                M = loadmat(row['path'], waitRetry=3)

                firstTTLIdx = np.where(M['allSignals'][:, 0] > 2)[0][0]
                timeAlign = np.round(firstTTLIdx / 100) / FPS_RAW

                for trialDirection in ['L', 'R']:
                    for perf in ['Correct', 'Mistake']:
                        key = '_'.join(['Trial', trialDirection+'Whole', perf])

                        if key not in M.keys():
                            print('No trials found for', key, 'skipping')
                        else:

                            timeIntervals = M[key]
                            if timeIntervals.ndim == 1:
                                timeIntervals = timeIntervals[None, :]

                            nTrial, nEvent = timeIntervals.shape

                            if nEvent != len(self.metaDataFrames['interval_maps'][perf]):
                                raise ValueError("Unexpected", nEvent)

                            timeIntervalsAligned = timeIntervals - timeAlign
                            timeIntervalsDownsampled = timeIntervalsAligned / 2

                            self.behaviorStateFrames += [np.round(timeIntervalsDownsampled * 10).astype(int)]
                            behaviorStateKeys += [[row["mousekey"], trialDirection, perf]]

                progBar.value += 1

            self.metaDataFrames['behaviorStates'] = pd.DataFrame(behaviorStateKeys, columns=["mousekey", "trialDirection", "performance"])

        else:
            print("No Neuro files loaded, skipping reading part")


    def get_interval_session_fromrow(self, idxNeuro, rowNeuro, idxBehavior, rowBehavior, startState, endState):
        dataNeuroThis = self.dataNeuronal[idxNeuro]
        framesArrThis = self.behaviorStateFrames[idxBehavior]

        startFrameIdxs = framesArrThis[:, startState]
        endFrameIdxs = framesArrThis[:, endState]

        print(list(zip(startFrameIdxs, endFrameIdxs)))

        return [dataNeuroThis[start:end] for start, end in zip(startFrameIdxs, endFrameIdxs)]


    # Return list of shape [nTrial, [nFrame, nCell]] for selected interval
    # The number of frames vary over trial, so numpy array cannot be used
    def get_interval_session(self, mousekey, startState, endState, direction, performance):
        idxNeuro, rowNeuro = get_one_row(filter_rows_colvals(self.metaDataFrames['neuro'], {"mousekey" : mousekey}))
        
        idxBehavior, rowBehavior = get_one_row(filter_rows_colvals(self.metaDataFrames['behaviorStates'], {
            "mousekey" : mousekey, "trialDirection" : direction, "performance" : performance}
        ))

        if (idxNeuro is None) or (idxBehavior is None):
            return None
        else:
            return self.get_interval_session_fromrow(idxNeuro, rowNeuro, idxBehavior, rowBehavior, startState, endState)


    def get_phase_session(self, mousekey, phase, direction, performance):
        startState, endState = self._phase_to_interval(phase, performance)
        return self.get_interval_session(mousekey, startState, endState, direction, performance)


    def get_interval_mouse(self, mousename, startState, endState, direction, performance):
        rez = []

        rowsNeuro = filter_rows_colval(self.metaDataFrames['neuro'], "mousename", mousename)
        for idxNeuro, rowNeuro in rowsNeuro.iterrows():
            thisMouseKey = rowNeuro["mousekey"]
            thisBehaviorKey = {"mousekey": thisMouseKey, "trialDirection": direction, "performance": performance}

            idxBehavior, rowBehavior = get_one_row(filter_rows_colvals(self.metaDataFrames['behaviorStates'], thisBehaviorKey))

            if idxBehavior is None:
                print("No behaviour found for", thisBehaviorKey, "; skipping")
            else:
                rez += self.get_interval_session_fromrow(idxNeuro, rowNeuro, idxBehavior, rowBehavior, startState, endState)

        return rez


    def get_phase_mouse(self, mousename, phase, direction, performance):
        startState, endState = self._phase_to_interval(phase, performance)
        return self.get_interval_mouse(mousename, startState, endState, direction, performance)


    def get_interval_all(self, startState, endState, direction, performance):
        rez = []
        for mousename in self.mice:
            rez += self.get_interval_mouse(mousename, startState, endState, direction, performance)
        return rez


    def get_phase_all(self, phase, direction, performance):
        startState, endState = self._phase_to_interval(phase, performance)
        return self.get_interval_all(startState, endState, direction, performance)