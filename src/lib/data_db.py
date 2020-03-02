import json
import numpy as np
import pandas as pd
from copy import deepcopy
from collections import OrderedDict

from os.path import basename, dirname, join, isfile, splitext

from src.lib.sys_lib import strlst2date
import src.lib.pandas_lib as pandas_lib
from src.lib.matlab_lib import loadmat
from src.lib.os_lib import getfiles_walk
from src.lib.baseline_lib import zscore, crop_quantile

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
            self.phases = {
                "Correct" : OrderedDict.fromkeys(self.metaDataFrames['interval_maps']["Correct"]["phase"]),
                "Mistake" : OrderedDict.fromkeys(self.metaDataFrames['interval_maps']["Mistake"]["phase"])
            }

        ##################################
        # Find and parse data files
        ##################################
        self._find_parse_data_files(param["root_path_data"])


    def _phase_to_interval(self, phase, performance):
        thisMap = self.metaDataFrames['interval_maps'][performance]

        if phase is None:  # Equivalently all phases
            thisPhaseMap = thisMap
        else:
            thisPhaseMap = thisMap[thisMap['phase'] == phase]

        nStates = len(thisMap)
        idxStart = min(thisPhaseMap.index)

        # Note 1: upper bound is inclusive
        # Note 2: phase next to each state denotes interval starting at that state and ending at next state
        # Note 3: phase indicator at last state is effectively useless, because trials stops at that state
        idxEnd = np.min([nStates-1, max(thisPhaseMap.index) + 1])

        return idxStart, idxEnd


    def _extract_high_activity(self, data, p):
        return np.array([crop_quantile(data[:, i], p) for i in range(data.shape[1])]).T


    def _find_parse_data_files(self, rootPath):
        # Remove non-digits from a string
        drop_non_digit = lambda s: ''.join([i for i in s if i.isdigit()])

        # # Convert string "YYYYMMDD" into date object
        # def _digits2date(sd):
        #     return strlst2date([sd[:4], sd[4:6], sd[6:]])

        dataWalk = getfiles_walk(rootPath, [".mat", "AcceptedCells"])
        behavWalk = getfiles_walk(rootPath, [".mat", "Behavior"])

        # Fill neuro dict
        self.metaDataFrames['neuro'] = pd.DataFrame([], columns=["mousename", "date", "mousekey", "path"])
        for path, name in dataWalk:
            sp = splitext(name)[0].split('_')
            mousename = sp[0].lower()  # Convert mousename to lowercase
            datestr = drop_non_digit(sp[-1])
            mousekey = mousename + '_' + datestr
            pathname = join(path, name)
            self.metaDataFrames['neuro'] = pandas_lib.add_list_as_row(self.metaDataFrames['neuro'], [mousename, datestr, mousekey, pathname])

        # Fill behaviour dict
        self.metaDataFrames['behavior'] = pd.DataFrame([], columns=["mousename", "date", "mousekey", "path"])
        for path, name in behavWalk:
            sp = splitext(name)[0].split('_')
            mousename = sp[1].lower()  # Convert mousename to lowercase
            datestr = drop_non_digit(sp[-1])
            mousekey = mousename + '_' + datestr
            pathname = join(path, name)
            self.metaDataFrames['behavior'] = pandas_lib.add_list_as_row(self.metaDataFrames['behavior'], [mousename, datestr, mousekey, pathname])

        self.mice = set(self.metaDataFrames['neuro']['mousename'])
        self.mice.update(set(self.metaDataFrames['behavior']['mousename']))


    def read_neuro_files(self):
        if 'neuro' in self.metaDataFrames.keys():
            nNeuroFiles = self.metaDataFrames['neuro'].shape[0]

            self.dataNeuronal = {"raw" : [], "high" : []}
            progBar = IntProgress(min=0, max=nNeuroFiles, description='Read Neuro Data:')
            display(progBar)  # display the bar
            for idx, filepath in enumerate(self.metaDataFrames['neuro']['path']):
                tracesRaw = list(loadmat(filepath, waitRetry=3).values())[0]
                self.dataNeuronal["raw"] += [tracesRaw]
                self.dataNeuronal["high"] += [self._extract_high_activity(tracesRaw, p=0.95)]
                # self.dataNeuronal["zscore"] += zscore(tracesRaw, axis=0)
                #
                # self._preprocess_neuronal(tracesRaw)

                progBar.value += 1

        else:
            print("No Neuro files loaded, skipping reading part")


    def read_behavior_files(self):
        FPS_TTL = 2000.0         # Microelectronics sampling timescale
        FPS_RAW = 20.0           # Framerate of raw calcium signal
        FPS_DOWNSAMPLED = 10.0   # Framerate of downsampled calcium signal

        if 'behavior' in self.metaDataFrames.keys():
            self.metaDataFrames['behaviorStates'] = pd.DataFrame([], columns=["mousekey", "direction", "performance"])
            self.behaviorStateFrames = []

            nBehaviorFiles = self.metaDataFrames['behavior'].shape[0]
            progBar = IntProgress(min=0, max=nBehaviorFiles, description='Read Neuro Data:')
            display(progBar)  # display the bar
            for idx, row in self.metaDataFrames['behavior'].iterrows():
                print("Reading", row['path'])
                M = loadmat(row['path'], waitRetry=3)

                firstTTLIdx = np.where(M['allSignals'][:, 0] > 2)[0][0]   # Trials start with first pulse
                nFramesRaw = np.round(firstTTLIdx / FPS_TTL * FPS_RAW)
                timeAlign = nFramesRaw / FPS_RAW                          # Compute corresponding event in calcium timescale

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

                            # Align raw calcium time intervals to the first pulse
                            timeIntervalsAligned = timeIntervals - timeAlign

                            # Convert time intervals into frames.
                            # Due to downsampling, framerate is lower than original
                            idxsFrameIntervalsAligned = np.round(timeIntervalsAligned * FPS_DOWNSAMPLED).astype(int)

                            self.behaviorStateFrames += [idxsFrameIntervalsAligned]
                            self.metaDataFrames['behaviorStates'] = pandas_lib.add_list_as_row(
                                self.metaDataFrames['behaviorStates'], [row["mousekey"], trialDirection, perf])

                progBar.value += 1

        else:
            print("No Neuro files loaded, skipping reading part")


    def get_data_session_interval_fromindex(self, idxNeuro, idxBehavior, startState, endState, dataType='raw'):
        if (idxNeuro is None) or (idxBehavior is None):
            raise ValueError("Unexpected indices", idxNeuro, idxBehavior)

        dataNeuroThis = self.dataNeuronal[dataType][idxNeuro]
        framesArrThis = self.behaviorStateFrames[idxBehavior]

        startFrameIdxs = framesArrThis[:, startState]
        endFrameIdxs = framesArrThis[:, endState]

        return [dataNeuroThis[start:end].T for start, end in zip(startFrameIdxs, endFrameIdxs)]


    def get_data_from_interval(self, startState, endState, queryDict):
        rez = []

        dataType = queryDict["datatype"] if "datatype" in queryDict.keys() else "raw"

        queryDictNeuro = {k : v for k,v in queryDict.items() if k in ["mousename", "mousekey"]}
        queryDictBehav = {k: v for k, v in queryDict.items() if k not in ["mousename", "datatype"]}

        rowsNeuro = pandas_lib.get_rows_colvals(self.metaDataFrames['neuro'], queryDictNeuro)
        for idxNeuro, rowNeuro in rowsNeuro.iterrows():
            queryDictBehav["mousekey"] = rowNeuro["mousekey"]
            rowsBehavThis = pandas_lib.get_rows_colvals(self.metaDataFrames['behaviorStates'], queryDictBehav)

            if len(rowsBehavThis) == 0:
                print("No behaviour found for", queryDictBehav, "; skipping")
            else:
                for idxBehavior, rowBehavior in rowsBehavThis.iterrows():
                    rezThis = self.get_data_session_interval_fromindex(idxNeuro, idxBehavior, startState, endState, dataType)

                    if np.any([np.prod(d.shape) == 0 for d in rezThis]):
                        print("Warning: ", rowNeuro["mousekey"], "has zero interval", startState, endState)

                    rez += rezThis

        return rez


    def get_data_from_phase(self, phase, queryDict):
        if "performance" in queryDict.keys():
            startState, endState = self._phase_to_interval(phase, queryDict["performance"])
            return self.get_data_from_interval(startState, endState, queryDict)
        else:
            # If performance is not specified, merge trials for all performance measures
            rez = []
            queryDictCopy = deepcopy(queryDict)
            for perf in set(self.metaDataFrames['behaviorStates']['performance']):
                queryDictCopy["performance"] = perf
                rez += self.get_data_from_phase(phase, queryDictCopy)
            return rez
