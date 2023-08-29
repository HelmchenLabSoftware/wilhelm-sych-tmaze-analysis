import json
import numpy as np
import pandas as pd
from copy import deepcopy
from collections import OrderedDict

from os.path import join, splitext, dirname

from mesostat.utils.matlab_helper import loadmat
from mesostat.utils.system import getfiles_walk
import mesostat.utils.pandas_helper as pandas_helper
from mesostat.utils.signals.filter import zscore

from pfc_mem_proj.lib.baseline_lib import crop_quantile

from IPython.display import display
from ipywidgets import IntProgress

class BehaviouralNeuronalDatabase :
    def __init__(self, param):
        # Set control parameters
        self.verbose = True    # Whether to print warnings

        # Find and parse Data filenames
        self.mice = set()
        self.sessions = set()
        self.metaDataFrames = {}

        ##################################
        # Read Phase Description
        ##################################
        self.read_interval_maps()

        ##################################
        # Find and parse data files
        ##################################
        self._find_parse_data_files('dff', param["root_path_dff"])
        self._find_parse_data_files('deconv', param["root_path_deconv"])
        self._find_parse_behaviour_files('behavior', param["root_path_dff"])

    # def phase_to_interval(self, phase, performance):
    #     thisMap = self.metaDataFrames['interval_maps'][performance]
    #
    #     if phase is None:  # Equivalently all phases
    #         thisPhaseMap = thisMap
    #     else:
    #         thisPhaseMap = thisMap[thisMap['phase'] == phase]
    #
    #     nStates = len(thisMap)
    #     idxStart = min(thisPhaseMap.index)
    #
    #     # Note 1: upper bound is inclusive
    #     # Note 2: phase next to each state denotes interval starting at that state and ending at next state
    #     # Note 3: phase indicator at last state is effectively useless, because trials stops at that state
    #     idxEnd = np.min([nStates-1, max(thisPhaseMap.index) + 1])
    #
    #     return idxStart, idxEnd

    # Remove non-digits from a string
    def _drop_non_digit(self, s):
        return ''.join([i for i in s if i.isdigit()])

    def _extract_high_activity(self, data, p):
        return np.array([crop_quantile(data[:, i], p) for i in range(data.shape[1])]).T

    def _find_parse_data_files(self, metaKey, rootPath):
        dataWalk = getfiles_walk(rootPath, [".mat", "AcceptedCells"])

        # Fill neuro dict
        self.metaDataFrames[metaKey] = pd.DataFrame([], columns=["mousename", "date", "session", "path"])
        for path, name in dataWalk:
            sp = splitext(name)[0].split('_')
            mousename = sp[0].lower()  # Convert mousename to lowercase
            datestr = self._drop_non_digit(sp[-1])
            session = mousename + '_' + datestr
            pathname = join(path, name)
            self.metaDataFrames[metaKey] = pandas_helper.pd_append_row(self.metaDataFrames[metaKey], [mousename, datestr, session, pathname])

        self.mice.update(set(self.metaDataFrames[metaKey]['mousename']))
        self.sessions.update(set(self.metaDataFrames[metaKey]['session']))

    def _find_parse_behaviour_files(self, metaKey, rootPath):
        behavWalk = getfiles_walk(rootPath, [".mat", "Behavior"])

        # Fill behaviour dict
        self.metaDataFrames[metaKey] = pd.DataFrame([], columns=["mousename", "date", "session", "path"])
        for path, name in behavWalk:
            sp = splitext(name)[0].split('_')
            mousename = sp[1].lower()  # Convert mousename to lowercase
            datestr = self._drop_non_digit(sp[-1])
            session = mousename + '_' + datestr
            pathname = join(path, name)
            self.metaDataFrames[metaKey] = pandas_helper.pd_append_row(self.metaDataFrames[metaKey], [mousename, datestr, session, pathname])

    def read_interval_maps(self):
        libPath = dirname(__file__)
        srcPath = dirname(libPath)
        phasesFilePath = join(srcPath, 'behaviour_phases.json')

        with open(phasesFilePath, 'r') as f:
            jFile = json.load(f)
            self.intervalsDict = {
                "Correct": {idx: label for idx, label in jFile['Intervals']['Correct']},
                "Mistake": {idx: label for idx, label in jFile['Intervals']['Mistake']},
            }

            self.phasesDict = jFile['Phases']
            self.semiphasesDict = jFile['Semiphases']

    def read_neuro_files(self):
        if 'dff' in self.metaDataFrames.keys():
            nNeuroFiles = self.metaDataFrames['dff'].shape[0]

            self.dataNeuronal = {"raw" : [], "high" : [], "zscore" : []}
            progBar = IntProgress(min=0, max=nNeuroFiles, description='Read DFF Data:')
            display(progBar)  # display the bar
            for idx, filepath in enumerate(self.metaDataFrames['dff']['path']):
                tracesRaw = list(loadmat(filepath, waitRetry=3).values())[0]

                self.dataNeuronal["raw"] += [tracesRaw]
                self.dataNeuronal["high"] += [self._extract_high_activity(tracesRaw, p=0.95)]
                self.dataNeuronal["zscore"] += [zscore(tracesRaw, axis=0)]

                progBar.value += 1
        else:
            print("No DFF files loaded, skipping reading part")

        if 'deconv' in self.metaDataFrames.keys():
            nNeuroFiles = self.metaDataFrames['deconv'].shape[0]

            self.dataNeuronal["deconv"] = []
            progBar = IntProgress(min=0, max=nNeuroFiles, description='Read DECONV Data:')
            display(progBar)  # display the bar
            for idx, filepath in enumerate(self.metaDataFrames['deconv']['path']):
                self.dataNeuronal["deconv"] += [loadmat(filepath, waitRetry=3)['Y_predict'].T]
                progBar.value += 1
        else:
            print("No DECONV files loaded, skipping reading part")

    def read_behavior_files(self):
        FPS_TTL = 2000.0         # Microelectronics sampling timescale
        FPS_RAW = 20.0           # Framerate of raw calcium signal
        FPS_DOWNSAMPLED = 10.0   # Framerate of downsampled calcium signal

        if 'behavior' in self.metaDataFrames.keys():
            self.metaDataFrames['behaviorStates'] = pd.DataFrame([], columns=["session", "direction", "performance"])
            self.behaviorStateFrames = []

            nBehaviorFiles = self.metaDataFrames['behavior'].shape[0]
            progBar = IntProgress(min=0, max=nBehaviorFiles, description='Read Neuro Data:')
            display(progBar)  # display the bar
            for idx, row in self.metaDataFrames['behavior'].iterrows():
                # print("Reading", row['path'])
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

                            if nEvent != len(self.intervalsDict[perf]):
                                raise ValueError("Unexpected", nEvent)

                            # Align raw calcium time intervals to the first pulse
                            timeIntervalsAligned = timeIntervals - timeAlign

                            # Convert time intervals into frames.
                            # Due to downsampling, framerate is lower than original
                            idxsFrameIntervalsAligned = np.round(timeIntervalsAligned * FPS_DOWNSAMPLED).astype(int)

                            self.behaviorStateFrames += [idxsFrameIntervalsAligned]
                            self.metaDataFrames['behaviorStates'] = pandas_helper.pd_append_row(
                                self.metaDataFrames['behaviorStates'],
                                [row["session"], trialDirection, perf]
                            )

                progBar.value += 1

        else:
            print("No Neuro files loaded, skipping reading part")

    def get_nchannel(self, mousename, datatype):
        dataFrameKey = datatype if datatype=="deconv" else "dff"
        rows = self.get_rows(dataFrameKey, {"mousename" : mousename})
        firstIdx = rows.index[0]
        return self.dataNeuronal[datatype][firstIdx].shape[1]

    def get_nsession(self, mousename, datatype):
        dataFrameKey = datatype if datatype=="deconv" else "dff"
        if mousename is None:
            return len(self.metaDataFrames[dataFrameKey])
        else:
            return len(self.get_rows(dataFrameKey, {"mousename": mousename}))

    def get_directions(self):
        return list(sorted(set(self.metaDataFrames['behaviorStates']['direction'])))

    def get_performances(self):
        return list(sorted(set(self.metaDataFrames['behaviorStates']['performance'])))

    def get_phasetype_keys(self, phaseType, performance, haveWaiting=True):
        if haveWaiting:
            if phaseType == 'interval':
                # TODO: Note that there are 1 less intervals than interval flags. Make flags and intervals distinct
                return np.array(list(self.intervalsDict[performance].keys()), dtype=int)[:-1]
            elif phaseType == 'phase':
                return list(self.phasesDict[performance].keys())
            elif phaseType == 'semiphase':
                return list(self.semiphasesDict[performance].keys())
            else:
                raise ValueError("Unexpected phase type", phaseType)
        else:
            phases = self.get_phasetype_keys('phase', performance, haveWaiting=True)[:-1]

            if phaseType == 'phase':
                return phases
            else:
                flatten2Dlist = lambda l: [item for sublist in l for item in sublist]
                return flatten2Dlist([self.get_phasetype_keys_from_phase(phase, phaseType, performance) for phase in phases])

    def get_phasetype_keys_from_phase(self, phase, phaseType, performance):
        if phaseType == 'phase':
            return [phase]
        if phaseType == 'interval':
            return list(range(*self.phasesDict[performance][phase]))
        elif phaseType == 'semiphase':
            semiphaseKey = phase[0]  # First letter of the phase word
            return [key for key in self.semiphasesDict[performance].keys() if semiphaseKey in key]
        else:
            raise ValueError('Unexpected phaseType', phaseType)

    def get_mouse_from_session(self, session, datatype):
        dataFrameKey = datatype if datatype == "deconv" else "dff"
        rows = self.get_rows(dataFrameKey, {"session": session})
        mice = list(rows["mousename"])
        assert len(mice) == 1, "Each session should be related to exactly one mouse"
        return mice[0]

    def get_rows(self, metaFrameName, query):
        return pandas_helper.pd_query(self.metaDataFrames[metaFrameName], query)

    def get_data_session_interval_fromindex(self, idxNeuro, idxBehavior, intervalStart, intervalEnd, dataType='raw'):
        if (idxNeuro is None) or (idxBehavior is None):
            raise ValueError("Unexpected indices", idxNeuro, idxBehavior)

        dataNeuroThis = self.dataNeuronal[dataType][idxNeuro]
        framesArrThis = self.behaviorStateFrames[idxBehavior]

        idxIntervalStart = intervalStart - 1
        idxIntervalEnd = intervalEnd - 1

        startFrameIdxs = framesArrThis[:, idxIntervalStart]
        endFrameIdxs = framesArrThis[:, idxIntervalEnd]

        return [np.copy(dataNeuroThis[start:end].T) for start, end in zip(startFrameIdxs, endFrameIdxs)]

    def get_data_from_interval(self, intervalStart, intervalEnd, queryDict):
        rez = []

        dataType = queryDict["datatype"] if "datatype" in queryDict.keys() else "raw"
        metaKey = "deconv" if dataType == "deconv" else "dff"

        queryDictNeuro = {k: v for k, v in queryDict.items() if k in ["mousename", "session"]}
        queryDictBehav = {k: v for k, v in queryDict.items() if k not in ["mousename", "datatype"]}

        rowsNeuro = self.get_rows(metaKey, queryDictNeuro)
        for idxNeuro, rowNeuro in rowsNeuro.iterrows():
            queryDictBehav["session"] = rowNeuro["session"]
            rowsBehavThis = self.get_rows('behaviorStates', queryDictBehav)

            if len(rowsBehavThis) == 0:
                if self.verbose:
                    print("No behaviour found for", queryDictBehav, "; skipping")
            else:
                for idxBehavior, rowBehavior in rowsBehavThis.iterrows():
                    rezThis = self.get_data_session_interval_fromindex(idxNeuro, idxBehavior, intervalStart, intervalEnd, dataType)

                    if np.any([np.prod(d.shape) == 0 for d in rezThis]):
                        if self.verbose:
                            print("Warning: ", rowNeuro["session"], "has zero interval", intervalStart, intervalEnd)

                    rez += rezThis

        return rez

    def get_data_from_phase(self, phase, queryDict):
        if "performance" in queryDict.keys():
            perf = queryDict["performance"]
            intervalStart, intervalEnd = self.phasesDict[perf][phase]
            return self.get_data_from_interval(intervalStart, intervalEnd, queryDict)
        else:
            # If performance is not specified, merge trials for all performance measures
            rez = []
            for perf in self.get_performances():
                rez += self.get_data_from_phase(phase, {**queryDict, **{"performance" : perf}})
            return rez

    def get_data_from_semiphase(self, semiphase, queryDict):
        if "performance" in queryDict.keys():
            perf = queryDict["performance"]
            intervalStart, intervalEnd = self.semiphasesDict[perf][semiphase]
            return self.get_data_from_interval(intervalStart, intervalEnd, queryDict)
        else:
            rez = []
            for perf in self.get_performances():
                if semiphase not in self.semiphasesDict[perf].keys():
                    raise ValueError('Semiphase', semiphase, 'does not exist for performance', perf)

                intervalStart, intervalEnd = self.semiphasesDict[perf][semiphase]
                rez += self.get_data_from_interval(intervalStart, intervalEnd, {**queryDict, **{"performance" : perf}})

            return rez

    # Wrapper for selecting phase or single interval
    def get_data_from_selector(self, selector, queryDict):
        selectorKey, selectorVal = *selector.keys(), *selector.values()

        if selectorKey == "interval":
            intervalStart = selectorVal
            return self.get_data_from_interval(intervalStart, intervalStart + 1, queryDict)
        elif selectorKey == "phase":
            return self.get_data_from_phase(selectorVal, queryDict)
        elif selectorKey == "semiphase":
            return self.get_data_from_semiphase(selectorVal, queryDict)
        elif selectorKey == "range":
            intervalStart, intervalEnd = selectorVal
            return self.get_data_from_interval(intervalStart, intervalEnd, queryDict)
        else:
            raise ValueError("Unexpected selector", selector)

    def get_phase_bounding_lines(self, phaseType, performance, haveWaiting=True):
        if phaseType == 'interval':
            intervBounds = []

            phases = self.get_phasetype_keys('phase', performance)
            if not haveWaiting:
                phases = phases[:-1]

            for phase in phases:
                intervBounds += list(self.phasesDict[performance][phase])
            return np.array(list(set(intervBounds))) - 0.5

        elif phaseType == 'semiphase':
            semiphases = self.get_phasetype_keys('semiphase', performance)
            lastE = [idx for idx, val in enumerate(semiphases) if 'E' in val][-1] + 1
            lastM = [idx for idx, val in enumerate(semiphases) if 'M' in val][-1] + 1
            lastR = [idx for idx, val in enumerate(semiphases) if 'R' in val][-1] + 1
            lastW = [idx for idx, val in enumerate(semiphases) if 'W' in val][-1] + 1

            if haveWaiting:
                return np.array([0, lastE, lastM, lastR, lastW]) + 0.5
            else:
                return np.array([0, lastE, lastM, lastR]) + 0.5
        else:
            raise ValueError("Unexpected phasetype", phaseType)
