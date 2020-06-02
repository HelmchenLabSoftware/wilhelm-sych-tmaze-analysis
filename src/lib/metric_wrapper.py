import numpy as np

from mesostat.utils.arrays import numpy_nonelist_to_array, get_list_shapes
import mesostat.stat.resampling as resampling

from mesostat.metric.metric_non_uniform import MetricCalculatorNonUniform

# Extract data list of shape [nTrial, nChannel, nTime] based on query
# Calculate metric as a list over repetiations, return.
# Note that if the loop is over repetitions, still a list of 1 result is returned
def _metric_single(selector, dataDB, queryDict, metricName, metricCalculator, dimOrderTrg, settings):
    queryMouseName = "mousename" in queryDict
    queryMouseKey = "session" in queryDict
    loopOverChannels = "p" in dimOrderTrg
    if (not queryMouseName) and (not queryMouseKey) and loopOverChannels:
        raise ValueError("Can't iterate over channels if more than 1 mouse selected")

    if "interval" in selector:
        iInterv = selector["interval"]
        dataLst = dataDB.get_data_from_interval(iInterv, iInterv + 1, queryDict)
    elif "phase" in selector:
        dataLst = dataDB.get_data_from_phase(selector["phase"], queryDict)
    else:
        raise ValueError("Unexpected selector", selector)

    if len(dataLst) == 0:
        return None
    else:
        # Extract settings
        zscoreChannel = None if "zscoreChannel" not in settings.keys() else settings["zscoreChannel"]
        metricSettings = None if "metricSettings" not in settings.keys() else settings["metricSettings"]
        sweepSettings = None if "sweepSettings" not in settings.keys() else settings["sweepSettings"]

        # Preprocess data if
        if "nTest" not in settings.keys():
            metricCalculator.set_data(dataLst, zscoreChannel)
            return metricCalculator.metric3D(metricName, dimOrderTrg, metricSettings=metricSettings, sweepSettings=sweepSettings)
        else:
            rezLst = []
            if settings["testMethod"] == "permChannel":
                for iTest in range(settings["nTest"]):
                    dataLst = [resampling.permutation(data, axis=0) for data in dataLst]
                    metricCalculator.set_data(dataLst, zscoreChannel)
                    rezLst += [metricCalculator.metric3D(metricName, dimOrderTrg, metricSettings=metricSettings, sweepSettings=sweepSettings)]

            if len(get_list_shapes(rezLst)) > 1:
                raise ValueError("Testing currently not available for metrics of variable shape")

            # Set testing as last dimension, so it does not interfere with postprocessing
            rezArr = np.array(rezLst)
            return rezArr.transpose(np.hstack([np.arange(1, rezArr.ndim), [0]]))


# If multiplexKey is None, then this function is equivalent to _metric_single
# Otherwise, compute metric over sessions ("session") or animals ("mousename"), and stitch together as if they were trials
def _metric_single_multiplexed(selector, dataDB, queryDict, metricName, metricCalculator, dimOrderTrg, settings, multiplexKey):
    def _metric_single_wrapper(queryThis):
        return _metric_single(selector, dataDB, queryThis, metricName, metricCalculator, dimOrderTrg, settings)

    if multiplexKey is None:
        return _metric_single_wrapper(queryDict)
    elif multiplexKey == "mousename":
        multiplexValLst = dataDB.mice
    elif multiplexKey == "session":
        multiplexValLst = dataDB.sessions
    else:
        raise ValueError("Unexpected Multiplex Key", multiplexKey)

    rezLst = []
    for multiplexVal in multiplexValLst:
        queryThis = queryDict.copy()
        queryThis[multiplexKey] = multiplexVal
        rezLst += [_metric_single_wrapper(queryThis)]

    # if rezLst[0].ndim == 0:
    #     return np.array(rezLst)
    # else:
    #     return np.concatenate(rezLst, axis=0)

    # if "r" in dimOrderTrg and multiplexKey is not None:
    #     raise ValueError("It does not make sense to multiplex over", multiplexKey, "if already looping over repetitions")

    if "r" not in dimOrderTrg:
        # Multiplex will now be considered a repetitions dimension
        # Because there is no real repetitions dimension
        return numpy_nonelist_to_array(rezLst)
    else:
        # We must merge multiplex iterations with repetitions
        # We have to be able to iterate separately over animals and repetitions because of data preprocessing in mesostat
        # For example, we cannot zscore data by channel if data has more than one mouse
        # NOTE: None indicates lack of trials, so dropping None's does not affect the number of trials
        return np.concatenate([rez for rez in rezLst if rez is not None], axis=0)


def metric_by_interval(dataDB, queryDict, metricName, dimOrderTrg, settings, intervals=None, multiplexKey=None):
    serial = None if "serial" not in settings.keys() else settings["serial"]
    mc = MetricCalculatorNonUniform(serial=serial, verbose=False)

    if intervals is None:
        intervals = range(len(dataDB.metaDataFrames['interval_maps'][queryDict["performance"]]) - 1)

    metricValues = []   # [nInterv, nDataPoint]
    for iInterv in intervals:
        selector = {"interval" : iInterv}
        metricValues += [_metric_single_multiplexed(selector, dataDB, queryDict, metricName, mc, dimOrderTrg, settings, multiplexKey)]

    return np.array(metricValues)


def metric_by_phase(dataDB, queryDict, metricName, dimOrderTrg, settings, phases=None, multiplexKey=None):
    serial = None if "serial" not in settings.keys() else settings["serial"]
    mc = MetricCalculatorNonUniform(serial=serial, verbose=False)

    if phases is None:
        phases = dataDB.phases["Correct"]

    metricValues = []
    for phase in phases:
        selector = {"phase" : phase}
        metricValues += [_metric_single_multiplexed(selector, dataDB, queryDict, metricName, mc, dimOrderTrg, settings, multiplexKey)]

    return np.array(metricValues)


# A sweep over query parameters. If mice are not specified, then we iterate over all mice
# However, we do not stack data from different mice. Instead we stack results from each mouse
# Mice can have different number or meaning of channels, so stacking may behave unexpectedly or break down
def metric_by_sweep(dataDB, sweepDF, metricName, dimOrderTrg, selector, settings, multiplexKey=None):
    serial = None if "serial" not in settings.keys() else settings["serial"]
    mc = MetricCalculatorNonUniform(serial=serial, verbose=False)

    dict_drop_value = lambda d, val: {k: v for k, v in d.items() if v != val}

    metricValues = []
    for idx, row in sweepDF.iterrows():
        queryDict = dict_drop_value(row.to_dict(), "All")
        metricValues += [_metric_single_multiplexed(selector, dataDB, queryDict, metricName, mc, dimOrderTrg, settings, multiplexKey)]

    return metricValues
