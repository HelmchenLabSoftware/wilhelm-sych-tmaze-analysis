import numpy as np

from mesostat.utils.arrays import numpy_nonelist_to_array, get_list_shapes
import mesostat.stat.resampling as resampling

from mesostat.metric.metric_non_uniform import MetricCalculatorNonUniform

# Extract data list of shape [nTrial, nChannel, nTime] based on query
# Calculate metric as a list over repetiations, return.
# Note that if the loop is over repetitions, still a list of 1 result is returned
def _metric_single(selector, dataDB, queryDict, metricName, metricCalculator, dimOrderTrg, settings, channelFilter):
    queryMouseName = "mousename" in queryDict
    queryMouseKey = "session" in queryDict
    loopOverChannels = "p" in dimOrderTrg
    if (not queryMouseName) and (not queryMouseKey) and loopOverChannels:
        raise ValueError("Can't iterate over channels if more than 1 mouse selected")

    dataLst = dataDB.get_data_from_selector(selector, queryDict)
    if len(dataLst) == 0:
        return None
    else:
        # Extract settings
        zscoreChannel = False if "zscoreChannel" not in settings.keys() else settings["zscoreChannel"]
        metricSettings = None if "metricSettings" not in settings.keys() else settings["metricSettings"]
        sweepSettings = None if "sweepSettings" not in settings.keys() else settings["sweepSettings"]

        # Filter channels if requested
        if channelFilter is not None:
            if (not queryMouseName) and (not queryMouseKey):
                raise ValueError("Channel filter currently not implemented for sweeps over multiple mice")
            elif queryMouseName:
                mousename = queryDict["mousename"]
            else:
                mousename = dataDB.get_mouse_from_session(queryDict["session"], queryDict["datatype"])

            channelFilterThis = channelFilter[mousename]
            dataLst = [d[channelFilterThis] for d in dataLst]

        # Test if requested
        if "nTest" not in settings.keys():
            metricCalculator.set_data(dataLst, zscoreChannel=zscoreChannel)
            return metricCalculator.metric3D(metricName, dimOrderTrg, metricSettings=metricSettings, sweepSettings=sweepSettings)
        else:
            rezLst = []
            if settings["testMethod"] == "permChannel":
                for iTest in range(settings["nTest"]):
                    dataLst = [resampling.permutation(data, axis=0) for data in dataLst]
                    metricCalculator.set_data(dataLst, zscoreChannel=zscoreChannel)
                    rezLst += [metricCalculator.metric3D(metricName, dimOrderTrg, metricSettings=metricSettings, sweepSettings=sweepSettings)]

            if len(get_list_shapes(rezLst)) > 1:
                raise ValueError("Testing currently not available for metrics of variable shape")

            # Set testing as last dimension, so it does not interfere with postprocessing
            rezArr = np.array(rezLst)
            return rezArr.transpose(np.hstack([np.arange(1, rezArr.ndim), [0]]))


def metric_by_interval(dataDB, queryDict, metricName, dimOrderTrg, settings, intervals=None, channelFilter=None):
    serial = True if "serial" not in settings.keys() else settings["serial"]
    mc = MetricCalculatorNonUniform(serial=serial, verbose=False)

    if intervals is None:
        intervals = range(len(dataDB.metaDataFrames['interval_maps'][queryDict["performance"]]) - 1)

    metricValues = []   # [nInterv, nDataPoint]
    for iInterv in intervals:
        selector = {"interval" : iInterv}
        metricValues += [_metric_single(selector, dataDB, queryDict, metricName, mc, dimOrderTrg, settings, channelFilter)]

    return np.array(metricValues)


def metric_by_phase(dataDB, queryDict, metricName, dimOrderTrg, settings, phases=None, channelFilter=None):
    serial = True if "serial" not in settings.keys() else settings["serial"]
    mc = MetricCalculatorNonUniform(serial=serial, verbose=False)

    if phases is None:
        phases = dataDB.phases["Correct"]

    metricValues = []
    for phase in phases:
        selector = {"phase" : phase}
        metricValues += [_metric_single(selector, dataDB, queryDict, metricName, mc, dimOrderTrg, settings, channelFilter)]

    return np.array(metricValues)


def metric_by_selector(dataDB, queryDict, metricName, dimOrderTrg, selector, settings, channelFilter=None):
    serial = True if "serial" not in settings.keys() else settings["serial"]
    mc = MetricCalculatorNonUniform(serial=serial, verbose=False)
    return _metric_single(selector, dataDB, queryDict, metricName, mc, dimOrderTrg, settings, channelFilter)

# A sweep over query parameters. If mice are not specified, then we iterate over all mice
# However, we do not stack data from different mice. Instead we stack results from each mouse
# Mice can have different number or meaning of channels, so stacking may behave unexpectedly or break down
def metric_by_sweep(dataDB, sweepDF, metricName, dimOrderTrg, selector, settings, channelFilter=None):
    serial = True if "serial" not in settings.keys() else settings["serial"]
    mc = MetricCalculatorNonUniform(serial=serial, verbose=False)

    dict_drop_value = lambda d, val: {k: v for k, v in d.items() if v != val}

    metricValues = []
    for idx, row in sweepDF.iterrows():
        queryDict = dict_drop_value(row.to_dict(), "All")
        metricValues += [_metric_single(selector, dataDB, queryDict, metricName, mc, dimOrderTrg, settings, channelFilter)]

    return metricValues
