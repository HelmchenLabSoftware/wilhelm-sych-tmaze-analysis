import numpy as np
import pandas as pd
from IPython.display import display
from scipy.stats import combine_pvalues

from mesostat.utils.pandas_helper import outer_product_df, get_rows_colvals
from mesostat.stat.machinelearning import binary_classifier

from src.lib.metric_wrapper import metric_by_sweep, metric_by_phase, metric_by_selector
from src.lib.stat_lib import test_rank_wrapper


def _table_multiple_tests(sweepDF, metricValuesDict, multiplexKey):
    condValues = list(metricValuesDict.keys())
    nConditions = len(condValues)
    rezDF = sweepDF.copy()
    for idx, row in sweepDF.iterrows():
        # Print means for each condition
        for condVal in condValues:
            rezDF.at[idx, "mu("+condVal+")"] = np.nanmean(metricValuesDict[condVal][idx])

        # Perform tests for each pair of conditions
        nNoNanLst = []
        for i in range(nConditions):
            for j in range(i + 1, nConditions):
                testLabel = "-logp(" + condValues[i] + "->" + condValues[j] + ")"
                logp, nNoNan = test_rank_wrapper(
                    metricValuesDict[condValues[i]][idx],
                    metricValuesDict[condValues[j]][idx])
                rezDF.at[idx, testLabel] = np.round(logp, 2)
                nNoNanLst += [nNoNan]

        # rezDF.at[idx, "nTrial"] = str(nNoNanLst)

    rezDF = _consolidate_multiplexed_pvalues(rezDF, sweepDF, multiplexKey)
    display(rezDF)


def _consolidate_multiplexed_pvalues(testDF, sweepDF, multiplexKey):
    if multiplexKey is None:
        return testDF
    else:
        '''
        Algorithm:
        1. Drop multiplex from sweep
        2. Loop over rows of sweep
        3. Query each row of sweep to resultDF
        4. For result, assert all columns equal except pval. Merge pVal using Fisher's Method
        5. Stitch resulting rows
        '''

        rezDF = pd.DataFrame()
        sweepDF = sweepDF.drop([multiplexKey], axis=1).drop_duplicates()
        testDF = testDF.drop([multiplexKey], axis=1)

        for idxSweep, rowSweep in sweepDF.iterrows():
            queryRows = get_rows_colvals(testDF, dict(rowSweep))

            filterDict = {}
            for queryCol in queryRows:
                vals = list(set(queryRows[queryCol]))
                if "logp(" in queryCol:
                    # If this is a tested column, combine all pvalues using Fisher's Method
                    pVals = 10**(-np.array(vals))
                    pValTot = combine_pvalues(pVals)[1]
                    filterDict[queryCol] = -np.log10(pValTot)
                elif "mu(" in queryCol:
                    filterDict[queryCol] = -np.mean(vals)
                else:
                    # If this is not a tested column, all values there should be equal, so test and take the one
                    if len(vals) != 1:
                        raise ValueError("For key", queryCol, "expected 1 value, got", vals)
                    filterDict[queryCol] = vals[0]

            rezDF = rezDF.append(pd.DataFrame(filterDict, index=[0]))

        return rezDF.reset_index(drop=True)


def _outer_product_multiplexed(dataDB, sweepDict, multiplexKey):
    if multiplexKey is None:
        return outer_product_df(sweepDict)
    else:
        assert "mousename" not in sweepDict, "Pointless to multiplex if specific mouse selected"
        assert "session" not in sweepDict, "Pointless to multiplex if specific session selected"
        if multiplexKey == "mousename":
            return outer_product_df({**sweepDict, **{"mousename" : list(dataDB.mice)}})
        elif multiplexKey == "session":
            return outer_product_df({**sweepDict, **{"mousename": list(dataDB.sessions)}})
        else:
            raise ValueError("Unexpected multiplexKey", multiplexKey)


def table_discriminate_behavior(dataDB, selector, condition, sweepDict, metricName, trgDimOrder="r", settings=None, multiplexKey=None):
    if settings == None:
        settings = dict()
    sweepDF = _outer_product_multiplexed(dataDB, sweepDict, multiplexKey)
    condValues = list(set(dataDB.metaDataFrames['behaviorStates'][condition]))

    metricValuesDict = {}
    for i, condVal in enumerate(condValues):
        sweepDFThis = sweepDF.copy()
        sweepDFThis[condition] = condVal

        rez = metric_by_sweep(dataDB, sweepDFThis, metricName, trgDimOrder, selector, settings)

        # If metric is not scalar, just average over all values
        ndim = rez[0].ndim
        if ndim > 1:
            print("Warning: Using nonscalar metric of", ndim-1, "dimensions")
            axis = tuple(range(1, ndim))
            rez = [np.mean(r, axis=axis) for r in rez]

        metricValuesDict[condVal] = rez

    _table_multiple_tests(sweepDF, metricValuesDict, multiplexKey)


def table_discriminate_time(dataDB, sweepDict, selectors, metricName, trgDimOrder="r", settings=None, multiplexKey=None):
    if settings == None:
        settings = dict()
    sweepDF = _outer_product_multiplexed(dataDB, sweepDict, multiplexKey)

    metricValuesDict = {}
    for selectorKey, selectorVals in selectors.items():
        for selectorVal in selectorVals:
            selector = {selectorKey : selectorVal}
            rez = metric_by_sweep(dataDB, sweepDF, metricName, trgDimOrder, selector, settings)

            # If metric is not scalar, just average over all values
            ndim = rez[0].ndim
            if ndim > 1:
                print("Warning: Using nonscalar metric of", ndim-1, "dimensions")
                axis = tuple(range(1, ndim))
                rez = [np.mean(r, axis=axis) for r in rez]

            metricValuesDict[selectorVal] = rez

    _table_multiple_tests(sweepDF, metricValuesDict, multiplexKey)


def table_binary_classification(dataDB, phase, binaryDimension, metric, dimOrdTarget, queryDict, settings):
    paramValues = set(dataDB.metaDataFrames['behaviorStates'][binaryDimension])

    nMetricLabels = []
    metricLabels = []
    metricValsFlat = []
    for paramVal in paramValues:
        queryDict[binaryDimension] = paramVal
        metricVals = metric_by_selector(dataDB, queryDict, metric, dimOrdTarget, {"phase", phase}, settings)

        if metricVals[0].ndim > 1:
            metricValsFlat += [val.flatten() for val in metricVals]
        else:
            metricValsFlat += list(metricVals)
        metricLabels += [paramVal] * len(metricVals)
        nMetricLabels += [len(metricVals)]

    bc = binary_classifier(np.array(metricValsFlat), np.array(metricLabels), method="looc", balancing=True)
    bc['nLabels'] = tuple(nMetricLabels)

    return bc


def table_binary_classification_bymouse(dataDB, phase, binaryDimension, metric, dimOrdTarget, queryDict, settings):
    rezDF = pd.DataFrame(columns=["mousename", "nLabels", "acc_train", "acc_test", "acc_naive", "p-value"])

    for mousename in dataDB.mice:
        print("Doing mouse", mousename)
        queryDict["mousename"] = mousename

        bc = table_binary_classification(dataDB, phase, binaryDimension, metric, dimOrdTarget, queryDict, settings)
        bc['mousename'] = mousename
        rezDF = rezDF.append(bc, ignore_index=True)

    display(rezDF)
