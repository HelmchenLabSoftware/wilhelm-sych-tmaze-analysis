import numpy as np


# Calculate number of channels with non-zero STD over time
def num_non_zero_std(data2DLst, settings):
    for data2D in data2DLst:
        if data2D.shape[1] <= 1:
            raise ValueError("Must have sufficient timesteps to estimate variance, have", data2D.shape[1])

    # Mean over timesteps
    meanNonZeroStd = [np.nanmean(np.nanstd(data, axis=1) > 1.0e-10) for data in data2DLst]

    # Mean over trials
    return np.nanmean(meanNonZeroStd)


# Calculate mean number of timesteps in each trial
def num_sample(data2DLst, settings):
    return np.mean([data.shape[1] for data in data2DLst])
