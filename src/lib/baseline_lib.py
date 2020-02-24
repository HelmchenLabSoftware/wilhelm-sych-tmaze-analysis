import numpy as np

def zscore(x):
    return (x - np.mean(x)) / np.std(x)

# Convert data to probability distribution
# By cropping lowest percentile
def crop_percentile(x, p):
    # xnorm = x - np.min(x)         # Subtract minimum
    xnorm = x - np.percentile(x, p)  # Subtract percentile to arrive at baseline
    xnorm[xnorm < 0] = 0
    return xnorm / np.sum(xnorm)  # Normalize
