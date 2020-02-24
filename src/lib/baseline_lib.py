import numpy as np

def zscore(x):
    return (x - np.mean(x)) / np.std(x)

# Convert data to probability distribution
# By cropping lowest percentile
def crop_quantile(x, p):
    xnorm = x - np.quantile(x, p)  # Subtract percentile to arrive at baseline
    xnorm[xnorm < 0] = 0
    return xnorm / np.mean(xnorm)  # Normalize
