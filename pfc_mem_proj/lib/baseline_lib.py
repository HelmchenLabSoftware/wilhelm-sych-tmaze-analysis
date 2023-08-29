import numpy as np

def zscore(x):
    return (x - np.mean(x)) / np.std(x)

# Crop lowest percentile, then z-score
def crop_quantile(x, p):
    xnorm = x - np.quantile(x, p)  # Subtract percentile to arrive at baseline
    xnorm[xnorm < 0] = 0
    return zscore(xnorm)