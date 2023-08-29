import numpy as np

# Cycle data backwards by n steps
def cycle_data(x, n):
    return np.hstack((x[n:], x[:n]))



