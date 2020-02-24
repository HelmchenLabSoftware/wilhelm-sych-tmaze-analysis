import numpy as np

def set_diag_zero(M):
    M[np.eye(M.shape[0], dtype=bool)] = 0
    return M