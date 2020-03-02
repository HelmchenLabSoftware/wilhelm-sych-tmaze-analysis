import numpy as np

def get_off_diag(M):
    return M[~np.eye(M.shape[0], dtype=bool)]

def set_diag_zero(M):
    M[np.eye(M.shape[0], dtype=bool)] = 0
    return M