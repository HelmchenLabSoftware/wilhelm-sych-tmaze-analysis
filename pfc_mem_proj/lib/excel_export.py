import numpy as np
import pandas as pd
from typing import Union


def write_excel_1D(data: np.ndarray, trg: Union[str, pd.ExcelWriter], sheetName: str = None,
                   rownames: np.ndarray = None, colname: str = None):
    ds = pd.Series(data, index=rownames, name=colname)
    haveColNames = colname is not None
    haveRowNames = rownames is not None
    ds.to_excel(trg, sheet_name=sheetName, header=haveColNames, index=haveRowNames)


def write_excel_2D(data: np.ndarray, trg: Union[str, pd.ExcelWriter], sheetName: str = None,
                   rownames: np.ndarray = None, colnames: np.ndarray = None):
    ds = pd.DataFrame(data, index=rownames, columns=colnames)
    haveColNames = colnames is not None
    haveRowNames = rownames is not None
    ds.to_excel(trg, sheet_name=sheetName, header=haveColNames, index=haveRowNames)
