import numpy as np
import pandas as pd

from mesostat.utils.pandas_helper import get_rows_colvals

class SignificantCells:
    def __init__(self, filename):
        self.df = pd.read_hdf(filename, key='df')

    def get_cells_by_mouse(self):
        return { mousename :
            self._get_cells({"mousename" : mousename, "direction" : "All", "performance" : "All"})
                 for mousename in set(self.df["mousename"])
         }

    def _get_cells(self, queryDict):
        queryDF = get_rows_colvals(self.df, queryDict)
        if len(queryDF) == 0:
            raise ValueError("Nothing found for query", queryDict)

        cellsUnion = set()
        for cells in queryDF['cells']:
            cellsUnion |= set(cells)

        return np.array(list(cellsUnion))
