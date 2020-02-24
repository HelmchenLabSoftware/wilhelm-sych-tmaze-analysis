import numpy as np
import pandas as pd

# G
def get_rows_colval(df, colname, val):
    return df[df[colname] == val]


# Get rows for which several columns have some exact values
def get_rows_colvals(df, queryDict):
    # Query likes strings to be wrapped in quotation marks for later evaluation
    strwrap = lambda val: '"' + val + '"' if isinstance(val, str) else str(val)
    query = ' and '.join([colname+'=='+strwrap(val) for colname, val in queryDict.items()])
    return df.query(query)

# Return all rows for which values match the list exactly
def get_rows_colvals_exact(df, lst):
    return get_rows_colvals(df, dict(zip(df.columns, lst)))

# Check if there is at least 1 row that matches the list exactly
def row_exists(df, lst):
    return len(get_rows_colvals_exact(df, lst)) > 0


# Add new row to dataframe, unless such a row is already present
def add_list_as_row(df, lst, skip_repeat=True):
    if skip_repeat and row_exists(df, lst):
        print("Skipping existing row", lst)
        return df
    else:
        newRow = pd.DataFrame([lst], columns=df.columns)
        return df.append(newRow, ignore_index=True)


def get_one_row(rows):
    nRows = rows.shape[0]
    if nRows == 0:
        return None, None
    elif nRows > 1:
        print(rows)
        raise ValueError("Expected 1 match, got", nRows)

    for idx, row in rows.iterrows():
        return idx, row