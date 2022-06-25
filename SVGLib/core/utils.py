import pandas as pd
from libpysal.weights import W as libpysal_W


def pysal_to_pandas(W: libpysal_W) -> pd.DataFrame:
    return pd.DataFrame(*W.full()).astype(int).T
