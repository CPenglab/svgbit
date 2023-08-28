from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from libpysal.weights import W as libpysal_W
from scipy.stats import norm


def pysal_to_pandas(W: libpysal_W) -> pd.DataFrame:
    return pd.DataFrame(*W.full()).astype(int).T


def plot_gmm(gmm, hist_data, save_path) -> None:
    """Plot a GMM given by sklearn."""
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.hist(hist_data, bins=100, density=True, color="#f0a1a8")
    x = np.linspace(min(hist_data), max(hist_data), 5000)
    ax.plot(
        x,
        np.exp(gmm.score_samples(x.reshape(-1, 1))),
        lw=4,
        label="GMM",
    )
    for i in range(gmm.get_params()["n_components"]):
        ax.plot(
            x,
            norm.pdf(
                x,
                gmm.means_[i, 0],
                gmm.covariances_[i, 0]**(1 / 2),
            ) * gmm.weights_[i],
            lw=4,
            ls="--",
            label=f"Gaussian {i}, weight {gmm.weights_[i]:.3f}",
        )
    ax.legend()
    fig.savefig(save_path)


def get_cmap(length):
    return plt.colormaps["hsv"].resampled(length)
