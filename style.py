"""Oded's house plotting style (reconstructed to the same spec as ~/Desktop/regev/style.py):
plain matplotlib defaults, square panel, Blues colormap, LogNorm counts,
n= / r= as bare text in the upper-left corner, filenames as a two-line centred title.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def hist2d_compare(x, y, file1, file2, xlabel=None, ylabel=None, bins=60, note=None, ax=None):
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(4.2, 4.2))
    m = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x)[m], np.asarray(y)[m]
    h = ax.hist2d(x, y, bins=bins, cmap="Blues", norm=LogNorm())
    r = np.corrcoef(x, y)[0, 1]
    txt = f"n={len(x)}\nr={r:.3f}" + (f"\n{note}" if note else "")
    ax.text(0.02, 0.98, txt, transform=ax.transAxes, va="top", ha="left")
    ax.set_xlabel(xlabel or file1)
    ax.set_ylabel(ylabel or file2)
    ax.set_title(f"{file1} vs\n{file2}")
    ax.set_box_aspect(1)
    if own:
        fig.colorbar(h[3], ax=ax)
        return fig, ax
    return h, r


def square(ax):
    ax.set_box_aspect(1)
    return ax
