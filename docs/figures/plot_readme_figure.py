#!/usr/bin/env python3
# Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
"""Figure of the README: the effective potential of the Morris-Thorne wormhole,

    V(x) = L(L + 1)/(1 + x^2) + (1 - s^2)/(1 + x^2)^2,

as written in QNM-WKB3.mw, x being the proper radial distance from the throat in units of the
throat radius, on both sides of the throat. From this folder,
    python3 plot_readme_figure.py
writes morristhorne-potential.png (Python 3 with NumPy and Matplotlib). The style is the stamped
copy manuscript_style.py beside it.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import manuscript_style as hs  # noqa: E402


def V(x, L, s):
    return L * (L + 1) / (1 + x**2) + (1 - s**2) / (1 + x**2) ** 2


hs.configure(family="libertine")
plt.rcParams["savefig.dpi"] = 260
plt.rcParams["figure.facecolor"] = hs.PAPER
plt.rcParams["savefig.facecolor"] = hs.PAPER

fig, (left, right) = plt.subplots(1, 2, figsize=(hs.WIDTH, 2.6))
x = np.linspace(-6.0, 6.0, 2401)
for axis, cases, title in ((left, [(2, 0), (2, 1), (2, 2)], r"(a) $L = 2$, the three spins"),
                           (right, [(0, 0), (1, 1), (2, 2)], r"(b) the lowest multipole $L = s$")):
    for k, (L, s) in enumerate(cases):
        style = hs.series(k)
        style["marker"] = None
        axis.plot(x, V(x, L, s), lw=1.1, label=rf"$s = {s}$, $L = {L}$", **style)
        print(f"L = {L}, s = {s}: V(0) = {V(0.0, L, s):g}")
    axis.axvline(0.0, color=hs.SLATE, lw=0.6, linestyle=(0, (1.6, 1.6)), zorder=0)
    axis.set_xlim(-6.0, 6.0)
    axis.set_ylim(0.0, 7.6)
    axis.set_xlabel(r"$x$ (proper distance from the throat, in throat radii)")
    axis.set_ylabel(r"$V(x)$")
    hs.panel(axis, title)
    hs.legend_outside(axis, where="below", ncol=3)

out = HERE / "morristhorne-potential.png"
hs.save(fig, out)
# The background is opaque: without its alpha channel, and compressed a second time without
# loss, the file is smaller.
Image.open(out).convert("RGB").save(out, optimize=True)
print("written", out)
