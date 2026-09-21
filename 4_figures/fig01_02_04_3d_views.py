#!/usr/bin/env python3
"""Figures 1, 2 and 4: three-dimensional views of the fracture networks. Unlike the other figure scripts these
need the raw simulation output under $WORK (meshes and VTK snapshots are not distributed) and pyvista.

Figure 1  log10 of the Darcy velocity magnitude in DFN6, DFN9, DFN15 and DFN19
Figure 2  initial calcite and pyrite enrichment relative to the network mean in DFN9 (anti-aligned / aligned)
Figure 4  pyrite volume fraction at 1,500 yr relative to its initial uniform value in DFN9 (best / homogeneous / worst)
"""
import os, sys
import numpy as np
import pyvista as pv
import matplotlib.image as mpimg
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "3_extract"))
from paper import *
from common import WORK, read_vtk_field, cell_volumes, darcy_speed

pv.OFF_SCREEN = True; pv.set_plot_theme("document")
SEED = NETWORKS["seed"]; DOMAIN = (-3.0, 3.0, -1.75, 1.75, -1.0, 1.0); TMP = os.path.join(OUT, "_panels"); os.makedirs(TMP, exist_ok=True)


def mesh(seed):
    """Triangulated fracture surfaces from the AVS UCD file written by dfnWorks (values live on the nodes = PFLOTRAN cells)."""
    L = open(os.path.join(WORK, "networks", "s%d" % seed, "full_mesh.inp")).read().split("\n"); n, c = (int(x) for x in L[0].split()[:2])
    pts = np.array([[float(x) for x in L[1 + i].split()[1:4]] for i in range(n)])
    tri = np.array([[int(x) - 1 for x in L[1 + n + i].split()[3:6]] for i in range(c) if L[1 + n + i].split()[2] == "tri"])
    return pv.PolyData(pts, np.hstack([np.full((len(tri), 1), 3), tri]).ravel())


def render(m, values, cmap, clim, name):
    m = m.copy(); m.point_data["f"] = values
    pl = pv.Plotter(off_screen=True, window_size=[1600, 1360])
    pl.add_mesh(m, scalars="f", cmap=cmap, clim=clim, show_scalar_bar=False)
    pl.add_mesh(pv.Box(DOMAIN), style="wireframe", color="black", line_width=3)
    pl.view_vector((0.6, -1.0, 0.5), viewup=(0, 0, 1)); pl.reset_camera(); pl.camera.zoom(1.15)
    out = os.path.join(TMP, name + ".png"); pl.screenshot(out, scale=2, transparent_background=True); pl.close(); return out


def compose(panels, titles, cmap, clim, label, name, ncol):
    nrow = -(-len(panels) // ncol); fig, ax = plt.subplots(nrow, ncol, figsize=(4.2 * ncol + 0.8, 3.4 * nrow), squeeze=False)
    for a, png, t in zip(ax.ravel(), panels, titles): a.imshow(mpimg.imread(png)); a.axis("off"); a.set_title(t, fontsize=10)
    sm = plt.cm.ScalarMappable(norm=plt.Normalize(*clim), cmap=cmap); sm.set_array([]); fig.colorbar(sm, ax=ax, shrink=0.8, pad=0.01).set_label(label)
    save(fig, name)


def snapshot(seed, bp, bc, k, field, n):
    return np.clip(read_vtk_field(os.path.join(WORK, "rxn", "factorial", "s%d" % seed, "pyr%+d_cal%+d" % (bp, bc), "deck-%03d.vtk" % k), field, n), 0, None)


# ---- Figure 1
panels = []
for dfn in (6, 9, 15, 19):
    n, vol = cell_volumes(SEED[dfn]); v = darcy_speed(SEED[dfn], n)
    panels.append(render(mesh(SEED[dfn]), np.log10(np.clip(v, v[v > 0].min(), None)), "turbo", [-3, 0], "fig01_dfn%d" % dfn))
compose(panels, ["(%s) DFN%d, $\\sigma_{\\log|q|}$ = %.2f, %d fractures" % (l, d, SIGMA[d], NETWORKS["n_fractures"][d]) for l, d in zip("abcd", (6, 9, 15, 19))],
        "turbo", [-3, 0], "$\\log_{10}|v|$ (m yr$^{-1}$)", "fig01_networks", 4)

# ---- Figures 2 and 4 (DFN9)
seed = SEED[9]; n, vol = cell_volumes(seed); m9 = mesh(seed)
panels, titles = [], []
for field, mineral, (lo, hi) in (("Calcite_VF", "calcite", ((0, -2), (0, 2))), ("Pyrite_VF", "pyrite", ((-1, 0), (1, 0)))):
    for (bp, bc), word in zip((lo, hi), ("anti-aligned", "aligned")):
        vf = snapshot(seed, bp, bc, 0, field, n); e = vf / ((vol * vf).sum() / vol.sum())
        panels.append(render(m9, np.log10(np.clip(e, 1e-3, None)), "RdBu_r", [-1.2, 1.2], "fig02_%s_%s" % (mineral, word))); titles.append("%s, %s" % (mineral, word))
compose(panels, titles, "RdBu_r", [-1.2, 1.2], "$\\log_{10}$(VF / network-mean VF)", "fig02_arrangement", 2)

panels, titles = [], []
for (bp, bc), word in zip((BEST, HOMOGENEOUS, WORST), ("best", "homogeneous", "worst")):
    p0, p = snapshot(seed, 0, 0, 0, "Pyrite_VF", n), snapshot(seed, bp, bc, 15, "Pyrite_VF", n)
    rem = (vol * p).sum() / (vol * snapshot(seed, bp, bc, 0, "Pyrite_VF", n)).sum()
    panels.append(render(m9, p / p0, "RdBu_r", [0, 2], "fig04_%s" % word)); titles.append("(%s) %s: %.0f%% of pyrite remaining" % ("abc"[len(titles)], word, 100 * rem))
compose(panels, titles, "RdBu_r", [0, 2], "VF$_{pyr}$ / VF$_{pyr,0}$", "fig04_pyrite_remaining", 3)
