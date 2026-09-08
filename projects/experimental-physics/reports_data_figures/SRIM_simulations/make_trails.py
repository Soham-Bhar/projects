#!/usr/bin/env python
r"""
make_trails.py -- Ion-trail / collision-cascade visualisations from the low-N
trajectory runs produced by run_viz.py (results\viz\<run>\).

Per run it makes 4 figures -> results\plots\viz\<run>\ :
  1. exyz_trails_2d  : primary-ion paths (EXYZ.txt), depth-vs-lateral, coloured by energy
  2. exyz_trails_3d  : the same paths in true 3D
  3. cascade_2d      : ion tracks + full recoil spray (COLLISON.txt), depth-vs-lateral
  4. cascade_3d      : recoil cascade in 3D

Requires numpy + matplotlib only.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection

ROOT = os.path.dirname(os.path.abspath(__file__))
VIZ  = os.path.join(ROOT, "results", "viz")
OUT  = os.path.join(ROOT, "results", "plots", "viz")

RUNS = {
    "01_H_10MeV_Be":    dict(tag="H 10 MeV → Be",   E0=10000., atoms={4: "Be"}),
    "02_H_2MeV_Si":     dict(tag="H 2 MeV → Si",    E0=2000.,  atoms={14: "Si"}),
    "03_B_200keV_Si":   dict(tag="B 200 keV → Si",  E0=200.,   atoms={14: "Si", 5: "B"}),
    "07_Ar_10keV_GaAs": dict(tag="Ar 10 keV → GaAs",E0=10.,    atoms={31: "Ga", 33: "As", 18: "Ar"}),
}
ATOM_COLOR = {"Si": "C1", "Be": "C1", "Ga": "C0", "As": "C3", "B": "C2", "Ar": "C4"}

plt.rcParams.update({"figure.dpi": 130, "savefig.dpi": 160, "font.size": 10})


def pick_unit(max_A):
    """Return (divisor, label) so the depth axis reads in sensible units."""
    if max_A >= 1e4:
        return 1e4, "µm"
    if max_A >= 100:
        return 10., "nm"
    return 1., "Å"


# ---------------------------------------------------------------- EXYZ ---------
def parse_exyz(path):
    """ion# -> Nx4 array [X, Y, Z, E(keV)] (in Å)."""
    ions = {}
    with open(path, "r", encoding="latin-1") as f:
        for ln in f:
            t = ln.split()
            if len(t) < 5 or not t[0].isdigit() or len(t[0]) < 5:
                continue
            try:
                num = int(t[0]); E = float(t[1]); X = float(t[2]); Y = float(t[3]); Z = float(t[4])
            except ValueError:
                continue
            ions.setdefault(num, []).append((X, Y, Z, E))
    return {k: np.array(v) for k, v in ions.items() if len(v) > 1}


# ------------------------------------------------------------- COLLISON --------
def parse_collison(path):
    """Return (ion_paths, recoils).
       ion_paths: ion# -> Nx3 [X,Y,Z] (Å) from cascade-start (ion collision) points.
       recoils:   Nx4 [X,Y,Z,Zatom] for every recoil atom in every cascade."""
    ion_paths = {}
    recoils = []
    cur = None
    with open(path, "r", encoding="latin-1") as f:
        for raw in f:
            # strip the box-drawing delimiters (any non-ASCII byte) to spaces
            ln = "".join(c if 32 <= ord(c) < 127 else " " for c in raw)
            t = ln.split()
            if "Start of New Cascade" in ln:
                # ion collision point:  ion# E X Y Z Se Atom recoilE
                try:
                    cur = int(t[0]); X = float(t[2]); Y = float(t[3]); Z = float(t[4])
                    ion_paths.setdefault(cur, []).append((X, Y, Z))
                except (ValueError, IndexError):
                    pass
                continue
            if "Summary" in ln or "Recoil Atom Energy" in ln:
                continue
            # recoil row:  recoil# Zatom energy X Y Z vac repl
            if len(t) >= 6 and t[0].isdigit() and t[1].isdigit() and len(t[1]) <= 3:
                try:
                    Zat = int(t[1]); E = float(t[2]); X = float(t[3]); Y = float(t[4]); Z = float(t[5])
                except ValueError:
                    continue
                recoils.append((X, Y, Z, Zat))
    ion_paths = {k: np.array(v) for k, v in ion_paths.items() if len(v) > 1}
    return ion_paths, np.array(recoils) if recoils else np.empty((0, 4))


# --------------------------------------------------------------- plots ---------
def save(run, name):
    d = os.path.join(OUT, run); os.makedirs(d, exist_ok=True)
    p = os.path.join(d, name + ".png"); plt.savefig(p, bbox_inches="tight"); plt.close(); return p


def exyz_2d(run, meta, ions):
    allX = np.concatenate([a[:, 0] for a in ions.values()])
    div, unit = pick_unit(allX.max())
    segs, cols = [], []
    for a in ions.values():
        xy = np.column_stack([a[:, 0] / div, a[:, 1] / div])
        segs.extend([[xy[i], xy[i + 1]] for i in range(len(xy) - 1)])
        cols.extend(a[:-1, 3])
    lc = LineCollection(segs, cmap="plasma", array=np.array(cols), linewidths=0.7, alpha=0.8)
    lc.set_clim(0, meta["E0"])
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.add_collection(lc); ax.autoscale()
    ax.set_xlabel(f"Depth  X ({unit})"); ax.set_ylabel(f"Lateral  Y ({unit})")
    ax.set_title(f"{meta['tag']}: ion trails  ({len(ions)} ions, coloured by energy)")
    ax.grid(alpha=0.25)
    fig.colorbar(lc, ax=ax, label="Ion energy (keV)")
    return save(run, "1_exyz_trails_2d")


def exyz_3d(run, meta, ions):
    allX = np.concatenate([a[:, 0] for a in ions.values()])
    div, unit = pick_unit(allX.max())
    fig = plt.figure(figsize=(8, 6)); ax = fig.add_subplot(projection="3d")
    for a in ions.values():
        p = a[:, :3] / div
        segs = np.stack([p[:-1], p[1:]], axis=1)
        lc = Line3DCollection(segs, cmap="plasma", linewidths=0.6, alpha=0.8)
        lc.set_array(a[:-1, 3]); lc.set_clim(0, meta["E0"])
        ax.add_collection3d(lc)
    P = np.concatenate([a[:, :3] for a in ions.values()]) / div
    ax.set_xlim(P[:, 0].min(), P[:, 0].max())
    ax.set_ylim(P[:, 1].min(), P[:, 1].max())
    ax.set_zlim(P[:, 2].min(), P[:, 2].max())
    ax.set_xlabel(f"Depth X ({unit})"); ax.set_ylabel(f"Y ({unit})"); ax.set_zlabel(f"Z ({unit})")
    ax.set_title(f"{meta['tag']}: ion trails in 3D  ({len(ions)} ions)")
    ax.view_init(elev=18, azim=-70)
    return save(run, "2_exyz_trails_3d")


def cascade_2d(run, meta, ipaths, recoils, exyz):
    src = exyz if exyz else None
    allX = recoils[:, 0] if len(recoils) else np.concatenate([a[:, 0] for a in ipaths.values()])
    div, unit = pick_unit(allX.max() if len(allX) else 1)
    fig, ax = plt.subplots(figsize=(8, 5))
    # recoil spray, coloured by atom species
    if len(recoils):
        for Zat in np.unique(recoils[:, 3]).astype(int):
            m = recoils[:, 3] == Zat
            name = meta["atoms"].get(Zat, f"Z{Zat}")
            ax.scatter(recoils[m, 0] / div, recoils[m, 1] / div, s=1.5, alpha=0.25,
                       color=ATOM_COLOR.get(name, "0.5"), edgecolors="none",
                       label=f"{name} recoils ({m.sum()})")
    # ion tracks on top (smooth EXYZ path if available, else COLLISON ion points)
    tracks = exyz if exyz else ipaths
    for a in tracks.values():
        ax.plot(a[:, 0] / div, a[:, 1] / div, color="k", lw=0.35, alpha=0.6)
    ax.set_xlabel(f"Depth  X ({unit})"); ax.set_ylabel(f"Lateral  Y ({unit})")
    ax.set_title(f"{meta['tag']}: collision cascade (ion tracks + recoil spray)")
    ax.grid(alpha=0.25)
    lg = ax.legend(fontsize=8, markerscale=4, framealpha=0.9)
    for h in lg.legend_handles:
        try: h.set_alpha(1)
        except Exception: pass
    return save(run, "3_cascade_2d")


def cascade_3d(run, meta, ipaths, recoils, exyz):
    allX = recoils[:, 0] if len(recoils) else np.array([1.0])
    div, unit = pick_unit(allX.max())
    fig = plt.figure(figsize=(8, 6)); ax = fig.add_subplot(projection="3d")
    if len(recoils):
        # subsample for a legible 3D cloud
        n = len(recoils)
        idx = np.random.default_rng(0).choice(n, size=min(n, 12000), replace=False)
        for Zat in np.unique(recoils[:, 3]).astype(int):
            m = recoils[idx, 3] == Zat
            name = meta["atoms"].get(Zat, f"Z{Zat}")
            ax.scatter(recoils[idx][m, 0] / div, recoils[idx][m, 1] / div, recoils[idx][m, 2] / div,
                       s=2, alpha=0.15, color=ATOM_COLOR.get(name, "0.5"), edgecolors="none",
                       label=f"{name} recoils")
    tracks = exyz if exyz else ipaths
    for a in list(tracks.values()):
        ax.plot(a[:, 0] / div, a[:, 1] / div, a[:, 2] / div, color="k", lw=0.3, alpha=0.5)
    ax.set_xlabel(f"Depth X ({unit})"); ax.set_ylabel(f"Y ({unit})"); ax.set_zlabel(f"Z ({unit})")
    ax.set_title(f"{meta['tag']}: cascade in 3D")
    ax.view_init(elev=18, azim=-70)
    lg = ax.legend(fontsize=8, markerscale=3)
    for h in lg.legend_handles:
        try: h.set_alpha(1)
        except Exception: pass
    return save(run, "4_cascade_3d")


def main():
    total = 0
    for run, meta in RUNS.items():
        d = os.path.join(VIZ, run)
        if not os.path.isdir(d):
            print(f"{run:20s} -- missing, skip"); continue
        exyz = parse_exyz(os.path.join(d, "EXYZ.txt"))
        ipaths, recoils = parse_collison(os.path.join(d, "COLLISON.txt"))
        made = []
        if exyz:
            made.append(exyz_2d(run, meta, exyz))
            made.append(exyz_3d(run, meta, exyz))
        made.append(cascade_2d(run, meta, ipaths, recoils, exyz))
        made.append(cascade_3d(run, meta, ipaths, recoils, exyz))
        total += len(made)
        print(f"{run:20s} -> {len(made)} figs  (EXYZ ions={len(exyz)}, recoils={len(recoils)})")
    print(f"\nTOTAL {total} figures under {OUT}")


if __name__ == "__main__":
    main()
