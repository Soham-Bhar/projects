#!/usr/bin/env python
r"""make_p5bc_plots.py -- Projects 5(b) and 5(c): compare the sputter yield of the
optimum argon beam with protons and gold at the same energy (100 keV) on GaAs.

Reads results\{5b_H_100keV_GaAs, 5c_Au_100keV_GaAs, 09_Ar_100keV_GaAs}\
and writes figures into results\plots\P5bc\.
"""
import os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = r"D:\SRIM_sim"
OUT = os.path.join(ROOT, "results", "plots", "P5bc")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"figure.dpi": 130, "savefig.dpi": 150, "font.size": 10,
                     "axes.grid": True, "grid.alpha": 0.3, "figure.autolayout": True})

# label, job, projectile mass (amu), colour
CASES = [("H",  "5b_H_100keV_GaAs",   1.008,  "C0"),
         ("Ar", "09_Ar_100keV_GaAs",  39.962, "C2"),
         ("Au", "5c_Au_100keV_GaAs",  196.967, "C1")]


def tdata(job, pat):
    t = open(os.path.join(ROOT, "results", job, "TDATA.txt"), encoding="latin-1").read()
    m = re.search(pat, t)
    return m.group(1).strip() if m else None


def sputter(job):
    """returns (Z array, energy eV array, cosX array) for ejected atoms only"""
    Z, E, cx = [], [], []
    p = os.path.join(ROOT, "results", job, "SPUTTER.txt")
    for ln in open(p, encoding="latin-1"):
        t = ln.split()
        if len(t) == 10 and t[0] == "S":
            try:
                Z.append(int(float(t[2]))); E.append(float(t[3])); cx.append(float(t[7]))
            except ValueError:
                continue
    return np.array(Z), np.array(E), np.array(cx)


rows = []
for lab, job, mass, col in CASES:
    N = int(tdata(job, r'Total Ions calculated =\s*([0-9]+)'))
    rp = float(tdata(job, r'Average Range\s+=\s+([0-9.E+-]+) Angstroms'))
    vac = float(tdata(job, r'Average Vacancy/Ion   =\s+([0-9.E+-]+)'))
    bs = int(tdata(job, r'Total Backscattered Ions=\s*([0-9]+)'))
    tr = int(tdata(job, r'Total Transmitted Ions  =\s*([0-9]+)'))
    Z, E, cx = sputter(job)
    nGa, nAs = int((Z == 31).sum()), int((Z == 33).sum())
    rows.append(dict(lab=lab, job=job, mass=mass, col=col, N=N, rp=rp, vac=vac,
                     bs=bs, tr=tr, Y=(nGa + nAs) / N, YGa=nGa / N, YAs=nAs / N,
                     E=E, Z=Z, cx=cx))

print(f"{'ion':>4} {'mass':>8} {'N':>7} {'Y':>7} {'Y_Ga':>6} {'Y_As':>6} {'As/Ga':>6} "
      f"{'Rp(A)':>8} {'vac/ion':>8} {'backsc':>7} {'transm':>7}")
for r in rows:
    print(f"{r['lab']:>4} {r['mass']:>8.3f} {r['N']:>7} {r['Y']:>7.3f} {r['YGa']:>6.3f} "
          f"{r['YAs']:>6.3f} {r['YAs']/max(r['YGa'],1e-9):>6.2f} {r['rp']:>8.1f} "
          f"{r['vac']:>8.1f} {r['bs']:>7} {r['tr']:>7}")

# ---- Fig 1: yield comparison bar chart --------------------------------------
fig, ax = plt.subplots(figsize=(7.0, 4.3))
x = np.arange(len(rows)); w = 0.35
ax.bar(x - w/2, [r["YGa"] for r in rows], w, label="Ga", color="C0")
ax.bar(x + w/2, [r["YAs"] for r in rows], w, label="As", color="C3")
for i, r in enumerate(rows):
    ax.text(i, max(r["YGa"], r["YAs"]) + 0.25, f"total\n{r['Y']:.3f}",
            ha="center", fontsize=9, weight="bold")
ax.set_xticks(x)
ax.set_xticklabels([f"{r['lab']}\n({r['mass']:.0f} amu)" for r in rows])
ax.set_ylabel("Sputter yield (atoms/ion)")
ax.set_title("Sputter yield from GaAs at 100 keV: H vs Ar vs Au")
ax.set_ylim(0, max(max(r["YAs"] for r in rows), 1) * 1.35)
ax.legend(fontsize=9)
plt.savefig(os.path.join(OUT, "1_yield_comparison.png")); plt.close()

# ---- Fig 2: yield vs projectile mass (log-log) ------------------------------
fig, ax = plt.subplots(figsize=(7.0, 4.3))
m = np.array([r["mass"] for r in rows]); Y = np.array([r["Y"] for r in rows])
ax.loglog(m, Y, "o-", color="C4", ms=8)
for r in rows:
    ax.annotate(f" {r['lab']}: {r['Y']:.3f}", (r["mass"], r["Y"]), fontsize=9,
                textcoords="offset points", xytext=(8, -4))
ax.set_xlabel("Projectile mass (amu)")
ax.set_ylabel("Sputter yield (atoms/ion)")
ax.set_title("Sputter yield vs projectile mass at fixed 100 keV")
plt.savefig(os.path.join(OUT, "2_yield_vs_mass.png")); plt.close()

# ---- Fig 3: ejected-atom energy spectra -------------------------------------
fig, ax = plt.subplots(figsize=(7.0, 4.3))
bins = np.linspace(0, 60, 60)
for r in rows:
    if len(r["E"]):
        ax.hist(r["E"], bins=bins, histtype="step", lw=1.6, color=r["col"],
                label=f"{r['lab']} (Y = {r['Y']:.3f})", density=True)
ax.set_yscale("log")
ax.set_xlabel("Ejection energy (eV)")
ax.set_ylabel("Normalised counts")
ax.set_title("Energy spectrum of sputtered atoms (100 keV, GaAs)")
ax.legend(fontsize=8)
plt.savefig(os.path.join(OUT, "3_energy_spectra.png")); plt.close()

# ---- Fig 4: range and damage vs projectile ----------------------------------
fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8))
axes[0].bar([r["lab"] for r in rows], [r["rp"] for r in rows],
            color=[r["col"] for r in rows])
axes[0].set_yscale("log"); axes[0].set_ylabel("Projected range $R_p$ (Å)")
axes[0].set_title("Range at 100 keV")
for i, r in enumerate(rows):
    axes[0].text(i, r["rp"] * 1.15, f"{r['rp']:.0f} Å", ha="center", fontsize=8.5)
axes[1].bar([r["lab"] for r in rows], [r["vac"] for r in rows],
            color=[r["col"] for r in rows])
axes[1].set_ylabel("Vacancies per ion"); axes[1].set_title("Lattice damage at 100 keV")
for i, r in enumerate(rows):
    axes[1].text(i, r["vac"] * 1.02, f"{r['vac']:.0f}", ha="center", fontsize=8.5)
plt.savefig(os.path.join(OUT, "4_range_damage.png")); plt.close()

print("\nfigures ->", OUT)
