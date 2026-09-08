#!/usr/bin/env python
r"""make_p6_plots.py -- Project 6 figures (proton therapy, layered head model).

Reads results\P6_*MeV\{IONIZ,VACANCY,TDATA}.txt -> results\plots\P6\.
Geometry: skin 0-1, skull 1-5, brain 5-8, TUMOUR 8-10 mm (centre 9 mm).
The 39 MeV run has an extra 4 mm brain backing (14 mm total) so it stops inside.
"""
import os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = r"D:\SRIM_sim"
OUT = os.path.join(ROOT, "results", "plots", "P6")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"figure.dpi": 130, "savefig.dpi": 150, "font.size": 10,
                     "axes.grid": True, "grid.alpha": 0.3, "figure.autolayout": True})

BOUND = [0, 1, 5, 8, 10]
LAYNAME = ["skin", "skull", "brain", "tumour"]
TUM = (8.0, 10.0)
CENTRE = 9.0

RUNS = [(33, "P6_33MeV"), (34, "P6_34MeV"), (35, "P6_35MeV"), (35.5, "P6_35.5MeV"),
        (35.65, "P6_35.65MeV"), (36, "P6_36MeV"), (36.5, "P6_36.5MeV"),
        (37, "P6_37MeV"), (38, "P6_38MeV")]
EXT = (39, "P6_39MeV_ext")


def ragged(p):
    b = {}
    for ln in open(p, encoding="latin-1"):
        t = ln.split()
        if len(t) < 2:
            continue
        try:
            v = [float(x) for x in t]
        except ValueError:
            continue
        b.setdefault(len(v), []).append(v)
    n = max(b, key=lambda k: len(b[k]))
    return np.array(b[n])


def dose(job):
    """depth (mm), energy deposited (MeV/mm per ion)"""
    I = ragged(os.path.join(ROOT, "results", job, "IONIZ.txt"))
    return I[:, 0] / 1e7, (I[:, 1] + I[:, 2]) * 10.0


def damage(job):
    V = ragged(os.path.join(ROOT, "results", job, "VACANCY.txt"))
    return V[:, 0] / 1e7, V[:, 2:].sum(axis=1)


def tdata(job, pat):
    t = open(os.path.join(ROOT, "results", job, "TDATA.txt"), encoding="latin-1").read()
    m = re.search(pat, t)
    return m.group(1).strip() if m else None


def shade(ax, xmax=10.0):
    ax.axvspan(TUM[0], TUM[1], color="C3", alpha=0.13, zorder=0)
    ax.axvline(CENTRE, color="C3", ls="--", lw=1.1, zorder=1)
    for b in BOUND[1:-1]:
        ax.axvline(b, color="0.5", ls=":", lw=0.8, zorder=0)
    top = ax.get_ylim()[1]
    for (a, b), n in zip(zip(BOUND[:-1], BOUND[1:]), LAYNAME):
        if a < xmax:
            ax.text((a + min(b, xmax)) / 2, top * 0.96, n, ha="center", va="top",
                    fontsize=7.5, color="0.35")


sel = [33, 35, 36, 37, 38]

# Fig 1 -- depth-dose curves
plt.figure(figsize=(7.2, 4.4))
ax = plt.gca()
for E, j in RUNS:
    if E not in sel:
        continue
    d, y = dose(j)
    lab = f"{E:g} MeV" + (" (exits target)" if E == 38 else "")
    ax.plot(d, y, lw=1.4, label=lab)
ax.set_ylim(0, None)
ax.set_xlim(0, 10)
shade(ax)
ax.set_xlabel("Depth (mm)")
ax.set_ylabel("Energy deposited (MeV/mm per ion)")
ax.set_title("Depth-dose curves: Bragg peak walking through the tumour")
ax.legend(fontsize=8, loc="upper left")
plt.savefig(os.path.join(OUT, "1_depth_dose.png"))
plt.close()

# Fig 2 -- damage profiles
plt.figure(figsize=(7.2, 4.4))
ax = plt.gca()
for E, j in RUNS:
    if E not in sel:
        continue
    d, y = damage(j)
    ax.plot(d, y, lw=1.4, label=f"{E:g} MeV")
ax.set_ylim(0, None)
ax.set_xlim(0, 10)
shade(ax)
ax.set_xlabel("Depth (mm)")
ax.set_ylabel("Vacancies / (Angstrom.ion)")
ax.set_title("Damage (recoil-atom) profiles vs depth")
ax.legend(fontsize=8, loc="upper left")
plt.savefig(os.path.join(OUT, "2_damage.png"))
plt.close()

# Fig 3 -- 36 MeV detail with dose ratio
d, y = dose("P6_36MeV")
i9 = int(np.argmin(abs(d - CENTRE)))
trav = y[(d > 0) & (d < 8.0)].mean()
plt.figure(figsize=(7.2, 4.4))
ax = plt.gca()
ax.plot(d, y, color="C0", lw=1.7, label="36 MeV")
ax.fill_between(d, y, where=(d >= TUM[0]) & (d <= TUM[1]), color="C3", alpha=0.35,
                label="tumour 8-10 mm")
ax.axhline(trav, color="0.4", ls="-.", lw=1, label=f"mean traversed 0-8 mm = {trav:.2f}")
ax.plot(CENTRE, y[i9], "o", color="C3", ms=7)
ax.set_ylim(0, None)
ax.set_xlim(0, 10)
ax.annotate(f"tumour centre {y[i9]:.2f} MeV/mm\nratio {y[i9]/trav:.2f}x",
            (CENTRE, y[i9]), textcoords="offset points", xytext=(-112, -34),
            fontsize=8.5, color="C3")
shade(ax)
ax.set_xlabel("Depth (mm)")
ax.set_ylabel("Energy deposited (MeV/mm per ion)")
ax.set_title("Optimum energy: 36 MeV places the Bragg peak on the tumour")
ax.legend(fontsize=8, loc="upper left")
plt.savefig(os.path.join(OUT, "3_optimum_36MeV.png"))
plt.close()

# Fig 4 -- dose at tumour centre vs energy
Ex, Dx, Rx = [], [], []
for e, j in RUNS:
    if int(tdata(j, r'Total Transmitted Ions  =\s*([0-9]+)')) > 50:
        continue
    d, y = dose(j)
    k = int(np.argmin(abs(d - CENTRE)))
    Ex.append(e)
    Dx.append(y[k])
    Rx.append(y[k] / y[(d > 0) & (d < 8.0)].mean())
plt.figure(figsize=(7.2, 4.3))
plt.plot(Ex, Dx, "o-", color="C0", label="dose at tumour centre")
if 35.65 in Ex:
    _k = Ex.index(35.65)
    plt.plot(35.65, Dx[_k], "D", color="C3", ms=8, zorder=5,
             label="35.65 MeV: fitted $R_p$ = 9.00 mm")
kk = int(np.argmax(Dx))
plt.axvline(Ex[kk], color="C2", ls=":")
plt.annotate(f"optimum {Ex[kk]:g} MeV\n{Dx[kk]:.2f} MeV/mm (ratio {Rx[kk]:.2f}x)",
             (Ex[kk], Dx[kk]), textcoords="offset points", xytext=(12, -26),
             fontsize=9, color="C2")
plt.xlabel("Proton energy (MeV)")
plt.ylabel("Dose at tumour centre (MeV/mm per ion)")
plt.title("Dose delivered to the tumour centre vs beam energy")
plt.legend(fontsize=8)
plt.savefig(os.path.join(OUT, "4_dose_vs_energy.png"))
plt.close()

# Fig 5 -- range-energy calibration
Ee, Rp = [], []
for e, j in RUNS + [EXT]:
    if int(tdata(j, r'Total Transmitted Ions  =\s*([0-9]+)')) > 50:
        continue
    Ee.append(e)
    Rp.append(float(tdata(j, r'Average Range\s+=\s+([0-9.E+-]+) Angstroms')) / 1e7)
Ee = np.array(Ee, float)
Rp = np.array(Rp, float)
m, c = np.polyfit(Ee, Rp, 1)
Eopt = (CENTRE - c) / m
plt.figure(figsize=(7.2, 4.3))
plt.plot(Ee, Rp, "o", color="C0", label="simulated $R_p$")
xs = np.linspace(Ee.min() - 0.5, Ee.max() + 0.5, 50)
plt.plot(xs, m * xs + c, "-", color="C1", lw=1.2, label=f"fit: {m*1000:.0f} um/MeV")
plt.axhline(CENTRE, color="C3", ls="--", lw=1)
plt.axvline(Eopt, color="C2", ls=":")
plt.annotate(f"$R_p$ = 9.00 mm at {Eopt:.2f} MeV", (Eopt, CENTRE),
             textcoords="offset points", xytext=(-150, 12), fontsize=9, color="C2")
plt.xlabel("Proton energy (MeV)")
plt.ylabel("Projected range $R_p$ (mm)")
plt.title("Range-energy calibration for the layered head model")
plt.legend(fontsize=8)
plt.savefig(os.path.join(OUT, "5_range_energy.png"))
plt.close()

# Fig 6 -- 39 MeV pass-through on the extended target
d, y = dose(EXT[1])
fig, ax = plt.subplots(figsize=(7.2, 4.4))
ax.plot(d, y, color="C0", lw=1.6, label="39 MeV dose")
ax.axvspan(TUM[0], TUM[1], color="C3", alpha=0.13)
ax.axvline(CENTRE, color="C3", ls="--", lw=1.1)
for b in BOUND[1:] + [14.0]:
    ax.axvline(b, color="0.5", ls=":", lw=0.8)
ax.set_xlim(0, 14)
ax.set_ylim(0, None)
ipk = int(np.argmax(y))
ax.plot(d[ipk], y[ipk], "o", color="C2", ms=6)
ax.annotate(f"Bragg peak {d[ipk]:.2f} mm\nbeyond the tumour",
            (d[ipk], y[ipk]), textcoords="offset points", xytext=(-42, -46),
            fontsize=8.5, color="C2")
top = ax.get_ylim()[1]
ax.text(9.0, top * 0.93, "tumour", ha="center", fontsize=7.5, color="0.35")
ax.text(12.0, top * 0.93, "brain backing\n(added)", ha="center", fontsize=7.5, color="0.35")
ax.set_xlabel("Depth (mm)")
ax.set_ylabel("Energy deposited (MeV/mm per ion)")
ax.set_title("39 MeV: protons pass straight through the tumour (14 mm target)")
ax.legend(fontsize=8, loc="upper left")
plt.savefig(os.path.join(OUT, "6_39MeV_passthrough.png"))
plt.close()

d, y = dose("P6_36MeV")
i9 = int(np.argmin(abs(d - CENTRE)))
tr = y[(d > 0) & (d < 8.0)].mean()
print(f"36 MeV : centre={y[i9]:.3f}  entrance={y[0]:.3f}  meanTraversed={tr:.3f} MeV/mm")
print(f"         ratio centre/traversed = {y[i9]/tr:.2f}")
print(f"         ratio centre/entrance  = {y[i9]/y[0]:.2f}")
print(f"fit: {m*1000:.0f} um/MeV,  Rp = 9.00 mm at {Eopt:.3f} MeV")
d39, y39 = dose(EXT[1])
print(f"39 MeV : Bragg peak at {d39[int(np.argmax(y39))]:.2f} mm (tumour ends at 10 mm)")
print("figures ->", OUT)
