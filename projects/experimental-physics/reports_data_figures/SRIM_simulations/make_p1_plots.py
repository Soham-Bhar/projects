#!/usr/bin/env python
r"""
make_p1_plots.py -- Project 1 figures: 1 GeV protons in air (run_1, manual GUI run).

Data sources (the ONLY files written by that run):
  run_1\SRIM Outputs\RANGE_3D.txt  -- 5000 final ion positions  -> range & straggling
  run_1\SRIM Outputs\EXYZ.txt      -- 2.3 GB trajectory file, pre-aggregated by awk into
                                      results\P1_H_1GeV_air\{bragg,trails}.csv
NOTE: the root run_1\*.txt summary tables are STALE (10 MeV H -> Be) and are NOT used.
"""
import os, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = r"D:\SRIM_sim"
R3D  = os.path.join(ROOT, "results", "P1_H_1GeV_air", "RANGE_3D.txt")
AGG  = os.path.join(ROOT, "results", "P1_H_1GeV_air")
OUT  = os.path.join(ROOT, "results", "plots", "P1_H_1GeV_air")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"figure.dpi":130,"savefig.dpi":150,"font.size":10,
                     "axes.grid":True,"grid.alpha":0.3,"figure.autolayout":True})
TAG = "H 1 GeV → air"

# ---------- RANGE_3D: final ion positions (Angstrom -> metres) ----------------
rows=[]
for ln in open(R3D, encoding="latin-1"):
    t=ln.split()
    if len(t)==4 and t[0].isdigit():
        try: rows.append([float(x) for x in t])
        except ValueError: pass
a=np.array(rows)
X=a[:,1]/1e10; Y=a[:,2]/1e10; Z=a[:,3]/1e10          # metres
rad=np.hypot(Y,Z)
Rp, dRp = X.mean(), X.std()
rad_mean, rad_std = rad.mean(), rad.std()
lat_std = np.std(np.concatenate([Y,Z]))
N=len(X)

print(f"N ions                     = {N}")
print(f"(a) Projected range Rp     = {Rp:.1f} m   ({Rp/1000:.4f} km)")
print(f"    Longitudinal straggling= {dRp:.2f} m   ({100*dRp/Rp:.3f} % of Rp)")
print(f"(b) Mean radial displ.     = {rad_mean:.2f} m")
print(f"    Radial straggling (sd) = {rad_std:.2f} m")
print(f"    Lateral 1-sigma (Y,Z)  = {lat_std:.2f} m")
print(f"    95th pct radius        = {np.percentile(rad,95):.2f} m")
print(f"(c) Rp + 3*dRp             = {Rp+3*dRp:.1f} m")
print(f"    deepest single ion     = {X.max():.1f} m")
# One ion in 5000 stops far short (a rare hard-scatter outlier). It barely moves the
# mean but inflates the straggling, so report both and keep the full-sample value as
# the primary (conservative) number for the safe-distance question.
core = X[X > Rp - 5*dRp]
n_out = N - len(core)
print(f"    outliers (<Rp-5sigma)  = {n_out}  (shallowest {X.min():.1f} m)")
print(f"    robust Rp / sigma      = {core.mean():.1f} m / {core.std():.2f} m "
      f"({100*core.std()/core.mean():.3f} %) excluding them")

# ---------- Fig 1: stopping-depth histogram ----------------------------------
plt.figure(figsize=(7,4.3))
plt.hist(X[X > Rp - 5*dRp], bins=70, color="C0", alpha=.85)
plt.axvline(Rp, color="k", ls="--", label=f"$R_p$ = {Rp:.0f} m")
for s,lab in ((1,rf"$\pm1\sigma$ = {dRp:.1f} m"),(3,None)):
    plt.axvline(Rp+s*dRp, color="C3", ls=":", alpha=.8, label=lab)
    plt.axvline(Rp-s*dRp, color="C3", ls=":", alpha=.8)
plt.axvline(Rp+3*dRp, color="C2", ls="-.", alpha=.9, label=rf"$R_p+3\sigma$ = {Rp+3*dRp:.0f} m")
lo = min(core.min(), Rp-4*dRp); hi = max(X.max(), Rp+4*dRp)
plt.xlim(lo - 0.02*(hi-lo), hi + 0.02*(hi-lo))
if n_out:
    plt.annotate(f"{n_out} ion outside axis (shallowest {X.min():.0f} m)",
                 xy=(0.02, 0.72), xycoords="axes fraction", fontsize=8, color="0.35")
plt.xlabel("Final ion depth (m)"); plt.ylabel("Ions")
plt.title(f"{TAG}: stopped-ion depth distribution (N={N})")
plt.legend(fontsize=8, loc="upper left"); plt.savefig(os.path.join(OUT,"1_range_hist.png")); plt.close()

# ---------- Fig 2: radial displacement histogram ------------------------------
plt.figure(figsize=(7,4.3))
plt.hist(rad, bins=70, color="C1", alpha=.85)
plt.axvline(rad_mean, color="k", ls="--", label=f"mean radius = {rad_mean:.1f} m")
plt.axvline(np.percentile(rad,95), color="C3", ls=":", label=f"95th pct = {np.percentile(rad,95):.1f} m")
plt.xlabel("Radial displacement at end of range (m)"); plt.ylabel("Ions")
plt.title(f"{TAG}: beam radius at end of range")
plt.legend(fontsize=8); plt.savefig(os.path.join(OUT,"2_radial_hist.png")); plt.close()

# ---------- Fig 3: stopping cloud (X-Y) --------------------------------------
plt.figure(figsize=(7,4.4))
plt.scatter(X, Y, s=3, alpha=.25, color="C0", edgecolors="none")
plt.xlabel("Depth  X (m)"); plt.ylabel("Lateral  Y (m)")
plt.title(f"{TAG}: ion stopping cloud")
plt.savefig(os.path.join(OUT,"3_stopping_cloud.png")); plt.close()

# ---------- Fig 4+5: from the aggregated EXYZ --------------------------------
b=np.genfromtxt(os.path.join(AGG,"bragg.csv"), delimiter=",", names=True)
d=b["bin_m"]; se=b["mean_Se_eV_per_A"]*1e4          # eV/A -> MeV/m
dep=b["dE_keV_total"]/1000.0/5.0/1000.0             # keV over 1000 ions per 5 m -> MeV/m/ion
m=d<=Rp+6*dRp

plt.figure(figsize=(7,4.3))
plt.plot(d[m], dep[m], color="C3", lw=1.4)
ipk=int(np.argmax(dep[m]))
plt.axvline(d[m][ipk], color="C2", ls=":", alpha=.8)
plt.annotate(f"Bragg peak @ {d[m][ipk]:.0f} m\n{dep[m][ipk]:.3f} MeV/m",
             (d[m][ipk], dep[m][ipk]), textcoords="offset points", xytext=(-105,-6),
             fontsize=8, color="C2")
plt.axhline(dep[m][0], color="k", ls="--", alpha=.5, lw=.9)
plt.annotate(f"entrance {dep[m][0]:.3f} MeV/m", (d[m][2], dep[m][0]),
             textcoords="offset points", xytext=(10,-30), fontsize=8)
plt.xlabel("Depth (m)"); plt.ylabel("Energy deposited (MeV/m per ion)")
plt.title(f"{TAG}: energy-deposition (Bragg) curve")
plt.savefig(os.path.join(OUT,"4_bragg_curve.png")); plt.close()
print(f"(d) entrance dose          = {dep[m][0]:.4f} MeV/m ; peak = {dep[m][ipk]:.4f} MeV/m "
      f"at {d[m][ipk]:.0f} m  -> peak/entrance = {dep[m][ipk]/dep[m][0]:.1f}x")

plt.figure(figsize=(7,4.3))
plt.plot(d[m], se[m], color="C0", lw=1.4)
plt.xlabel("Depth (m)"); plt.ylabel("Electronic stopping power (MeV/m)")
plt.title(f"{TAG}: mean electronic stopping power vs depth")
plt.savefig(os.path.join(OUT,"5_stopping_power.png")); plt.close()

# ---------- Fig 6: ion trails -------------------------------------------------
t=np.genfromtxt(os.path.join(AGG,"trails.csv"), delimiter=",")
from matplotlib.collections import LineCollection
segs,cols=[],[]
for i in np.unique(t[:,0]):
    p=t[t[:,0]==i]
    xy=np.column_stack([p[:,1]/1e10, p[:,2]/1e10])
    segs+= [[xy[k],xy[k+1]] for k in range(len(xy)-1)]
    cols+= list(p[:-1,4]/1000.0)                     # keV -> MeV
lc=LineCollection(segs, cmap="plasma", array=np.array(cols), linewidths=.7, alpha=.85)
lc.set_clim(0,1000)
fig,ax=plt.subplots(figsize=(7.6,4.5)); ax.add_collection(lc); ax.autoscale()
ax.set_xlabel("Depth  X (m)"); ax.set_ylabel("Lateral  Y (m)")
ax.set_title(f"{TAG}: ion trails (40 ions, coloured by energy)"); ax.grid(alpha=.25)
fig.colorbar(lc, ax=ax, label="Ion energy (MeV)")
plt.savefig(os.path.join(OUT,"6_ion_trails.png")); plt.close()
print("\nplots ->", OUT)
