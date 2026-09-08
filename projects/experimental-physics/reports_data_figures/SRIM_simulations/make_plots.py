#!/usr/bin/env python
r"""
make_plots.py -- Generate every reasonable plot from the harvested SRIM/TRIM
outputs in results\<run>\.  Deliberately over-generates; we narrow down to the
report figures later.

Per-run plots  -> results\plots\<run>\*.png
Cross-run plots -> results\plots\_summary\*.png

Pure parsing of SRIM's fixed-format .txt tables (mantissa.E+exp tokens parse
straight to float).  Std numpy + matplotlib only.
"""
import os, csv, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.abspath(__file__))
RES  = os.path.join(ROOT, "results")
OUT  = os.path.join(RES, "plots")

# ---- run metadata (ion, energy label, target, target element names) ----------
RUNS = {
    "01_H_10MeV_Be":    dict(ion="H",  E="10 MeV",  keV=10000, target="Be",   elems=["Be"]),
    "02_H_2MeV_Si":     dict(ion="H",  E="2 MeV",   keV=2000,  target="Si",   elems=["Si"]),
    "03_B_200keV_Si":   dict(ion="B",  E="200 keV", keV=200,   target="Si",   elems=["Si"]),
    "04_H_200keV_Si":   dict(ion="H",  E="200 keV", keV=200,   target="Si",   elems=["Si"]),
    "05_Ar_1keV_GaAs":  dict(ion="Ar", E="1 keV",   keV=1,     target="GaAs", elems=["Ga","As"]),
    "06_Ar_3keV_GaAs":  dict(ion="Ar", E="3 keV",   keV=3,     target="GaAs", elems=["Ga","As"]),
    "07_Ar_10keV_GaAs": dict(ion="Ar", E="10 keV",  keV=10,    target="GaAs", elems=["Ga","As"]),
    "08_Ar_30keV_GaAs": dict(ion="Ar", E="30 keV",  keV=30,    target="GaAs", elems=["Ga","As"]),
    "09_Ar_100keV_GaAs":dict(ion="Ar", E="100 keV", keV=100,   target="GaAs", elems=["Ga","As"]),
    "10_Ar_300keV_GaAs":dict(ion="Ar", E="300 keV", keV=300,   target="GaAs", elems=["Ga","As"]),
    "11_Ar_1MeV_GaAs":  dict(ion="Ar", E="1 MeV",   keV=1000,  target="GaAs", elems=["Ga","As"]),
}
GAAS_RUNS = [k for k, v in RUNS.items() if v["target"] == "GaAs"]

plt.rcParams.update({"figure.dpi": 130, "savefig.dpi": 150, "font.size": 10,
                     "axes.grid": True, "grid.alpha": 0.3, "figure.autolayout": True})


# ---- generic parsers ---------------------------------------------------------
def _lines(path):
    with open(path, "r", encoding="latin-1") as f:
        return f.readlines()

def read_table(path, ncols=None):
    """Return numeric data rows (any line whose whitespace tokens ALL parse as
    float and number >= 2, optionally == ncols)."""
    rows = []
    for ln in _lines(path):
        toks = ln.split()
        if len(toks) < 2:
            continue
        try:
            vals = [float(t) for t in toks]
        except ValueError:
            continue
        if ncols is not None and len(vals) != ncols:
            continue
        rows.append(vals)
    return np.array(rows) if rows else np.empty((0, 0))

def read_ragged(path):
    """Like read_table but keeps the widest consistent numeric block."""
    blocks = {}
    for ln in _lines(path):
        toks = ln.split()
        if len(toks) < 2:
            continue
        try:
            vals = [float(t) for t in toks]
        except ValueError:
            continue
        blocks.setdefault(len(vals), []).append(vals)
    if not blocks:
        return np.empty((0, 0))
    n = max(blocks, key=lambda k: len(blocks[k]))   # most common width
    return np.array(blocks[n])


def savefig(run, name):
    d = os.path.join(OUT, run)
    os.makedirs(d, exist_ok=True)
    p = os.path.join(d, name + ".png")
    plt.savefig(p)
    plt.close()
    return p


# ---- per-run plotting --------------------------------------------------------
def to_um(a):  # Angstrom -> micron
    return a / 1e4

def plot_run(run):
    meta = RUNS[run]
    d = os.path.join(RES, run)
    ion, E, elems = meta["ion"], meta["E"], meta["elems"]
    tag = f"{ion} {E} → {meta['target']}"
    made = []

    # 1) RANGE: implanted-ion distribution + target recoil distribution --------
    p = os.path.join(d, "RANGE.txt")
    if os.path.exists(p):
        t = read_ragged(p)
        if t.size:
            depth = to_um(t[:, 0])
            # ion implant profile
            plt.figure(figsize=(7, 4.3))
            plt.plot(depth, t[:, 1], color="C0")
            plt.fill_between(depth, t[:, 1], alpha=0.25, color="C0")
            plt.xlabel("Depth (µm)"); plt.ylabel("Ion distribution\n(atoms/cm³)/(atoms/cm²)")
            plt.title(f"{tag}: implanted-ion depth distribution")
            made.append(savefig(run, "01_range_ion_profile"))
            # target recoil distribution(s)
            if t.shape[1] >= 3:
                plt.figure(figsize=(7, 4.3))
                for j in range(2, t.shape[1]):
                    lab = elems[j - 2] if (j - 2) < len(elems) else f"target {j-1}"
                    plt.plot(depth, t[:, j], label=f"{lab} recoils")
                plt.xlabel("Depth (µm)"); plt.ylabel("Recoil-atom distribution")
                plt.title(f"{tag}: final target-recoil distribution")
                plt.legend()
                made.append(savefig(run, "02_range_recoil_profile"))

    # 2) VACANCY: damage profile (total + per element) + knock-ons -------------
    p = os.path.join(d, "VACANCY.txt")
    if os.path.exists(p):
        t = read_ragged(p)
        if t.size and t.shape[1] >= 3:
            depth = to_um(t[:, 0])
            knock = t[:, 1]
            vac_cols = t[:, 2:]
            total = vac_cols.sum(axis=1)
            plt.figure(figsize=(7, 4.3))
            plt.plot(depth, total, color="k", lw=1.8, label="total vacancies")
            if vac_cols.shape[1] > 1:
                for j in range(vac_cols.shape[1]):
                    lab = elems[j] if j < len(elems) else f"elem {j+1}"
                    plt.plot(depth, vac_cols[:, j], lw=1.0, alpha=0.9, label=f"{lab} vacancies")
            plt.plot(depth, knock, "--", color="C3", alpha=0.7, label="ion knock-ons")
            # mark peak-damage depth
            ipk = int(np.argmax(total))
            plt.axvline(depth[ipk], color="C2", ls=":", alpha=0.8)
            plt.annotate(f"peak @ {depth[ipk]:.3g} µm", (depth[ipk], total[ipk]),
                         textcoords="offset points", xytext=(6, -4), fontsize=8, color="C2")
            plt.xlabel("Depth (µm)"); plt.ylabel("Vacancies / (Å·ion)")
            plt.title(f"{tag}: damage (vacancy) depth profile")
            plt.legend(fontsize=8)
            made.append(savefig(run, "03_vacancy_profile"))

    # 3) IONIZ: electronic energy loss (dose) by ions & recoils ----------------
    p = os.path.join(d, "IONIZ.txt")
    if os.path.exists(p):
        t = read_table(p, ncols=3)
        if t.size:
            depth = to_um(t[:, 0])
            plt.figure(figsize=(7, 4.3))
            plt.plot(depth, t[:, 1], label="by ions", color="C0")
            plt.plot(depth, t[:, 2], label="by recoils", color="C1")
            plt.plot(depth, t[:, 1] + t[:, 2], label="total", color="k", lw=1.6)
            ipk = int(np.argmax(t[:, 1] + t[:, 2]))
            plt.axvline(depth[ipk], color="C2", ls=":", alpha=0.7)
            plt.annotate(f"Bragg peak @ {depth[ipk]:.3g} µm",
                         (depth[ipk], (t[:, 1] + t[:, 2])[ipk]),
                         textcoords="offset points", xytext=(6, -4), fontsize=8, color="C2")
            plt.xlabel("Depth (µm)"); plt.ylabel("Ionization  eV/(Å·ion)")
            plt.title(f"{tag}: ionization energy loss vs depth")
            plt.legend()
            made.append(savefig(run, "04_ionization_profile"))

    # 4) PHONON --------------------------------------------------------------
    p = os.path.join(d, "PHONON.txt")
    if os.path.exists(p):
        t = read_table(p, ncols=3)
        if t.size:
            depth = to_um(t[:, 0])
            plt.figure(figsize=(7, 4.3))
            plt.plot(depth, t[:, 1], label="by ions")
            plt.plot(depth, t[:, 2], label="by recoils")
            plt.plot(depth, t[:, 1] + t[:, 2], label="total", color="k", lw=1.6)
            plt.xlabel("Depth (µm)"); plt.ylabel("Phonons / (Å·ion)")
            plt.title(f"{tag}: phonon production vs depth")
            plt.legend()
            made.append(savefig(run, "05_phonon_profile"))

    # 5) E2RECOIL ------------------------------------------------------------
    p = os.path.join(d, "E2RECOIL.txt")
    if os.path.exists(p):
        t = read_ragged(p)
        if t.size and t.shape[1] >= 2:
            depth = to_um(t[:, 0])
            plt.figure(figsize=(7, 4.3))
            plt.plot(depth, t[:, 1], label="energy from ions")
            for j in range(2, t.shape[1]):
                lab = elems[j - 2] if (j - 2) < len(elems) else f"target {j-1}"
                plt.plot(depth, t[:, j], label=f"absorbed by {lab}", alpha=0.8)
            plt.xlabel("Depth (µm)"); plt.ylabel("Energy  eV/(Å·ion)")
            plt.title(f"{tag}: energy transferred to recoils vs depth")
            plt.legend(fontsize=8)
            made.append(savefig(run, "06_e2recoil_profile"))

    # 6) NOVAC replacement collisions ---------------------------------------
    p = os.path.join(d, "NOVAC.txt")
    if os.path.exists(p):
        t = read_table(p, ncols=2)
        if t.size:
            depth = to_um(t[:, 0])
            plt.figure(figsize=(7, 4.3))
            plt.plot(depth, t[:, 1], color="C4")
            plt.xlabel("Depth (µm)"); plt.ylabel("Replacement collisions / (Å·ion)")
            plt.title(f"{tag}: replacement collisions vs depth")
            made.append(savefig(run, "07_replacement_collisions"))

    # 7) LATERAL spread vs depth --------------------------------------------
    # NOTE: SRIM writes 0.0 in LATERAL.txt for depth bins where it has no ion
    # statistics -- those are "no data", NOT a physical zero lateral spread.
    # Plotting through them draws a line along y=0 with single-ion spikes, which
    # is misleading, so mask them and show only populated bins as markers.
    p = os.path.join(d, "LATERAL.txt")
    if os.path.exists(p):
        t = read_table(p, ncols=5)
        if t.size:
            m = t[:, 1] != 0                     # populated bins only
            npop, ntot = int(m.sum()), len(m)
            if npop >= 5:
                depth = to_um(t[m, 0])
                proj, projs = to_um(t[m, 1]), to_um(t[m, 2])
                rad, rads = to_um(t[m, 3]), to_um(t[m, 4])
                plt.figure(figsize=(7, 4.3))
                plt.errorbar(depth, proj, yerr=projs, fmt="o-", ms=3, lw=1,
                             capsize=2, label="projected lateral range", color="C0")
                plt.errorbar(depth, rad, yerr=rads, fmt="s-", ms=3, lw=1,
                             capsize=2, label="radial range", color="C1")
                plt.xlabel("Depth (µm)"); plt.ylabel("Lateral extent (µm)")
                plt.title(f"{tag}: lateral / radial spread vs depth\n"
                          f"({npop}/{ntot} depth bins populated)", fontsize=10)
                plt.legend(fontsize=8)
                made.append(savefig(run, "08_lateral_spread"))

    # 8) Energy partition (ioniz vs phonon, ion vs recoil) ------------------
    pi, pp = os.path.join(d, "IONIZ.txt"), os.path.join(d, "PHONON.txt")
    if os.path.exists(pi) and os.path.exists(pp):
        ti, tp = read_table(pi, 3), read_table(pp, 3)
        if ti.size and tp.size and ti.shape[0] == tp.shape[0]:
            depth = to_um(ti[:, 0])
            plt.figure(figsize=(7, 4.3))
            plt.plot(depth, ti[:, 1], label="ionization (ion)", color="C0")
            plt.plot(depth, ti[:, 2], label="ionization (recoil)", color="C0", ls="--")
            plt.plot(depth, tp[:, 1], label="phonon (ion)", color="C3")
            plt.plot(depth, tp[:, 2], label="phonon (recoil)", color="C3", ls="--")
            plt.xlabel("Depth (µm)"); plt.ylabel("Energy loss  eV/(Å·ion)")
            plt.title(f"{tag}: energy-loss partition vs depth")
            plt.legend(fontsize=8)
            made.append(savefig(run, "09_energy_partition"))

    # 9) RANGE_3D: final-position histogram + stopping cloud ----------------
    p = os.path.join(d, "RANGE_3D.txt")
    if os.path.exists(p):
        t = read_table(p, ncols=4)
        if t.size:
            X = t[:, 1]; Y = t[:, 2]; Z = t[:, 3]
            rad = np.sqrt(Y**2 + Z**2)
            Xu, radu, Yu, Zu = to_um(X), to_um(rad), to_um(Y), to_um(Z)
            # depth histogram
            plt.figure(figsize=(7, 4.3))
            plt.hist(Xu, bins=80, color="C0", alpha=0.85)
            mu, sd = Xu.mean(), Xu.std()
            plt.axvline(mu, color="k", ls="--", label=f"mean Rp = {mu:.4g} µm")
            plt.axvline(mu + sd, color="C3", ls=":", alpha=0.7, label=f"straggle = {sd:.3g} µm")
            plt.axvline(mu - sd, color="C3", ls=":", alpha=0.7)
            plt.xlabel("Final ion depth (µm)"); plt.ylabel("Ions")
            plt.title(f"{tag}: stopped-ion depth histogram (N={len(Xu)})")
            plt.legend(fontsize=8)
            made.append(savefig(run, "10_range3d_depth_hist"))
            # stopping cloud (subsample)
            n = len(Xu); idx = np.random.default_rng(0).choice(n, size=min(n, 5000), replace=False)
            plt.figure(figsize=(7, 4.6))
            plt.scatter(Xu[idx], Yu[idx], s=3, alpha=0.25, color="C0", edgecolors="none")
            plt.xlabel("Depth  X (µm)"); plt.ylabel("Lateral  Y (µm)")
            plt.title(f"{tag}: ion stopping cloud (X–Y projection)")
            plt.axis("equal")
            made.append(savefig(run, "11_range3d_cloud"))
            # radial histogram at stop
            plt.figure(figsize=(7, 4.3))
            plt.hist(radu, bins=80, color="C1", alpha=0.85)
            plt.axvline(radu.mean(), color="k", ls="--", label=f"mean radial = {radu.mean():.3g} µm")
            plt.xlabel("Radial displacement at stop (µm)"); plt.ylabel("Ions")
            plt.title(f"{tag}: radial straggle distribution")
            plt.legend(fontsize=8)
            made.append(savefig(run, "12_range3d_radial_hist"))

    # 10) SPUTTER (GaAs only): energy spectrum, angular dist, counts ---------
    p = os.path.join(d, "SPUTTER.txt")
    if os.path.exists(p):
        typ, Z, energy, cosX = [], [], [], []
        for ln in _lines(p):
            toks = ln.split()
            if len(toks) == 10 and toks[0] in ("S", "B", "T"):
                try:
                    Zv = int(float(toks[2])); Ev = float(toks[3]); cx = float(toks[7])
                except ValueError:
                    continue
                typ.append(toks[0]); Z.append(Zv); energy.append(Ev); cosX.append(cx)
        typ = np.array(typ); Z = np.array(Z); energy = np.array(energy); cosX = np.array(cosX)
        sput = typ == "S"
        if sput.any():
            Zs, Es, cxs = Z[sput], energy[sput], cosX[sput]
            isGa, isAs = Zs == 31, Zs == 33
            # energy spectrum (clip to the informative low-energy Thompson peak,
            # log-y to show the high-energy tail too)
            plt.figure(figsize=(7, 4.3))
            hi = np.percentile(Es, 97) if Es.size else 1
            hi = min(max(hi, 10), 80)          # sputtered atoms are near-surface, few eV
            bins = np.linspace(0, hi, 60)
            plt.hist(Es[isGa], bins=bins, alpha=0.6, label=f"Ga (n={isGa.sum()})", color="C0")
            plt.hist(Es[isAs], bins=bins, alpha=0.6, label=f"As (n={isAs.sum()})", color="C1")
            plt.yscale("log")
            plt.xlabel(f"Ejection energy (eV, shown to {hi:.0f} eV)"); plt.ylabel("Sputtered atoms")
            plt.title(f"{tag}: sputtered-atom energy spectrum")
            plt.legend(fontsize=8)
            made.append(savefig(run, "13_sputter_energy_spectrum"))
            # angular distribution (polar angle from surface normal)
            theta = np.degrees(np.arccos(np.clip(-cxs, -1, 1)))  # cosX<0 = leaving
            plt.figure(figsize=(7, 4.3))
            plt.hist(theta[isGa], bins=45, range=(0, 90), alpha=0.6, label="Ga", color="C0")
            plt.hist(theta[isAs], bins=45, range=(0, 90), alpha=0.6, label="As", color="C1")
            plt.xlabel("Ejection polar angle from surface normal (°)"); plt.ylabel("Sputtered atoms")
            plt.title(f"{tag}: sputtered-atom angular distribution")
            plt.legend(fontsize=8)
            made.append(savefig(run, "14_sputter_angular"))

    return made


# ---- cross-run summary plots -------------------------------------------------
def read_summary():
    rows = {}
    with open(os.path.join(RES, "SUMMARY.csv"), newline="") as f:
        for r in csv.DictReader(f):
            rows[r["job"]] = r
    return rows

def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan

def summary_plots():
    s = read_summary()
    made = []
    order = GAAS_RUNS
    E   = np.array([RUNS[k]["keV"] for k in order], float)
    Y   = np.array([fnum(s[k]["sputter_yield"]) for k in order])
    YGa = np.array([fnum(s[k]["Y_Ga"]) for k in order])
    YAs = np.array([fnum(s[k]["Y_As"]) for k in order])
    N   = np.array([fnum(s[k]["ions"]) for k in order])
    rng = np.array([fnum(s[k]["range_A"]) for k in order])
    bsc = np.array([fnum(s[k]["backscatt"]) for k in order])
    vac = np.array([fnum(s[k]["vac_per_ion"]) for k in order])

    def save(name):
        d = os.path.join(OUT, "_summary"); os.makedirs(d, exist_ok=True)
        p = os.path.join(d, name + ".png"); plt.savefig(p); plt.close(); return p

    # S1 yield vs energy (with realistic ~2/sqrt(N) error bars)
    err = Y * 2.0 / np.sqrt(N)
    plt.figure(figsize=(7.2, 4.6))
    plt.errorbar(E, Y, yerr=err, marker="o", capsize=3, label="total", color="k")
    plt.plot(E, YGa, "s--", label="Ga", color="C0")
    plt.plot(E, YAs, "^--", label="As", color="C1")
    ipk = int(np.nanargmax(Y))
    plt.axvline(E[ipk], color="C2", ls=":", alpha=0.7)
    plt.annotate(f"peak {Y[ipk]:.2f} @ {order[ipk].split('_')[2]}",
                 (E[ipk], Y[ipk]), textcoords="offset points", xytext=(6, 8), color="C2", fontsize=9)
    plt.xscale("log"); plt.xlabel("Ar energy (keV)"); plt.ylabel("Sputter yield (atoms/ion)")
    plt.title("Ar → GaAs: sputter yield vs beam energy")
    plt.legend(); made.append(save("S1_yield_vs_energy"))

    # S2 preferential sputtering ratio
    plt.figure(figsize=(7.2, 4.3))
    plt.plot(E, YAs / YGa, "o-", color="C3")
    plt.axhline(1.0, color="k", ls=":", alpha=0.5)
    plt.xscale("log"); plt.xlabel("Ar energy (keV)"); plt.ylabel("As/Ga yield ratio")
    plt.title("Ar → GaAs: preferential sputtering (As/Ga)")
    made.append(save("S2_preferential_ratio"))

    # S3 range vs energy
    plt.figure(figsize=(7.2, 4.3))
    plt.loglog(E, to_um(rng), "o-", color="C0")
    plt.xlabel("Ar energy (keV)"); plt.ylabel("Projected range Rp (µm)")
    plt.title("Ar → GaAs: projected range vs energy")
    made.append(save("S3_range_vs_energy"))

    # S4 backscatter fraction vs energy
    plt.figure(figsize=(7.2, 4.3))
    plt.semilogx(E, 100 * bsc / N, "o-", color="C4")
    plt.xlabel("Ar energy (keV)"); plt.ylabel("Backscattered ions (%)")
    plt.title("Ar → GaAs: ion backscattering vs energy")
    made.append(save("S4_backscatter_vs_energy"))

    # S5 vacancies/ion vs energy
    plt.figure(figsize=(7.2, 4.3))
    plt.loglog(E, vac, "o-", color="C5")
    plt.xlabel("Ar energy (keV)"); plt.ylabel("Vacancies / ion")
    plt.title("Ar → GaAs: lattice damage vs energy")
    made.append(save("S5_vacancies_vs_energy"))

    # S6/S7 Project 4: B vs H at 200 keV in Si -- overlays
    def load_range_ion(run):
        t = read_ragged(os.path.join(RES, run, "RANGE.txt"))
        return to_um(t[:, 0]), t[:, 1]
    def load_vac_total(run):
        t = read_ragged(os.path.join(RES, run, "VACANCY.txt"))
        return to_um(t[:, 0]), t[:, 2:].sum(axis=1)

    xb, yb = load_range_ion("03_B_200keV_Si"); xh, yh = load_range_ion("04_H_200keV_Si")
    plt.figure(figsize=(7.2, 4.4))
    plt.plot(xb, yb, label="B 200 keV", color="C3")
    plt.plot(xh, yh, label="H 200 keV", color="C0")
    plt.xlabel("Depth (µm)"); plt.ylabel("Implanted-ion distribution")
    plt.title("200 keV implant in Si: B vs H depth distribution")
    plt.legend(); made.append(save("S6_B_vs_H_range"))

    xb, yb = load_vac_total("03_B_200keV_Si"); xh, yh = load_vac_total("04_H_200keV_Si")
    plt.figure(figsize=(7.2, 4.4))
    plt.plot(xb, yb, label="B 200 keV", color="C3")
    plt.plot(xh, yh, label="H 200 keV", color="C0")
    plt.xlabel("Depth (µm)"); plt.ylabel("Vacancies / (Å·ion)")
    plt.title("200 keV implant in Si: B vs H damage profile")
    plt.legend(); made.append(save("S7_B_vs_H_vacancy"))

    return made


def main():
    os.makedirs(OUT, exist_ok=True)
    total = 0
    for run in RUNS:
        made = plot_run(run)
        total += len(made)
        print(f"{run:22s} -> {len(made):2d} plots")
    sm = summary_plots()
    total += len(sm)
    print(f"{'_summary':22s} -> {len(sm):2d} plots")
    print(f"\nTOTAL: {total} plots under {OUT}")

if __name__ == "__main__":
    main()
