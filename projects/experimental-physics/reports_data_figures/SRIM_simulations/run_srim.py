#!/usr/bin/env python3
r"""
run_srim.py -- Configure TRIM.IN files and run SRIM/TRIM 2013 Monte-Carlo
ion-transport simulations in parallel across cloned SRIM root directories.

Why cloned roots: TRIM writes its output .txt files (RANGE.txt, SPUTTER.txt, ...)
to *hardcoded* names in the SRIM root. Two concurrent runs in the same root would
clobber each other, so each concurrent job runs in its own copy (run_1 .. run_6),
and results are harvested into a clean per-job folder afterward.

Batch mode: TRIMAUTO line1 = "1" makes TRIM.exe run from TRIM.IN with no keyboard
input and self-terminate after all ions. We launch TRIM.exe (NOT SRIM.exe, which is
the GUI menu) with cwd set to the run directory.

Usage:
    python run_srim.py            # run the full 11-job queue over 6 workers
    python run_srim.py --write    # only (re)write TRIM.IN files, do not execute
    python run_srim.py --jobs 5   # run only the first N jobs of the queue
"""

import argparse
import ctypes
import os
import queue
import shutil
import subprocess
import sys
import threading
import time
from ctypes import wintypes

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
SRIM_ROOT   = r"D:\SRIM_sim"
RUN_DIRS    = [os.path.join(SRIM_ROOT, f"run_{i}") for i in range(1, 7)]  # 6 workers
RESULTS_DIR = os.path.join(SRIM_ROOT, "results")
TRIM_EXE    = "TRIM.exe"

# Known TRIM output files to harvest (hardcoded names TRIM rewrites each run).
OUTPUT_FILES = [
    "TRIM.IN", "TDATA.txt", "RANGE.txt", "RANGE_3D.txt", "VACANCY.txt",
    "IONIZ.txt", "PHONON.txt", "E2RECOIL.txt", "NOVAC.txt", "LATERAL.txt",
    "COLLISON.txt", "BACKSCAT.txt", "TRANSMIT.txt", "SPUTTER.txt", "EXYZ.txt",
]

# --------------------------------------------------------------------------- #
# Element data:  symbol -> (Z, atomic mass amu of dominant isotope, SRIM default)
# --------------------------------------------------------------------------- #
ELEMENTS = {
    "H":  (1,  1.008),
    "B":  (5,  11.009),
    "Ar": (18, 39.962),
    "Be": (4,  9.012),
    "Si": (14, 28.086),
    "Ga": (31, 69.723),
    "As": (33, 74.922),
    "Au": (79, 196.967),   # Project 5c projectile
}

# Per-element target properties: (displacement Ed, lattice El, surface Es) in eV.
# Surface binding energy (Es) is the key driver of sputter yield -> use SRIM DB values.
BINDING = {
    "Be": (25.0, 3.0, 3.38),
    "Si": (15.0, 2.0, 4.70),
    "Ga": (25.0, 3.0, 2.82),
    "As": (25.0, 3.0, 1.26),
}


# --------------------------------------------------------------------------- #
# TRIM.IN builder
# --------------------------------------------------------------------------- #
def _fmt_stoich(v):
    if abs(v - 1.0) < 1e-9:
        return "1"
    if abs(v) < 1e-9:
        return "0"
    return f"{v:.6f}"


def build_trimin(cfg):
    """Return the full text of a TRIM.IN file for one job.

    cfg keys:
        ion       : ion symbol (e.g. "H")
        energy    : energy in keV
        ions      : total ions to simulate (the line-3 "Number" field)
        autosave  : autosave interval (>= ions => no intermediate save)
        cascade   : 2 = Full Damage Cascades (needed for damage + sputtering)
        desc      : quoted target description string
        elements  : list of target element symbols, in column order
        layers    : list of (name, width_Ang, density_gcc, {sym: stoich})
        disk      : dict of diskfile flags {ranges, backscat, transmit,
                    sputter, collisions, exyz}
        plottype  : 5 = no on-screen plots (fastest); txt outputs still written
    """
    Z1, M1 = ELEMENTS[cfg["ion"]]
    els = cfg["elements"]
    n_el = len(els)
    layers = cfg["layers"]
    n_layer = len(layers)
    d = cfg["disk"]

    lines = []
    lines.append("==> SRIM-2011.00 This file controls TRIM Calculations.")
    lines.append("Ion: Z1 ,  M1,  Energy (keV), Angle,Number,Bragg Corr,AutoSave Number.")
    # Z1, M1, Energy, Angle, Number(=total ions), Bragg Corr, AutoSave  (free format)
    lines.append(f"     {Z1}   {M1:g}   {cfg['energy']:g}   0   {cfg['ions']}   0   {cfg['autosave']}")
    lines.append("Cascade Type[1=None;2=Full;3=Mono/Sputt;4-5=Ions(TRIM.dat);6-7=Recoils(TRIM.dat)], Random Number Seed, Reminders")
    lines.append(f"                      {cfg['cascade']}                                   0       0")
    lines.append("Diskfiles (0=no,1=yes): Ranges, Backscatt, Transmit, Sputtered, Collisions(1=Ion;2=Ion+Recoils), Special EXYZ.txt file")
    lines.append(f"                          {d['ranges']}       {d['backscat']}           "
                 f"{d['transmit']}       {d['sputter']}               {d['collisions']}                               {d['exyz']}")
    lines.append("Target material : Number of Elements & Layers")
    lines.append(f'"{cfg["desc"]}"       {n_el}               {n_layer}')
    lines.append("PlotType (0-5); Plot Depths: Xmin, Xmax(Ang.) [=0 0 for Viewing Full Target]")
    lines.append(f"       {cfg['plottype']}                         0            0")
    lines.append("Target Elements:    Z   Mass(amu)")
    for i, sym in enumerate(els, start=1):
        Z, M = ELEMENTS[sym]
        lines.append(f"Atom {i} = {sym} =        {Z}     {M:g}")

    # Column header line for the layer table (element stoich columns).
    hdr_cols = "".join(f"{sym+'('+str(ELEMENTS[sym][0])+')':>8}" for sym in els)
    lines.append("Layer   Layer Name /               Width Density    " + hdr_cols)
    lines.append("Numb.   Description                (Ang) (g/cm3)   " + "  Stoich" * n_el)

    for n, (name, width, dens, stoich) in enumerate(layers, start=1):
        stoich_str = "  ".join(_fmt_stoich(stoich.get(sym, 0.0)) for sym in els)
        # Width must be a plain integer of Angstroms -- never scientific notation,
        # which SRIM's Val() parser can misread (e.g. "1.5e+07").
        lines.append(f' {n}      "{name}"           {int(round(width))}  {dens:g}       {stoich_str}')

    lines.append("0  Target layer phases (0=Solid, 1=Gas)")
    lines.append(" ".join("0" for _ in layers) + " ")
    # Bragg compound-correction is ONE value PER LAYER (not per element). Getting
    # this per-element desyncs SRIM's parser -> VB "Type mismatch" (error 13) on a
    # later line whenever n_elements != n_layers (e.g. a 2-element single layer).
    lines.append("Target Compound Corrections (Bragg)")
    lines.append(" " + "   ".join("1" for _ in layers) + "  ")
    lines.append("Individual target atom displacement energies (eV)")
    lines.append("   " + "   ".join(f"{BINDING[s][0]:g}" for s in els))
    lines.append("Individual target atom lattice binding energies (eV)")
    lines.append("   " + "   ".join(f"{BINDING[s][1]:g}" for s in els))
    lines.append("Individual target atom surface binding energies (eV)")
    lines.append("   " + "   ".join(f"{BINDING[s][2]:g}" for s in els))
    lines.append("Stopping Power Version (1=2008, 0=2008)")
    lines.append(" 0 ")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# Job queue
# --------------------------------------------------------------------------- #
def si_layer(width):
    return ("Silicon", width, 2.321, {"Si": 1.0})


def make_jobs():
    """Return the ordered list of simulation jobs.

    Target thicknesses are chosen a few x the ion range so ions stop inside the
    target and the full end-of-range straggling is captured (never punch through).
    """
    jobs = []

    # --- Run 1: 10 MeV protons into Beryllium ---------------------------------
    # Rp(10 MeV H in Be) ~ 0.9 mm; use 1.5 mm target.
    jobs.append(dict(
        name="01_H_10MeV_Be", ion="H", energy=10000, ions=10000, autosave=10000,
        cascade=2, plottype=5, desc="10 MeV H into Beryllium",
        elements=["Be"],
        layers=[("Beryllium", 15_000_000, 1.848, {"Be": 1.0})],
        disk=dict(ranges=1, backscat=0, transmit=0, sputter=0, collisions=0, exyz=0),
    ))

    # --- Run 2: 2 MeV protons into Silicon (IBIC) -----------------------------
    # Rp(2 MeV H in Si) ~ 48 um; use 80 um target.
    jobs.append(dict(
        name="02_H_2MeV_Si", ion="H", energy=2000, ions=10000, autosave=10000,
        cascade=2, plottype=5, desc="2 MeV H into Silicon",
        elements=["Si"], layers=[si_layer(800_000)],
        disk=dict(ranges=1, backscat=0, transmit=0, sputter=0, collisions=0, exyz=0),
    ))

    # --- Run 3: 200 keV Boron into Silicon (implant doping) -------------------
    # Rp(200 keV B in Si) ~ 0.5 um; use 2 um target.
    jobs.append(dict(
        name="03_B_200keV_Si", ion="B", energy=200, ions=10000, autosave=10000,
        cascade=2, plottype=5, desc="200 keV B into Silicon",
        elements=["Si"], layers=[si_layer(20_000)],
        disk=dict(ranges=1, backscat=0, transmit=0, sputter=0, collisions=0, exyz=0),
    ))

    # --- Run 4: 200 keV protons into Silicon (for comparison with B) ----------
    # Rp(200 keV H in Si) ~ 1.9 um; use 3 um target.
    jobs.append(dict(
        name="04_H_200keV_Si", ion="H", energy=200, ions=10000, autosave=10000,
        cascade=2, plottype=5, desc="200 keV H into Silicon",
        elements=["Si"], layers=[si_layer(30_000)],
        disk=dict(ranges=1, backscat=0, transmit=0, sputter=0, collisions=0, exyz=0),
    ))

    # --- Run 5: Ar -> GaAs sputter-yield energy sweep -------------------------
    # Full cascade + Sputtered diskfile => SPUTTER.txt. Thickness scaled per energy.
    # More ions at low energy (low yield => needs more counts) than at high energy.
    sweep = [
        # energy_keV, thickness_Ang, ions
        (1,     500,  20000),
        (3,     600,  20000),
        (10,    900,  20000),
        (30,    1600, 20000),
        (100,   3200, 15000),
        (300,   6500, 10000),
        (1000, 16000, 8000),
    ]
    for e, thick, nions in sweep:
        etag = f"{e}keV" if e < 1000 else "1MeV"
        jobs.append(dict(
            name=f"{5 + sweep.index((e, thick, nions)):02d}_Ar_{etag}_GaAs",
            ion="Ar", energy=e, ions=nions, autosave=nions, cascade=2, plottype=5,
            desc=f"{etag} Ar into GaAs (sputter)",
            elements=["Ga", "As"],
            layers=[("GaAs", thick, 5.320, {"Ga": 0.5, "As": 0.5})],
            disk=dict(ranges=1, backscat=0, transmit=0, sputter=1, collisions=0, exyz=0),
        ))
    return jobs


# --------------------------------------------------------------------------- #
# Execution
# --------------------------------------------------------------------------- #
_print_lock = threading.Lock()

# Stagger + watchdog tuning. Launching several VB6 TRIM.exe instances in the same
# instant makes some hang at startup (blank window, 0 CPU) due to a shared-init
# race. We stagger initial launches and watchdog each one: a healthy TRIM burns
# CPU within a second or two; if it hasn't after WARMUP, it's wedged -> relaunch.
LAUNCH_STAGGER = 3.0     # seconds between successive worker first-launches
WARMUP         = 30.0    # seconds to wait for a process to prove it is computing
MIN_CPU        = 1.0     # CPU-seconds that count as "it started computing"
MAX_RETRIES    = 3       # relaunch attempts on a startup hang


def log(msg):
    with _print_lock:
        print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _proc_cpu_seconds(pid):
    """Total (kernel+user) CPU seconds a PID has consumed, or None if unavailable."""
    PROCESS_QUERY_LIMITED_INFORMATION = 0x1000
    k32 = ctypes.windll.kernel32
    h = k32.OpenProcess(PROCESS_QUERY_LIMITED_INFORMATION, False, pid)
    if not h:
        h = k32.OpenProcess(0x0400, False, pid)  # PROCESS_QUERY_INFORMATION
    if not h:
        return None
    try:
        ct, et, kt, ut = (wintypes.FILETIME() for _ in range(4))
        if not k32.GetProcessTimes(h, ctypes.byref(ct), ctypes.byref(et),
                                   ctypes.byref(kt), ctypes.byref(ut)):
            return None
        to100 = lambda ft: (ft.dwHighDateTime << 32) | ft.dwLowDateTime
        return (to100(kt) + to100(ut)) / 1e7
    finally:
        k32.CloseHandle(h)


def launch_with_watchdog(exe, run_dir, name):
    """Start TRIM.exe and make sure it is actually computing; relaunch if it hangs.

    Returns the live/finished Popen object of the healthy attempt (or the last
    attempt if all retries hung).
    """
    for attempt in range(1, MAX_RETRIES + 1):
        proc = subprocess.Popen([exe], cwd=run_dir)
        t0 = time.time()
        while time.time() - t0 < WARMUP:
            if proc.poll() is not None:
                return proc  # finished already (very fast job) -> healthy
            cpu = _proc_cpu_seconds(proc.pid)
            if cpu is not None and cpu >= MIN_CPU:
                return proc  # it's computing
            time.sleep(1.0)
        cpu = _proc_cpu_seconds(proc.pid)
        log(f"{name}: startup HANG (attempt {attempt}/{MAX_RETRIES}, "
            f"CPU={cpu}s after {WARMUP:.0f}s) -> killing & relaunching")
        try:
            proc.kill()
            proc.wait(timeout=10)
        except Exception:  # noqa: BLE001
            pass
        time.sleep(2.0)
    return proc  # give up; caller will harvest whatever exists


def set_trimauto(run_dir, mode=1):
    """Write TRIMAUTO so TRIM runs in batch mode (line1=1)."""
    with open(os.path.join(run_dir, "TRIMAUTO"), "w") as f:
        f.write(f"{mode}\n\n")


def clean_outputs(run_dir):
    """Remove stale output files so a job never harvests a previous run's data."""
    for fn in OUTPUT_FILES:
        if fn == "TRIM.IN":
            continue
        for base in (run_dir, os.path.join(run_dir, "SRIM Outputs")):
            p = os.path.join(base, fn)
            if os.path.exists(p):
                try:
                    os.remove(p)
                except OSError:
                    pass


def harvest(run_dir, dest):
    """Copy all known TRIM output files from run_dir (root or SRIM Outputs) to dest."""
    os.makedirs(dest, exist_ok=True)
    copied = []
    for fn in OUTPUT_FILES:
        for base in (run_dir, os.path.join(run_dir, "SRIM Outputs")):
            src = os.path.join(base, fn)
            if os.path.exists(src):
                shutil.copy2(src, os.path.join(dest, fn))
                copied.append(fn)
                break
    return copied


def run_job(run_dir, job, write_only=False):
    """Configure and (optionally) execute one job in run_dir; harvest results."""
    name = job["name"]
    dest = os.path.join(RESULTS_DIR, name)

    # 1. Write TRIM.IN
    trimin = build_trimin(job)
    with open(os.path.join(run_dir, "TRIM.IN"), "w") as f:
        f.write(trimin)
    if write_only:
        log(f"{name}: TRIM.IN written into {os.path.basename(run_dir)} (write-only)")
        return

    # 2. Batch mode + clean stale outputs
    set_trimauto(run_dir, 1)
    clean_outputs(run_dir)

    # 3. Execute TRIM.exe (blocks until the batch run self-terminates)
    log(f"{name}: START in {os.path.basename(run_dir)}  "
        f"({job['ion']} {job['energy']} keV, {job['ions']} ions)")
    t0 = time.time()
    try:
        proc = launch_with_watchdog(os.path.join(run_dir, TRIM_EXE), run_dir, name)
        rc = proc.wait()  # block until the batch run self-terminates
    except Exception as exc:  # noqa: BLE001
        log(f"{name}: ERROR launching TRIM.exe: {exc}")
        return
    dt = time.time() - t0

    # 4. Harvest outputs
    copied = harvest(run_dir, dest)
    log(f"{name}: DONE rc={rc} in {dt:5.0f}s -> results/{name}/ "
        f"({len(copied)} files: {'SPUTTER.txt ' if 'SPUTTER.txt' in copied else ''}"
        f"{'RANGE.txt ' if 'RANGE.txt' in copied else ''}...)")


def worker(run_dir, job_q, write_only, start_delay=0.0):
    # Stagger first launches so multiple TRIM.exe don't init in the same instant.
    if start_delay:
        time.sleep(start_delay)
    while True:
        try:
            job = job_q.get_nowait()
        except queue.Empty:
            return
        try:
            run_job(run_dir, job, write_only)
        finally:
            job_q.task_done()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true",
                    help="only write TRIM.IN files, do not execute")
    ap.add_argument("--jobs", type=int, default=0,
                    help="run only the first N jobs of the queue (0 = all)")
    args = ap.parse_args()

    os.makedirs(RESULTS_DIR, exist_ok=True)
    jobs = make_jobs()
    if args.jobs > 0:
        jobs = jobs[:args.jobs]

    log(f"{len(jobs)} jobs, {len(RUN_DIRS)} parallel workers "
        f"(run_1..run_{len(RUN_DIRS)}). Results -> {RESULTS_DIR}")
    for j in jobs:
        thick = j["layers"][0][1]
        log(f"  queued {j['name']:20s} target={thick:>10,} Ang  ions={j['ions']:>6,}")

    job_q = queue.Queue()
    for j in jobs:
        job_q.put(j)

    t0 = time.time()
    threads = [threading.Thread(target=worker,
                                args=(rd, job_q, args.write, i * LAUNCH_STAGGER),
                                daemon=True)
               for i, rd in enumerate(RUN_DIRS)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    log(f"ALL DONE in {time.time() - t0:.0f}s. Outputs organized under {RESULTS_DIR}\\<job>\\")


if __name__ == "__main__":
    main()
