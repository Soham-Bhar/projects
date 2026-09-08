#!/usr/bin/env python
r"""
p6_manage.py -- helper for running Project 6 MANUALLY in rounds of 6.

Run from a terminal (PowerShell, cmd or Git Bash) inside D:\SRIM_sim:

    python p6_manage.py prep                        # BEFORE every round
    python p6_manage.py layers                      # verify layer widths/units BEFORE running
    python p6_manage.py check                       # did the tables actually get written?
    python p6_manage.py harvest 30 31 32 33 34 35   # AFTER every round
                                                    # (energies in run_1..run_N order)

Why prep matters: TRIM overwrites hardcoded filenames, and an interactive run can
silently fail to flush its summary tables. With stale files present, a failed run
looks identical to a good one. Clearing them first means each file is afterwards
either fresh or absent -- a silent failure becomes an obvious one.
"""
import os, re, sys, shutil

ROOT = os.path.dirname(os.path.abspath(__file__))
RUNS = [os.path.join(ROOT, f"run_{i}") for i in range(1, 7)]
OUTS = ["TDATA","RANGE","RANGE_3D","VACANCY","IONIZ","PHONON","E2RECOIL","NOVAC",
        "LATERAL","COLLISON","BACKSCAT","TRANSMIT","SPUTTER","EXYZ"]
# SRIM writes layer widths in Angstrom but may use scientific notation (1E+07),
# so the width group must accept both plain integers and E-notation.
LAYER_RE = re.compile(r'^\s*(\d+)\s+"([^"]*)"\s+([\d.]+(?:[Ee][+-]?\d+)?)\s+([\d.]+)')

def _bases(d):
    return (d, os.path.join(d, "SRIM Outputs"))

def _title(p):
    try:
        with open(p, encoding="latin-1") as f:
            return f.readline().replace("=", "").strip()
    except OSError:
        return None

def prep():
    print("PREP: clearing scratch + setting interactive mode in run_1..run_6")
    for i, d in enumerate(RUNS, 1):
        if not os.path.isdir(d):
            continue
        with open(os.path.join(d, "TRIMAUTO"), "w") as f:
            f.write("0\n\n")
        n = 0
        for name in OUTS:
            for b in _bases(d):
                p = os.path.join(b, name + ".txt")
                if os.path.isfile(p):
                    try: os.remove(p); n += 1
                    except OSError: pass
        print(f"  run_{i} : TRIMAUTO=0, {n:2d} stale output file(s) removed")
    print("Ready. Set up each run in its own folder, then launch TRIM.")

def layers():
    """Show layer widths as SRIM actually STORED them. TRIM.IN is always in Angstrom,
    so this proves whether the mm/um unit dropdown was applied on save."""
    print("LAYERS: widths as written to TRIM.IN (always Angstrom -- unit dropdown already applied)")
    for i, d in enumerate(RUNS, 1):
        f = os.path.join(d, "TRIM.IN")
        if not os.path.isfile(f):
            print(f"  run_{i} : no TRIM.IN"); continue
        lines = open(f, encoding="latin-1").read().splitlines()
        desc = lines[8].split('"')[1].strip() if len(lines) > 8 and '"' in lines[8] else "?"
        print(f"  run_{i} : {desc}")
        tot = 0
        for ln in lines:
            m = LAYER_RE.match(ln)
            if not m:
                continue
            n, name, ang, rho = m.group(1), m.group(2), float(m.group(3)), m.group(4)
            tot += ang
            flag = "   <-- SUSPICIOUS (sub-micron)" if ang < 10000 else ""
            print(f"        layer {n}: {name:26s} {ang:>15,.0f} A = {ang/1e7:8.4f} mm  rho={rho}{flag}")
        if tot:
            print(f"        TOTAL {tot:,.0f} A = {tot/1e7:.4f} mm")

def check():
    print("CHECK: which run does each folder currently hold?")
    for i, d in enumerate(RUNS, 1):
        t = _title(os.path.join(d, "TDATA.txt"))
        print(f"  run_{i} : " + (t if t else "** no TDATA.txt -- tables NOT written yet **"))

def harvest(energies):
    ok = bad = 0
    for i, E in enumerate(energies, 1):
        d = RUNS[i-1]
        dest = os.path.join(ROOT, "results", f"P6_{E}MeV")
        if not os.path.isfile(os.path.join(d, "TDATA.txt")):
            print(f"  run_{i} -> P6_{E}MeV : ** FAIL: no TDATA.txt (tables never written) **")
            bad += 1; continue
        os.makedirs(dest, exist_ok=True)
        for name in OUTS:
            for b in _bases(d):
                p = os.path.join(b, name + ".txt")
                if os.path.isfile(p):
                    shutil.copy2(p, dest); break
        ti = os.path.join(d, "TRIM.IN")
        if os.path.isfile(ti): shutil.copy2(ti, dest)
        ions = "?"
        for ln in open(os.path.join(dest, "TDATA.txt"), encoding="latin-1"):
            if "Total Ions calculated" in ln:
                ions = ln.split("=")[-1].strip(); break
        print(f"  run_{i} -> P6_{E}MeV : OK  [{_title(os.path.join(dest,'TDATA.txt'))}]  ions={ions}")
        ok += 1
    print(f"harvested {ok}, failed {bad}")
    if bad:
        print("!! Re-run the failed ones -- their tables were never flushed.")
    print("CHECK the titles above: each must name the energy you intended for that folder.")

if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if   cmd == "prep":   prep()
    elif cmd == "layers": layers()
    elif cmd == "check":  check()
    elif cmd == "harvest" and len(sys.argv) > 2: harvest(sys.argv[2:])
    else:
        print("usage: python p6_manage.py prep | layers | check | harvest E1 [E2 ...]")
        sys.exit(1)
