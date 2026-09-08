#!/usr/bin/env python
r"""
run_p6.py -- Project 6 energy scan in BATCH mode (TRIMAUTO=1).

Why batch: interactive TRIM in this build does NOT flush the depth-resolved summary
tables (TDATA/IONIZ/VACANCY/E2RECOIL) -- confirmed over three attempts. Batch mode
self-terminates through TRIM's own end-of-run path and has written complete tables
on 11/11 previous runs. Project 6's answers ARE those tables, so batch is required.

The TRIM.IN built interactively in run_1 is used verbatim as the template (correct
4-layer skin/bone/brain/tumour target, densities, full cascades, Ion-Ranges diskfile);
only the ion energy and the description string are substituted per job.

    python run_p6.py --write        # generate TRIM.IN files only, run nothing
    python run_p6.py                # run the full scan across run_1..run_6
    python run_p6.py --ions 5000    # override ion count (default: template's 10000)
"""
import argparse, os, queue, re, threading, time
import run_srim as rs

TEMPLATE = os.path.join(rs.SRIM_ROOT, "results", "P6_round1_partial", "TRIM_30MeV.IN")
ENERGIES = list(range(30, 45))          # 30..44 MeV, 15 jobs
DEST     = os.path.join(rs.SRIM_ROOT, "results")


def build(energy_mev, ions=None):
    """Return TRIM.IN text for one energy, from the known-good interactive template."""
    # newline="" keeps CRLF intact -- text mode would normalise it and break the split
    txt = open(TEMPLATE, encoding="latin-1", newline="").read()
    lines = txt.split("\r\n")
    f = lines[2].split()                                   # Z1 M1 E angle N bragg autosave
    f[2] = str(int(energy_mev * 1000))                     # keV
    if ions:
        f[4] = str(ions); f[6] = str(ions)
    lines[2] = "     " + "   ".join(f)
    # description string (keep the trailing element/layer counts intact)
    m = re.match(r'"(.*)"(\s+.*)$', lines[8])
    if m:
        lines[8] = f'"{energy_mev} MeV H into skin/skull/brain/tumour"{m.group(2)}'
    # PlotType 5 = no on-screen graphics: much faster for a 15-run batch, and the
    # .txt tables are written regardless (all 11 earlier batch runs used PlotType 5).
    lines[10] = re.sub(r"^(\s*)\d+", lambda mm: mm.group(1) + "5", lines[10])
    return "\r\n".join(lines)


def run_job(run_dir, job, write_only=False):
    name, E = job["name"], job["E"]
    with open(os.path.join(run_dir, "TRIM.IN"), "w", newline="") as fh:
        fh.write(job["trimin"])
    if write_only:
        rs.log(f"{name}: TRIM.IN written into {os.path.basename(run_dir)}")
        return
    rs.set_trimauto(run_dir, 1)                # batch: self-terminates AND flushes tables
    rs.clean_outputs(run_dir)
    rs.log(f"{name}: START in {os.path.basename(run_dir)}")
    t0 = time.time()
    proc = rs.launch_with_watchdog(os.path.join(run_dir, rs.TRIM_EXE), run_dir, name)
    rc = proc.wait()
    dest = os.path.join(DEST, name)
    copied = rs.harvest(run_dir, dest)
    ok = os.path.isfile(os.path.join(dest, "TDATA.txt"))
    rs.log(f"{name}: DONE rc={rc} in {time.time()-t0:5.0f}s -> results/{name}/ "
           f"({len(copied)} files) {'OK' if ok else '** NO TDATA -- FAILED **'}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--ions", type=int, default=None)
    ap.add_argument("--jobs", type=int, default=0)
    ap.add_argument("--energies", type=float, nargs="+", default=None,
                    help="explicit list of energies in MeV (overrides the default 30-44 scan)")
    a = ap.parse_args()

    energies = a.energies if a.energies else ENERGIES
    jobs = [dict(name=f"P6_{E:g}MeV", E=E, trimin=build(E, a.ions)) for E in energies]
    if a.jobs:
        jobs = jobs[:a.jobs]
    rs.log(f"{len(jobs)} Project-6 jobs ({energies[0]}-{energies[-1]} MeV) "
           f"over {len(rs.RUN_DIRS)} workers -> {DEST}")

    q = queue.Queue()
    for j in jobs:
        q.put(j)

    def worker(rd, delay):
        time.sleep(delay)
        while True:
            try: job = q.get_nowait()
            except queue.Empty: return
            try: run_job(rd, job, a.write)
            finally: q.task_done()

    ts = [threading.Thread(target=worker, args=(rd, i * rs.LAUNCH_STAGGER), daemon=True)
          for i, rd in enumerate(rs.RUN_DIRS)]
    for t in ts: t.start()
    for t in ts: t.join()
    rs.log("P6 SCAN COMPLETE")


if __name__ == "__main__":
    main()
