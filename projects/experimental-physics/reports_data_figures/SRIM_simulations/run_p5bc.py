#!/usr/bin/env python
r"""run_p5bc.py -- Projects 5(b) and 5(c): sputter-yield comparison at the optimum
argon energy (100 keV) found in Project 5(a).

  5b : H  (protons) -> GaAs @ 100 keV
  5c : Au (gold)    -> GaAs @ 100 keV

Same GaAs target, full damage cascades, 0 deg incidence, Sputtered diskfile on, so
the yields compare directly with argon's 3.487 atoms/ion.

Target thicknesses differ because the ranges differ by two orders of magnitude at the
same energy: a 100 keV proton penetrates ~0.4 um while a 100 keV gold ion stops within
a few hundred Angstrom. Both are set well beyond the expected range so no ion is
transmitted (verified after the run).

    python run_p5bc.py --write   # generate TRIM.IN only
    python run_p5bc.py           # run both across run_1 / run_2
"""
import os, queue, sys, threading, time
import run_srim as rs

IONS = 10000

JOBS = [
    dict(name="5b_H_100keV_GaAs", ion="H", energy=100, ions=IONS, autosave=IONS,
         cascade=2, plottype=5, desc="100 keV H into GaAs (sputter)",
         elements=["Ga", "As"],
         layers=[("GaAs", 20_000, 5.320, {"Ga": 0.5, "As": 0.5})],
         disk=dict(ranges=1, backscat=0, transmit=0, sputter=1, collisions=0, exyz=0)),
    dict(name="5c_Au_100keV_GaAs", ion="Au", energy=100, ions=IONS, autosave=IONS,
         cascade=2, plottype=5, desc="100 keV Au into GaAs (sputter)",
         elements=["Ga", "As"],
         layers=[("GaAs", 2_000, 5.320, {"Ga": 0.5, "As": 0.5})],
         disk=dict(ranges=1, backscat=0, transmit=0, sputter=1, collisions=0, exyz=0)),
]


def run_job(run_dir, job, write_only=False):
    name = job["name"]
    with open(os.path.join(run_dir, "TRIM.IN"), "w") as f:
        f.write(rs.build_trimin(job))
    if write_only:
        rs.log(f"{name}: TRIM.IN written into {os.path.basename(run_dir)}")
        return
    rs.set_trimauto(run_dir, 1)
    rs.clean_outputs(run_dir)
    rs.log(f"{name}: START in {os.path.basename(run_dir)} "
           f"({job['ion']} {job['energy']} keV, {job['ions']} ions)")
    t0 = time.time()
    rc = rs.launch_with_watchdog(os.path.join(run_dir, rs.TRIM_EXE), run_dir, name).wait()
    dest = os.path.join(rs.RESULTS_DIR, name)
    copied = rs.harvest(run_dir, dest)
    ok = os.path.isfile(os.path.join(dest, "TDATA.txt"))
    sp = os.path.isfile(os.path.join(dest, "SPUTTER.txt"))
    rs.log(f"{name}: DONE rc={rc} in {time.time()-t0:5.0f}s -> results/{name}/ "
           f"({len(copied)} files) {'OK' if ok else '** NO TDATA **'}"
           f"{'' if sp else '  ** NO SPUTTER.txt **'}")


def main():
    write_only = "--write" in sys.argv
    q = queue.Queue()
    for j in JOBS:
        q.put(j)
    rs.log(f"{len(JOBS)} jobs (Projects 5b, 5c) over 2 workers")

    def worker(rd, delay):
        time.sleep(delay)
        while True:
            try:
                job = q.get_nowait()
            except queue.Empty:
                return
            try:
                run_job(rd, job, write_only)
            finally:
                q.task_done()

    ts = [threading.Thread(target=worker, args=(rs.RUN_DIRS[i], i * rs.LAUNCH_STAGGER),
                           daemon=True) for i in range(2)]
    for t in ts:
        t.start()
    for t in ts:
        t.join()
    rs.log("5b / 5c COMPLETE")


if __name__ == "__main__":
    main()
