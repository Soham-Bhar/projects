#!/usr/bin/env python
r"""run_p6_extra.py -- three follow-up Project 6 jobs.

  35.5, 36.5 MeV : ORIGINAL 4-layer target (10 mm) -- identical geometry to the
                   33-38 MeV series so they are directly comparable.
  39 MeV         : 5-layer target with a 4 mm brain backing (14 mm total), so the
                   beam stops inside instead of being truncated. Confirms the ions
                   pass straight through the tumour at this energy.
"""
import os, queue, sys, threading, time
import run_srim as rs, run_p6, build5

JOBS = [
    ("P6_35.5MeV",     lambda: run_p6.build(35.5, 10000)),
    ("P6_36.5MeV",     lambda: run_p6.build(36.5, 10000)),
    ("P6_39MeV_ext",   lambda: build5.build5(39, backing_mm=4.0, ions=10000)),
]

def run_job(run_dir, name, maker, write_only=False):
    with open(os.path.join(run_dir, "TRIM.IN"), "w", newline="") as fh:
        fh.write(maker())
    if write_only:
        rs.log(f"{name}: TRIM.IN written into {os.path.basename(run_dir)}"); return
    rs.set_trimauto(run_dir, 1)
    rs.clean_outputs(run_dir)
    rs.log(f"{name}: START in {os.path.basename(run_dir)}")
    t0 = time.time()
    rc = rs.launch_with_watchdog(os.path.join(run_dir, rs.TRIM_EXE), run_dir, name).wait()
    dest = os.path.join(rs.SRIM_ROOT, "results", name)
    copied = rs.harvest(run_dir, dest)
    ok = os.path.isfile(os.path.join(dest, "TDATA.txt"))
    rs.log(f"{name}: DONE rc={rc} in {time.time()-t0:5.0f}s -> results/{name}/ "
           f"({len(copied)} files) {'OK' if ok else '** NO TDATA -- FAILED **'}")

def main():
    write_only = "--write" in sys.argv
    q = queue.Queue()
    for j in JOBS: q.put(j)
    rs.log(f"{len(JOBS)} follow-up jobs over 3 workers")
    def worker(rd, delay):
        time.sleep(delay)
        while True:
            try: name, maker = q.get_nowait()
            except queue.Empty: return
            try: run_job(rd, name, maker, write_only)
            finally: q.task_done()
    ts = [threading.Thread(target=worker, args=(rs.RUN_DIRS[i], i*rs.LAUNCH_STAGGER), daemon=True)
          for i in range(3)]
    for t in ts: t.start()
    for t in ts: t.join()
    rs.log("FOLLOW-UP JOBS COMPLETE")

if __name__ == "__main__":
    main()
