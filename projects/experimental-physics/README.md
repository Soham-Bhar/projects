# Experimental Physics

Laboratory experiments and simulation work from the experimental physics sequence at NUS.
Experiments A to E are from **PC2193 Experimental Physics**; the SRIM simulation study is from
**PC3193 Experimental Physics II**, as is everything added here from now on. Each report covers
the measurement, the error analysis, and the comparison against theory.

Alongside every report is the material it was built from: raw data, the spreadsheets the analysis
was carried out in, the analysis scripts, and the figures. Several of these experiments were
analysed entirely in Excel, so the workbooks contain the plots and the fitting as well as the
numbers.

| Report | Experiment | What was measured | Supporting material |
|---|---|---|---|
| [Experiment A](reports_data_figures/electron_spin_resonance/ExpA-electron-spin-resonance.pdf) | Electron Spin Resonance | The g-factor of a free radical, via resonant microwave absorption in a uniform magnetic field | [Data and analysis workbook](reports_data_figures/electron_spin_resonance/) |
| [Experiment B](reports_data_figures/gaussian_beam_optics/ExpB-gaussian-beam-optics.pdf) | Propagation of Laser Light: Gaussian Beam Optics | Beam waist, Rayleigh range and beam quality factor, from transverse intensity profiles fitted to the Gaussian beam propagation model | [Raw data, 4 Python scripts, apparatus photographs](reports_data_figures/gaussian_beam_optics/) |
| [Experiment C](reports_data_figures/hall_effect/ExpC-hall-effect-germanium.pdf) | Hall Effect in n- and p-type Germanium | Carrier concentration, mobility and Hall coefficient, confirming the sign reversal of the Hall voltage between carrier types | [14 CSV runs across three sub-experiments](reports_data_figures/hall_effect/) |
| [Experiment D](reports_data_figures/XRay_diffraction/ExpD-x-ray-diffraction.pdf) | X-Ray Diffraction of Crystals | Lattice constants of LiF and KBr by Bragg diffraction, and Planck's constant from the Duane-Hunt relation | [4 workbooks, 8 scan and analysis figures](reports_data_figures/XRay_diffraction/) |
| [Experiment E](reports_data_figures/magnetic_moments/ExpE-magnetic-moments.pdf) | Understanding Magnetic Moments | Torque on a current-carrying loop in a Helmholtz field, against coil current, turns, angle, sample current and diameter, with the Helmholtz constant compared to theory | [Raw data and modified-approach workbooks](reports_data_figures/magnetic_moments/) |
| [Simulation](reports_data_figures/SRIM_simulations/SRIM_Lab_Report.pdf) | Interaction of Ions with Matter: Computer Simulations (PC3193) | Ion ranges, damage profiles and sputtering yields from Monte-Carlo binary collision simulations in SRIM | [25 simulation runs, 10 scripts, 51 figures](reports_data_figures/SRIM_simulations/) |

## A note on Experiment E

Aarav Bindawala was my lab partner and part of the data collection was done together. All
analysis, discussion, interpretations and conclusions are independently my own.

## Approach

Each experiment follows the same structure: characterise the apparatus, take measurements across
a controlled parameter, propagate uncertainties through to the derived quantity, and compare the
result against the accepted value with a stated error budget. The reports show the fitting
procedure and the reasoning behind each source of uncertainty, not just the final number.

The SRIM work is the exception, being simulation rather than bench measurement. It follows the
same discipline in a different setting: the input decks (`TRIM.IN`) that define every run are
included, so any result in the report can be reproduced from scratch.

## Layout

```
reports_data_figures/
  electron_spin_resonance/   report, data and analysis workbook
  gaussian_beam_optics/      report, raw data, Python scripts, apparatus photographs
  hall_effect/               report, 14 CSV runs in exp1/ exp2/ exp3/
  XRay_diffraction/          report, 4 workbooks, 8 figures
  magnetic_moments/          report, raw data and analysis workbooks
  SRIM_simulations/          report, 10 scripts, 25 runs, 51 figures
```

### A note on the SRIM data

SRIM writes a per-ion record of every collision, and for these runs those records came to several
gigabytes. What is committed here is the reproducible subset: the `TRIM.IN` input deck for
every run, the compact summary outputs SRIM produces alongside them (`TDATA`, `RANGE`, `RANGE_3D`,
`IONIZ`, `PHONON`, `VACANCY`, `E2RECOIL`, `LATERAL`, `NOVAC`), the scripts that drive the runs and
build the plots, and the figures used in the report. The bulk per-ion dumps (`SPUTTER`,
`COLLISON`, `EXYZ`) are not included, since SRIM regenerates them from the input decks.

## A note on these reports

These are my own graded submissions, published as a portfolio. They are **not** licensed for
reuse and are not a study aid for anyone currently taking PC2193 or PC3193. Read
[DISCLAIMER.md](../../DISCLAIMER.md) first.

---
Reports and figures: CC BY-NC-ND 4.0. Code and scripts: MIT.
