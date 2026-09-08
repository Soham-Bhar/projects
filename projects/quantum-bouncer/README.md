# Quantum Bouncer

Time evolution of a quantum particle bouncing in a uniform gravitational field. Submitted for
**PC2130 Quantum Mechanics 1**.

## The problem

The "quantum bouncer" is a particle in a linear gravitational potential above a hard floor. It
is one of the few systems with an exact analytic solution in terms of Airy functions, which
makes it a clean setting to watch genuinely non-classical behaviour: a wave packet that spreads,
collapses, and then periodically reassembles itself.

## Method

- Energy eigenvalues obtained from the zeros of the Airy function, with the eigenstates
  normalised numerically.
- An initial Gaussian wave packet decomposed onto that eigenbasis.
- Time evolution by applying the phase factor of each eigenstate and resumming, so the
  propagation is exact rather than a finite-difference approximation.
- NumPy for the array operations, SciPy for numerical integration and the Airy special
  functions, Matplotlib for the visualisation.

## What it shows

The packet initially behaves like a classical bouncing ball, then decoheres into an
interference pattern with no classical analogue, and later undergoes **quantum revivals**,
reassembling close to its original shape after a characteristic revival time.

## Contents

```
code/quantum_bouncer.py   simulation and plotting
report/                   compiled report (PDF)
```

## Running it

```bash
pip install numpy scipy matplotlib
python code/quantum_bouncer.py
```

---
Code: MIT. Report and figures: CC BY-NC-ND 4.0. See [DISCLAIMER.md](../../DISCLAIMER.md).
