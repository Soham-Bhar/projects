# Projects

Coursework and personal projects by **Soham Bhar**, B.Sc. Physics with a second major in Data
Analytics at the National University of Singapore.

Everything here is currently physics: computational simulations and laboratory work. Other
projects will be added alongside it over time.

Each project folder holds the code, the generated figures, and the written report.

> **Reuse notice.** The reports here are my own graded coursework, published as a portfolio.
> They are **not** licensed for reuse, and copying any part of them into your own submission is
> plagiarism. See [DISCLAIMER.md](DISCLAIMER.md) before using anything in this repository.

---

## Projects

| Project | What it is | Contents |
|---|---|---|
| [Driven Quadruple Pendulum](projects/quadruple-pendulum) | Chaotic dynamics of four coupled masses, integrated with a hand-coded RK4 scheme and validated against SciPy | Notebook, Typst source, 7 figures, report |
| [Quantum Bouncer](projects/quantum-bouncer) | Time evolution of a quantum particle in a gravitational potential, including quantum revivals | Python simulation, report |
| [Experimental Physics](projects/experimental-physics) | Five experiments from PC2193 and one simulation study from PC3193: electron spin resonance, Gaussian beam optics, the Hall effect in germanium, X-ray diffraction, magnetic moments, and ion-matter interaction in SRIM | 6 reports, raw data, analysis workbooks, scripts, figures |

---

## Repository layout

```
projects/
  quadruple-pendulum/     code/ figures/ report/
  quantum-bouncer/        code/ report/
  experimental-physics/   reports_data_figures/<experiment>/
```

## Licensing

Two licences apply, because code and writing need different terms:

- **Code** (`.py`, `.ipynb`, scripts): [MIT](LICENSE). Use it freely, with attribution.
- **Reports, figures and written content**: [CC BY-NC-ND 4.0](LICENSE-CONTENT.md). You may share
  them with attribution, but not commercially, and not in modified form.

## Not included here

Work from my research internship is deliberately excluded. That data, the sample identifiers,
and the unpublished results belong to the host laboratory and are not mine to publish.
