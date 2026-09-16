# Finite-Volume Methods for Hyperbolic PDEs in MATLAB

**A progressive numerical study of conservation laws — from first-order transport schemes to MUSCL-Hancock, HLL Riemann solvers, fourth-order ADER, and an introductory CFD validation workflow in ANSYS Fluent.**

This repository contains the numerical work developed during a research internship at **CITMAGA (Centro de Investigación e Tecnoloxía Matemática de Galicia)**, under the supervision of **Saray Busto Ulloa**, in July–August 2026.

The project was designed as a structured progression in **scientific computing, numerical analysis, and computational fluid dynamics**. It starts with the 1D linear advection equation, moves to nonlinear scalar conservation laws through Burgers' equation, extends the finite-volume framework to the **1D compressible Euler equations**, and progresses toward **fourth-order ADER methods** for smooth periodic advection.

As a final applied extension, the internship also introduced a professional CFD workflow using **ANSYS Workbench, DesignModeler, Meshing, and Fluent**, including numerical validation of two-dimensional Poiseuille flow against its analytical solution.

The emphasis throughout the project is not only on obtaining numerical solutions, but on understanding and measuring:

- conservation and physically correct shock propagation;
- CFL stability and numerical diffusion;
- exact and approximate Riemann solvers;
- TVD slope limiting and high-resolution reconstruction;
- shock, contact, and rarefaction resolution;
- convergence rates and error norms;
- long-time dissipation and amplitude preservation;
- the trade-off between robustness, accuracy, and formal order;
- the connection between explicit numerical implementations and professional CFD software.

---

## Full internship report

The complete internship work is documented in a **69-page technical report**, covering the mathematical foundations, MATLAB implementations, validation experiments, convergence studies, interpretation of results, and the introductory CFD study with ANSYS Fluent.

[**Read the full CITMAGA Internship Report (Spanish)**](reports/CITMAGA_Internship_Report_ES.pdf)

---

## Project at a glance

| Stage | Governing problem | Numerical methods / tools | Main numerical question |
|---|---|---|---|
| **01** | Linear advection | Centered, Lax-Friedrichs, Upwind | How do stability and numerical diffusion emerge? |
| **02** | Burgers' equation | Godunov, Lax-Friedrichs, Q-scheme | Why is conservation essential for the correct shock speed? |
| **03** | Burgers' equation | MUSCL-Hancock + TVD limiters | How much accuracy can be gained without introducing spurious oscillations? |
| **04** | 1D Euler equations | HLL, MUSCL-Hancock-HLL | How are shocks, contacts, and rarefactions captured in a nonlinear hyperbolic system? |
| **05** | Smooth linear advection | MUSCL-Hancock vs ADER4 | What is gained by moving from second- to fourth-order accuracy? |
| **06** | Reconstruction experiments | Limited / centered reconstruction variants | How does slope reconstruction affect numerical behaviour? |
| **07** | 2D laminar flow | ANSYS Fluent, first- and second-order discretisation | How do the numerical concepts studied from first principles translate to a professional CFD environment? |

### Numerical progression

```text
Linear advection
      │
      ├── stability, CFL condition, numerical diffusion
      ▼
Burgers' equation
      │
      ├── entropy solutions, shocks, rarefactions, conservation
      ▼
MUSCL-Hancock
      │
      ├── second-order reconstruction + TVD limiters
      ▼
1D Euler equations
      │
      ├── HLL Riemann solver + Sod shock tube
      ▼
High-order finite-volume methods
      │
      ├── fourth-order ADER + convergence analysis
      ▼
Professional CFD workflow
      │
      └── ANSYS Fluent + analytical validation of Poiseuille flow
```

---

## Core mathematical setting

The main problems in this repository can be written as hyperbolic conservation laws of the form

```math
\frac{\partial \mathbf{U}}{\partial t}
+
\frac{\partial \mathbf{F}(\mathbf{U})}{\partial x}
=
0.
```

Here, **U** denotes the vector of conserved variables, while **F(U)** represents the corresponding physical flux.

The finite-volume update used throughout the project has the generic form

```math
\mathbf{U}_i^{n+1}
=
\mathbf{U}_i^n
-
\frac{\Delta t}{\Delta x}
\left(
\mathbf{F}_{i+1/2}
-
\mathbf{F}_{i-1/2}
\right).
```

The numerical work therefore revolves around one central question: **how should the interface fluxes be constructed so that the method remains conservative, stable, accurate, and physically meaningful?**

---

# Featured numerical results

## 1. Burgers' equation: first-order vs high-resolution reconstruction

For the smooth Burgers test at `N = 200`, Van Leer MUSCL-Hancock reduces the L1 error from **2.50e-3** for first-order Godunov to **4.17e-5**.

| Method | L1 error | Linf error | Final TV |
|---|---:|---:|---:|
| Godunov | 2.50051e-3 | 2.69449e-3 | ~0.9938 |
| Minmod | 8.32440e-5 | 5.66190e-4 | 0.99803 |
| **Van Leer** | **4.17380e-5** | **3.15540e-4** | 0.99889 |
| Van Albada | 5.73180e-5 | 4.31100e-4 | 0.99848 |
| Superbee | 7.61420e-5 | 4.95210e-4 | **0.99948** |

<p align="center">
  <img src="results/figures/03_burgers_muscl_hancock/shock_comparison.png" width="760" alt="Burgers shock comparison">
</p>

The shock and rarefaction experiments reproduce the expected limiter trade-off: **Minmod is the most diffusive**, while **Superbee produces the sharpest discontinuities** among the tested TVD reconstructions. For the smooth test, **Van Leer provides the lowest L1 error**.

The refinement study recovers approximately **first-order convergence for Godunov** and **second-order convergence for MUSCL-Hancock** while the solution remains smooth.

---

## 2. 1D Euler equations: HLL vs MUSCL-Hancock-HLL

The compressible Euler module solves **Sod's shock-tube problem** with `gamma = 1.4`, adaptive CFL time stepping, conservative finite-volume updates, reconstruction in primitive variables, and an HLL approximate Riemann solver.

<p align="center">
  <img src="results/figures/04_euler1d_muscl_hll/hll_vs_muscl.png" width="760" alt="HLL versus MUSCL-Hancock-HLL on the Sod shock tube">
</p>

The solution reproduces the three characteristic structures of Sod's problem:

- a rarefaction wave travelling to the left;
- a contact discontinuity travelling to the right;
- a shock travelling further to the right.

MUSCL-Hancock-HLL reduces numerical diffusion relative to first-order HLL, particularly around the contact discontinuity and the shock, while remaining non-oscillatory in the reported tests.

The implementation also includes:

- Minmod, Van Leer, and Superbee limiter comparisons;
- positivity safeguards for density and pressure;
- detailed contact and shock zooms;
- grid-refinement / self-convergence analysis.

<p align="center">
  <img src="results/figures/04_euler1d_muscl_hll/contact_shock_zooms.png" width="760" alt="Contact and shock zooms for the Sod problem">
</p>

The observed global convergence rate is close to **one**, despite the formally second-order MUSCL-Hancock reconstruction. This is expected for a solution containing shocks and a contact discontinuity: the limiter deliberately reduces the local order near non-smooth regions, and those regions dominate global error norms.

---

## 3. MUSCL-Hancock vs fourth-order ADER

The high-order module studies smooth periodic linear advection and compares a limited second-order MUSCL-Hancock solver with a **fourth-order ADER** method based on cubic reconstruction and the Cauchy-Kowalewski procedure.

### Long-time propagation: 50 periods, `N = 80`

| Method | L1 | L2 | Linf | Relative amplitude |
|---|---:|---:|---:|---:|
| MUSCL-Hancock | 1.65932e-1 | 2.11394e-1 | 3.60285e-1 | 0.639344 |
| **ADER4** | **1.78131e-4** | **1.97909e-4** | **2.79881e-4** | **1.000002** |

For this smooth long-time benchmark, ADER4 gives approximately **932× lower L1 error** than MUSCL-Hancock and preserves the wave amplitude almost exactly after fifty full periods.

<p align="center">
  <img src="results/figures/05_ader4/long_time_50_periods.png" width="760" alt="Long-time MUSCL-Hancock versus ADER4 comparison">
</p>

### Experimental convergence

| N | MUSCL L1 | MUSCL order | ADER4 L1 | ADER4 order |
|---:|---:|---:|---:|---:|
| 20 | 6.3069e-2 | — | 8.9445e-4 | — |
| 40 | 2.4450e-2 | 1.3671 | 5.6783e-5 | 3.9775 |
| 80 | 6.9791e-3 | 1.8087 | 3.5626e-6 | 3.9944 |
| 160 | 1.8980e-3 | 1.8786 | 2.2288e-7 | 3.9986 |
| 320 | 5.1622e-4 | 1.8784 | 1.3933e-8 | 3.9996 |
| 640 | 1.3680e-4 | 1.9159 | 8.7079e-10 | **4.0001** |

<p align="center">
  <img src="results/figures/05_ader4/convergence.png" width="760" alt="MUSCL-Hancock and ADER4 convergence study">
</p>

The ADER implementation reproduces the expected **fourth-order convergence** with high numerical consistency on the selected smooth periodic problem.

> **Scope note:** the ADER4 implementation uses unlimited cubic reconstruction and is deliberately evaluated on smooth solutions. It is **not** presented as a shock-capturing replacement for limited MUSCL schemes.

---

## 4. From explicit numerical methods to CFD: ANSYS Fluent

As a final extension of the internship, the numerical concepts studied through explicit MATLAB implementations were explored within a professional **Computational Fluid Dynamics (CFD)** environment using **ANSYS Workbench**.

The workflow included:

- geometry definition with **DesignModeler**;
- mesh generation and refinement concepts in **ANSYS Meshing**;
- boundary definition through Named Selections;
- physical-model and boundary-condition setup in **ANSYS Fluent**;
- iterative convergence monitoring through residuals;
- first- and second-order spatial discretisation;
- post-processing and direct comparison with an analytical solution.

### Validation case: 2D Poiseuille flow

A steady, incompressible, laminar flow between two parallel walls was used as a representative validation problem. The analytical solution produces the characteristic parabolic velocity profile, with zero velocity at the walls and a maximum at the centre of the channel.

The Fluent solution was evaluated on a transverse section at `x = 5` and compared directly with the analytical profile using a Custom Field Function. The numerical result reproduced the expected parabolic behaviour with close agreement to the analytical reference.

This final stage connected the methods implemented explicitly in MATLAB — discretisation, order of accuracy, boundary conditions, convergence, and quantitative validation — with the workflow of an integrated engineering CFD solver.

> The ANSYS component is documented in the full internship report rather than distributed as a standalone Fluent project in this repository.

---

# Methods implemented

## Linear advection

Introductory experiments establish the basic numerical concepts used throughout the project:

- explicit centered discretisation as an instability demonstration;
- Lax-Friedrichs;
- Upwind;
- finite-volume interpretation of upwind fluxes;
- CFL reasoning;
- numerical diffusion;
- comparison against the analytical transported profile.

## Burgers' equation — first order

The nonlinear scalar problem introduces:

- exact entropy solutions for Riemann problems;
- shock and rarefaction waves;
- Rankine-Hugoniot shock speed;
- Godunov flux;
- Lax-Friedrichs flux;
- Q-scheme;
- conservative vs non-conservative discretisations.

A central result of this module is that **a numerically stable method can still propagate a shock at the wrong speed if the discretisation is not conservative**.

## MUSCL-Hancock

The second-order extension combines:

1. piecewise-linear reconstruction;
2. TVD slope limiting;
3. a half-step Hancock predictor;
4. interface Riemann problems;
5. a conservative finite-volume update.

Implemented limiters:

- Minmod;
- Van Leer;
- Van Albada;
- Superbee.

The experiments quantify the trade-off between **diffusion, sharpness, total variation, and smooth-region accuracy**.

## HLL approximate Riemann solver

For the Euler equations, the HLL flux approximates the Riemann fan using two bounding signal speeds. It provides a robust baseline for compressible-flow calculations, although it does not explicitly represent the intermediate contact wave and is therefore more diffusive around the contact discontinuity.

## Fourth-order ADER

The ADER4 implementation combines:

- cubic reconstruction from cell averages;
- interface derivatives up to third order;
- the Cauchy-Kowalewski procedure;
- fourth-order time averaging of interface fluxes;
- a one-step finite-volume update.

For the selected smooth periodic problem, the implementation reaches approximately **fourth-order convergence**.

## CFD workflow with ANSYS Fluent

The final stage complements the MATLAB implementations with an introduction to professional CFD modelling:

- geometry generation;
- mesh construction;
- boundary-condition definition;
- solver configuration;
- first- and second-order discretisation choices;
- residual monitoring;
- post-processing;
- analytical validation with Poiseuille flow.

The purpose of this stage is not to replace the explicit numerical implementations, but to connect the underlying numerical concepts with their use in an integrated engineering simulation environment.

---

# Validation strategy

The implementations are not assessed from plots alone. Depending on the problem, validation includes:

- exact or analytical reference solutions;
- L1, L2, and Linf error norms;
- experimental orders of convergence;
- conservation / mass-balance residuals;
- Rankine-Hugoniot shock-speed verification;
- total-variation monitoring;
- positivity checks for Euler density and pressure;
- limiter comparisons;
- grid-refinement and self-convergence studies;
- long-time propagation tests;
- amplitude preservation and cumulative numerical diffusion;
- analytical validation of the CFD Poiseuille velocity profile;
- comparison of first- and second-order discretisation settings in ANSYS Fluent.

Detailed quantitative results are collected in [`docs/RESULTS.md`](docs/RESULTS.md), while the complete mathematical development and discussion are available in the [full internship report](reports/CITMAGA_Internship_Report_ES.pdf).

---

# Repository structure

```text
.
├── README.md
├── CITATION.cff
├── .gitignore
│
├── methods/
│   ├── 01_linear_advection/
│   ├── 02_burgers_first_order/
│   ├── 03_burgers_muscl_hancock/
│   ├── 04_euler1d_muscl_hll/
│   ├── 05_ader4_linear_advection/
│   └── 06_reconstruction_experiments/
│
├── results/
│   └── figures/
│       ├── 01_linear_advection/
│       ├── 02_burgers_first_order/
│       ├── 03_burgers_muscl_hancock/
│       ├── 04_euler1d_muscl_hll/
│       └── 05_ader4/
│
├── reports/
│   ├── CITMAGA_Internship_Report_ES.pdf
│   ├── 02_burgers_first_order_analysis_es.pdf
│   ├── 03_burgers_muscl_hancock_report_es.pdf
│   ├── 04_euler1d_muscl_hll_report_es.pdf
│   └── 05_muscl_vs_ader4_report_es.pdf
│
└── docs/
    ├── METHODS.md
    ├── RESULTS.md
    ├── RUN_GUIDE.md
    └── CURATION_NOTES.md
```

Each main method directory contains the implementation associated with one stage of the numerical progression. The `reports/` directory provides the extended technical documentation behind the code and numerical results.

---

# Reproducing the MATLAB experiments

## Requirements

- MATLAB
- No external datasets are required for the main numerical experiments.

The source uses modern MATLAB features such as `arguments`, string scalars, tables, and `exportgraphics`. The exact MATLAB release used during the internship was not recorded, so no unsupported version requirement is claimed.

## Representative entry points

### Burgers — MUSCL-Hancock

```matlab
cd methods/03_burgers_muscl_hancock
main
experiment_smooth_convergence
```

### Euler 1D — Sod shock tube

```matlab
cd methods/04_euler1d_muscl_hll
hll_primer_orden_sod
muscl_hancock_sod
comparacion_hll_muscl
comparacion_limitadores_sod
estudio_refinamiento_sod
```

### MUSCL-Hancock vs ADER4

```matlab
cd methods/05_ader4_linear_advection
main_comparison
convergence_study
error_vs_periods
```

For a complete execution guide, see [`docs/RUN_GUIDE.md`](docs/RUN_GUIDE.md).

---

# Technical reports

The [`reports/`](reports/) directory contains the technical documentation produced throughout the internship.

The main document is the complete **69-page final internship report**, which integrates the mathematical background, implementations, numerical experiments, convergence studies, interpretation of results, and introductory CFD work into a single document.

| Report | Topic |
|---|---|
| **[`CITMAGA_Internship_Report_ES.pdf`](reports/CITMAGA_Internship_Report_ES.pdf)** | **Complete internship report: finite-volume methods, Burgers, MUSCL-Hancock, Euler 1D, HLL, ADER4, and ANSYS Fluent** |
| [`02_burgers_first_order_analysis_es.pdf`](reports/02_burgers_first_order_analysis_es.pdf) | First-order conservative schemes for Burgers' equation |
| [`03_burgers_muscl_hancock_report_es.pdf`](reports/03_burgers_muscl_hancock_report_es.pdf) | MUSCL-Hancock, TVD limiters, and convergence |
| [`04_euler1d_muscl_hll_report_es.pdf`](reports/04_euler1d_muscl_hll_report_es.pdf) | HLL and MUSCL-Hancock-HLL for the 1D Euler equations |
| [`05_muscl_vs_ader4_report_es.pdf`](reports/05_muscl_vs_ader4_report_es.pdf) | Long-time accuracy and convergence of MUSCL-Hancock vs ADER4 |

The reports are written in **Spanish**, while the repository documentation is provided in English for international accessibility.

---

# Scope and limitations

This repository is intended as a rigorous educational/research implementation of classical finite-volume methods, not as a production CFD package.

- The centered transport scheme is included intentionally to demonstrate instability.
- HLL is robust but does not explicitly resolve the contact wave.
- MUSCL-Hancock is formally second order in smooth regions; limiting reduces the local order near discontinuities.
- Global convergence on Sod's problem is therefore close to first order.
- The ADER4 implementation is designed for **smooth linear advection** and uses unlimited cubic reconstruction.
- Applying the current ADER4 reconstruction directly across shocks would require a suitable non-oscillatory extension, such as ENO/WENO or an appropriate limiting strategy.
- The ANSYS section is an introductory CFD application and is documented in the report; the repository does not distribute a standalone Fluent project.
- CPU timings reported in the technical reports are illustrative and depend on machine and MATLAB JIT behaviour; error norms and convergence rates are the primary validation metrics.

---

# Research context

This work was developed during a research internship at **CITMAGA (Centro de Investigación e Tecnoloxía Matemática de Galicia)**, with the objective of building a rigorous foundation in:

- numerical analysis;
- finite-volume methods;
- hyperbolic partial differential equations;
- Riemann problems and shock-capturing methods;
- scientific programming in MATLAB;
- high-order numerical schemes;
- convergence and error analysis;
- computational fluid dynamics;
- professional CFD workflows with ANSYS Fluent.

The internship connected **mathematical modelling, numerical implementation, quantitative validation, and engineering simulation**, progressing from scalar model equations to nonlinear hyperbolic systems and high-order methods.

---

# Author

**Lucía Gómez Barros**  
Biomedical Engineering — Universidade de Vigo

---

# Acknowledgements

Supervision: **Saray Busto Ulloa — CITMAGA**.

The mathematical development of the internship followed, in particular, the finite-volume and hyperbolic-conservation-law framework presented in:

- M. E. Vázquez-Cendón, *Solving Hyperbolic Equations with Finite Volume Methods*, Springer, 2015.
- E. F. Toro, *Riemann Solvers and Numerical Methods for Fluid Dynamics*, 3rd ed., Springer, 2009.

The ANSYS Fluent component was developed as an introductory professional-CFD extension of the numerical work documented in the final internship report.

---

# Citation

If you reuse material from this repository, citation metadata are provided in [`CITATION.cff`](CITATION.cff).
