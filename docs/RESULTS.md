# Numerical results

This document collects the main quantitative results reported in the technical reports and generated during the internship.

## Burgers equation - MUSCL-Hancock

### Shock, N = 200, CFL = 0.8, t = 0.6

| Method | L1 error | Linf error | Final TV | Conservative residual |
|---|---:|---:|---:|---:|
| Godunov | 3.52435e-3 | 1.75630e-1 | 1.000 | 4.44e-16 |
| Minmod | 2.11880e-3 | 1.05930e-1 | 1.000 | 4.44e-16 |
| Van Leer | 1.85020e-3 | 9.25090e-2 | 1.000 | 2.22e-16 |
| Van Albada | 1.95970e-3 | 9.79810e-2 | 1.000 | 4.44e-16 |
| Superbee | **1.58900e-3** | **7.94500e-2** | 1.000 | 4.44e-16 |

All tested limited schemes remain TVD in this experiment. Superbee is the sharpest of the tested reconstructions for the discontinuity, while Minmod is the most diffusive.

### Rarefaction, N = 200, CFL = 0.8, t = 0.6

| Method | L1 error | Linf error |
|---|---:|---:|
| Godunov | 1.16278e-2 | 4.58831e-2 |
| Minmod | 2.42420e-3 | 2.01010e-2 |
| Van Leer | 1.76150e-3 | 1.60010e-2 |
| Van Albada | 1.93280e-3 | 1.71150e-2 |
| Superbee | **8.96470e-4** | **8.49830e-3** |

### Smooth solution, N = 200, CFL = 0.6, t = 0.5

| Method | L1 error | Linf error | Final TV |
|---|---:|---:|---:|
| Godunov | 2.50051e-3 | 2.69449e-3 | ~0.9938 |
| Minmod | 8.32440e-5 | 5.66190e-4 | 0.99803 |
| Van Leer | **4.17380e-5** | **3.15540e-4** | 0.99889 |
| Van Albada | 5.73180e-5 | 4.31100e-4 | 0.99848 |
| Superbee | 7.61420e-5 | 4.95210e-4 | **0.99948** |

The smooth-solution refinement study shows first-order behaviour for Godunov and near-second-order L1 convergence for the MUSCL-Hancock variants.

## Euler 1D - Sod shock tube

Parameters: gamma = 1.4, CFL = 0.45, t = 0.2.

The MUSCL-Hancock-HLL method reproduces the rarefaction, contact discontinuity, and shock without negative density or pressure in the reported tests.

### Self-convergence after conservative restriction

The public repository version restricts the fine-grid **conservative variables** `(rho, rho*u, E)` to the coarse grid before recovering velocity and pressure. This is more consistent than directly averaging primitive velocity and pressure.

| Grid comparison | Density L1 | Velocity L1 | Pressure L1 | Energy L1 |
|---|---:|---:|---:|---:|
| 100 -> 200 | 2.386821e-3 | 4.464160e-3 | 1.833210e-3 | 4.921107e-3 |
| 200 -> 400 | 1.246440e-3 | 2.185093e-3 | 8.997976e-4 | 2.453840e-3 |
| 400 -> 800 | 6.541981e-4 | 1.044192e-3 | 4.375542e-4 | 1.214176e-3 |

Observed orders:

| Triplet | Density | Velocity | Pressure | Energy |
|---|---:|---:|---:|---:|
| 100-200-400 | 0.9373 | 1.0307 | 1.0267 | 1.0039 |
| 200-400-800 | 0.9300 | 1.0653 | 1.0401 | 1.0151 |

The near-first-order global rate is expected for a solution containing shocks and a contact discontinuity: the limiter intentionally reduces the reconstruction order near non-smooth regions.

## MUSCL-Hancock vs ADER4 - smooth linear advection

### 50-period propagation, N = 80

| Method | L1 | L2 | Linf | Relative amplitude |
|---|---:|---:|---:|---:|
| MUSCL-Hancock | 1.65932e-1 | 2.11394e-1 | 3.60285e-1 | 0.639344 |
| ADER4 | **1.78131e-4** | **1.97909e-4** | **2.79881e-4** | **1.000002** |

For this smooth long-time test, ADER4 has approximately 932 times lower L1 error.

### Experimental L1 convergence

| N | MUSCL L1 | MUSCL order | ADER4 L1 | ADER4 order |
|---:|---:|---:|---:|---:|
| 20 | 6.3069e-2 | - | 8.9445e-4 | - |
| 40 | 2.4450e-2 | 1.3671 | 5.6783e-5 | 3.9775 |
| 80 | 6.9791e-3 | 1.8087 | 3.5626e-6 | 3.9944 |
| 160 | 1.8980e-3 | 1.8786 | 2.2288e-7 | 3.9986 |
| 320 | 5.1622e-4 | 1.8784 | 1.3933e-8 | 3.9996 |
| 640 | 1.3680e-4 | 1.9159 | 8.7079e-10 | 4.0001 |

The ADER implementation therefore reproduces the expected fourth-order convergence with high numerical consistency on the chosen smooth periodic problem.
