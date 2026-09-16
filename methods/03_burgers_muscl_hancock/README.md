# 03 - Burgers equation: MUSCL-Hancock

Second-order finite-volume reconstruction for Burgers' equation using:
- exact Godunov interface fluxes;
- Hancock predictor;
- Minmod, Van Leer, Van Albada, and Superbee TVD limiters;
- shock, rarefaction, and smooth periodic validation problems;
- mass-balance and total-variation diagnostics;
- mesh-convergence analysis.

Recommended entry point: `main.m`.

The main smooth-solution convergence study is in `experiment_smooth_convergence.m`.
