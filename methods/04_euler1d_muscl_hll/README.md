# 04 - 1D Euler equations: HLL and MUSCL-Hancock-HLL

Finite-volume solution of the compressible 1D Euler equations on Sod's shock-tube problem.

Main components:
- conservative/primitive variable conversion;
- CFL-controlled adaptive time step;
- HLL approximate Riemann solver;
- first-order HLL baseline;
- MUSCL-Hancock reconstruction;
- Minmod, Van Leer, and Superbee limiter comparison;
- density/pressure positivity safeguards;
- grid-refinement and self-convergence study.

Recommended run order:

```matlab
hll_primer_orden_sod
muscl_hancock_sod
comparacion_hll_muscl
comparacion_limitadores_sod
zooms_discontinuidades_sod
estudio_refinamiento_sod
```

The refinement script in this public version performs coarse-grid restriction on the conserved variables before recovering velocity and pressure.
