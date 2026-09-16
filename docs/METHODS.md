# Methods overview

## Finite-volume formulation

The core modules update cell averages through numerical fluxes at cell interfaces. This makes conservation explicit and provides a common framework for linear transport, Burgers' equation, and the Euler equations.

## First-order schemes

The introductory module compares centered, Lax-Friedrichs, and upwind discretisations for linear advection. The centered explicit scheme is retained as a deliberate instability example, while Lax-Friedrichs and upwind illustrate the roles of numerical diffusion and information propagation direction.

For Burgers' equation, the project then introduces entropy solutions, shocks, rarefactions, the exact Godunov flux, Lax-Friedrichs, a Q-scheme, and the consequences of using a non-conservative discretisation across shocks.

## MUSCL-Hancock

MUSCL-Hancock raises the spatial and temporal accuracy through:

1. limited linear reconstruction in each cell;
2. a half-step Hancock predictor;
3. an interface Riemann problem;
4. a conservative finite-volume update.

For Burgers' equation, Minmod, Van Leer, Van Albada, and Superbee limiters are compared quantitatively. For the Euler system, reconstruction is combined with an HLL approximate Riemann solver.

## HLL for the 1D Euler equations

The HLL solver approximates the Riemann fan using two bounding signal speeds. It is robust and simple but does not explicitly include the contact wave, which explains some diffusion at the Sod contact discontinuity.

## ADER4

The fourth-order ADER implementation is developed for smooth linear advection. A cubic reconstruction provides interface values and spatial derivatives; the Cauchy-Kowalewski procedure replaces time derivatives with spatial derivatives, enabling a fourth-order time-averaged numerical flux in a one-step update.

This implementation is deliberately limited to smooth solutions. An unlimited cubic reconstruction should not be used directly across shocks.
