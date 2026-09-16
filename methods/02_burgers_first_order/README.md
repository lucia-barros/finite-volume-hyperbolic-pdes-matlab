# 02 - Burgers equation: first-order schemes

This module moves from linear transport to the nonlinear scalar conservation law

`u_t + (u^2/2)_x = 0`.

It includes exact entropy solutions for Riemann problems and compares:
- Godunov;
- Lax-Friedrichs;
- Q-scheme;
- a non-conservative discretisation.

The central numerical lesson is that conservative discretisation is essential for obtaining the correct shock propagation speed.
