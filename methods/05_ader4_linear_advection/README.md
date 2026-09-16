# 05 - Fourth-order ADER for smooth linear advection

A fourth-order ADER finite-volume scheme is compared with MUSCL-Hancock on periodic smooth advection.

ADER4 uses:
- cubic reconstruction from cell averages;
- interface spatial derivatives up to third order;
- the Cauchy-Kowalewski procedure;
- fourth-order time averaging of the interface flux.

Experiments include:
- 50-period long-time propagation;
- L1/L2/Linf error analysis;
- amplitude preservation;
- convergence from N = 20 to N = 640;
- error growth with propagation time.

Important: the cubic reconstruction is **unlimited** and this implementation is therefore intended for smooth solutions, not shocks.
