# Repository curation notes

This repository is the public, MATLAB-focused version of a broader internship workspace. It contains the numerical-methods code, selected figures, and polished technical reports required to understand and reproduce the computational work.

Material outside that scope - administrative files, study notes, LaTeX build artefacts, ANSYS project/cache files, and cached MATLAB outputs - is intentionally omitted.

Two technical/documentation corrections were made during curation:

1. the Euler refinement study now restricts the conserved variables to the coarse grid before recovering primitive variables, rather than directly averaging velocity and pressure;
2. the MUSCL-Hancock Sod script now exports its own base-solution figure, preventing confusion with the first-order HLL figure in the original working directory.

The remaining numerical implementations were preserved from the internship workspace.
