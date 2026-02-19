# Computational Cardiovascular Modeling

This repository contains exploratory modeling work developed as part of a transition
from applied cardiovascular physiology and signal processing toward mechanistic and
computational cardiovascular modeling.

The goal of this repository is not to provide production-ready tools, but to document
learning, prototyping, and conceptual development related to:
- Cardiac electrophysiology (reduced-order and nonlinear models)
- Cardiovascular hemodynamics and lumped-parameter vascular models

All examples are implemented in Python and organized for reproducibility using a
Conda environment (see `environment.yml`). 

To understand the models conceptually, I began by manually coding the integration loops using Euler's method.
These scripts can be found in the Euler folder with separate folders for Python and Matlab files.

I have since moved to calling prior defined ODE solvers (e.g. solve_ips from scipy.integrate or ode15s in matlab). 
These scripts can be found in the ODE Solver folder.


