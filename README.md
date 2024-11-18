# Sheared Voronoi Model

This repository contains MATLAB code to simulate the shear-driven behavior of epithelial tissues using a Voronoi-based model. The code supports research published in:

1. Huang, Junxiang, James O. Cochran, Suzanne M. Fielding, M. Cristina Marchetti, and Dapeng Bi. "Shear-driven solidification and nonlinear elasticity in epithelial tissues." *Physical Review Letters* **128**, no. 17 (2022): 178001.

2. Nguyen, Anh Q., Junxiang Huang, and Dapeng Bi. "Origin of yield stress and mechanical plasticity in biological tissues." *arXiv preprint arXiv:2409.04383* (2024).

This model simulates cellular mechanics under simple shear deformation, capturing phenomena such as solidification, yield stress, and mechanical plasticity.

The primary script for the simulation is `run_shear_experiment.m`. This script initializes the system, applies shear strain incrementally, and saves outputs for each step. The simulation generates PNG images at each shear step, showing the state of the Voronoi tessellation. Images are saved using the naming convention:

qs_ka<ka_value>_N<N_cell>_p<p0_value>_gamma<gamma_value>.png

**Example**: `qs_ka0_N64_p3.78_gamma1.000.png`.

The code is tested with MATLAB R2021a and later versions.

## Features
- **Shear Simulation**: Incremental application of shear strain to a Voronoi tessellation representing cellular monolayers.
- **Visualization**: Generation of images for each shear step showing cell configurations and Voronoi tessellation structure.
- **Random Initialization with Reproducibility**: Configurable random seed to initialize cell positions consistently.
