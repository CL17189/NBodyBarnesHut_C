# NBodyBarnesHut_C
This project implements a Barnes-Hut algorithm for simulating the gravitational interactions between multiple particles (or stars) in 2D space. The implementation supports OpenMP for parallel processing and uses both Array of Structures (AoS) and Structure of Arrays (SoA) representations for performance comparison and spatial efficiency.

✨ Features
Quadtree-based spatial decomposition for $O(n \log n)$ gravitational computation.

Adaptive force approximation via the Barnes-Hut criterion ($\theta_{\text{max}}$).

Flexible configuration via command-line arguments.

Uses SoA memory layout for improved cache locality.

Accelerated using OpenMP (multi-threading).

Includes improved center-of-mass and force aggregation logic.

