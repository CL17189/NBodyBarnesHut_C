# NBodyBarnesHut_C

Efficient 2D N‑body simulation using the Barnes–Hut algorithm in C, with both Array‑of‑Structures (AoS) and Structure‑of‑Arrays (SoA) layouts and OpenMP acceleration.

## Overview
This project implements a quadtree‑based Barnes–Hut solver to approximate gravitational forces among \(n\) particles in \(O(n\log n)\) time. You can compare performance and memory behavior between AoS and SoA data layouts, and easily tune accuracy and parallelism.

## Features
- **\(O(n\log n)\) complexity** via quadtree spatial decomposition  
- **Barnes–Hut criterion**: configurable opening angle \(\theta_{\max}\) for adaptive force approximation  
- **AoS & SoA layouts**: switch at compile‑time to evaluate cache locality vs. ease of indexing  
- **OpenMP support**: multi‑threaded force calculation  
- **Command‑line interface**: specify particle count, time steps, \(\theta_{\max}\), layout, thread count, etc.  
- **Optimized COM aggregation**: improved center‑of‑mass and force‑accumulation routines  

## Requirements
- GCC 9.0+ or Clang with C11 support  
- OpenMP (–fopenmp)  
- Make (or compatible build tool)  

# Binaries produced:
# - `nbody_aos` (AoS layout)  
# - `nbody_soa` (SoA layout)  

# Build & Install
git clone https://github.com/your‑username/NBodyBarnesHut_C.git
cd NBodyBarnesHut_C
make            # builds both AoS and SoA executables

# Usage example
./nbody_soa \
  --particles 100000 \
  --steps 500 \
  --theta 0.5 \
  --threads 8 \
  --output sim.dat

# Options
# --particles <N>     Number of bodies
# --steps <T>         Simulation time steps
# --theta <θ>         Opening angle threshold (0.0 – 1.0)
# --threads <P>       OpenMP thread count
# --layout [aos|soa]  Choose data layout (AoS default)
# --output <file>     Write positions/velocities at each step

# Run `./nbody_aos --help` or `./nbody_soa --help` for full list of flags.

## Performance
| Layout | Threads | Time (s) | Speedup vs. serial |
|:------:|:-------:|:--------:|:------------------:|
| AoS    | 1       | 120.3    | 1×                 |
| AoS    | 8       |  20.8    | 5.8×               |
| SoA    | 1       | 110.6    | 1×                 |
| SoA    | 8       |  17.5    | 6.3×               |

*(Measured on Intel Xeon Gold, 2.4 GHz, 32 GB RAM.)*



