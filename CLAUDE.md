# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

Configuration is in `common.mk`. Default compiler is `gfortran`, platform is `wsl`.

```bash
make Release    # Optimized build (-O3)
make Debug      # Debug build with bounds checking and profiling
make clean      # Remove build artifacts
make install    # Install to $HOME/bin (after building)
make sharedLib  # Build libOceanWave3D.so instead of executable
```

**WSL prerequisites:**
```bash
sudo apt install libhdf5-dev libhdf5-fortran-dev libopenblas-dev
```
SPARSKIT2 (`libskit`) must be compiled from source and installed to `$HOME/lib/`.

To switch compiler or platform, edit `common.mk`:
- `FC = gfortran` or `FC = ifort`
- `PLATFORM = wsl`, `gbar`, or `jess`

## Running the Solver

```bash
cd /some/working/directory
cp /path/to/OceanWave3D-Fortran90/examples/inputfiles/OceanWave3D.inp ./
$HOME/bin/OceanWave3D-Fortran OceanWave3D.inp
```

Example input files for different scenarios are in `examples/inputfiles/` (stream function waves, Berkhoff 3D, wave breaking, random waves, etc.).

Key output files:
- `LOG.txt` — convergence and progress log
- `OceanWave3D.end` — final state; rename to `OceanWave3D.init` for hot restart
- Kinematics binary files (`01`–`20`) — wave field data at sampled locations

## Architecture

OceanWave3D solves the potential flow (Laplace) equation with nonlinear free-surface boundary conditions on 3D Cartesian or curvilinear grids.

**Execution flow:**
1. `src/main/OceanWave3D.f90` — entry point; calls setup then time-marches
2. `src/main/OceanWave3DT0Setup.f90` — reads input, builds grid, allocates arrays, sets ICs
3. `src/main/OceanWave3DTakeATimeStep.f90` — one time step: free-surface update → linear solve → smoothing → output

**Module dependency layers** (bottom to top):

| Layer | Directory | Role |
|-------|-----------|------|
| Foundation | `src/variabledefs/` | Precision, constants, derived types (`datatypes.f90`), all global variables (`globalvariables.f90`) |
| Grid & Init | `src/initialization/` | Grid setup, stencil generation, relaxation zones, initial conditions |
| Linear algebra | `src/multigrid/`, `src/iterative/` | Multigrid preconditioner (V-cycle, Gauss-Seidel smoothing), GMRES solver |
| PDE assembly | `src/functions/` | Build Laplace matrix (`BuildLinearSystem.f90`), stencil operators |
| Time integration | `src/timeintegration/` | RK4 and low-storage RK45 schemes, free-surface BC evaluation |
| Physics | `src/analyticalsolutions/`, `src/pressure/`, `src/curvilinear/` | Exact wave solutions, dynamic pressure, sigma-coordinate transforms |
| I/O | `src/IO/` | HDF5 (`h5_*.f90`), ASCII, VTK output; input file parser |
| Utilities | `src/utilities/`, `src/fft/`, `src/wrappers/` | Random waves, JONSWAP/PM spectra, FFT, external library wrappers |
| Coupling | `src/OpenFoam/` | Shared-library interface for OpenFOAM coupling |

**Key global state** is in `src/variabledefs/globalvariables.f90` — all solution arrays (free-surface elevation `eta`, velocity potential `phi`), solver parameters, grid data, and I/O settings live here as module-level variables.

**Each source subdirectory** contains a `makefile.inc` that lists its `.f90`/`.f` files; the top-level `makefile` includes them all.

**MATLAB post-processing utilities** are in `utils/` (read kinematics, check against analytical solutions, visualize results).
