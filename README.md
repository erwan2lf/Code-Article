# Chromatin Brownian Dynamics Simulation

A C implementation of overdamped Langevin (Brownian) dynamics for a coarse-grained
polymer chain, developed to model chromatin fiber dynamics as part of the study
*"Cell Cycle-Dependent Chromatin Motion: A Role for DNA Content Doubling Over
Cohesion"* (see [Citation](#citation) below).

The chain is confined in a spherical volume and evolves under harmonic spring
bonds, optional bending stiffness, Lennard-Jones excluded-volume interactions,
and thermal noise, integrated with an Euler–Maruyama scheme. Trajectories are
recorded to a compact binary format and analyzed for end-to-end distance,
center-of-mass motion, and mean squared displacement (MSD) via an FFT-based
O(T log T) algorithm (Calandrini / Kneller).

## Repository layout

```
.
├── main.c                  Program entry point / simulation driver
├── config.c/.h             CLI argument parsing and simulation parameters
├── simulation.c/.h         Main simulation loop
├── movement.c/.h           Brownian dynamics integrator
├── forces.c/.h             Spring, bending, confinement, Lennard-Jones forces
├── potentials.c/.h         Pairwise potential definitions
├── neighborlist.c/.h       Verlet neighbor list for LJ interactions
├── initial_structures.c/.h Initial polymer configurations (straight, solenoid,
│                           fractal globule, knot)
├── basic_functions.c/.h    Vector math, measurements, plotting utilities
├── file.c/.h                Output file management
├── traj_binary.c/.h        Binary trajectory writer/reader
├── msd.c/.h                 MSD computation from binary trajectory
├── vmd_export.c/.h         LAMMPS-format trajectory export (for VMD/OVITO)
├── mt19937ar.c/.h          Mersenne Twister PRNG
├── include/                Header files
└── plot_r2.py              Python post-processing: end-to-end distance,
                            autocorrelation, and Gaussian-chain analysis
```

## Building

Requires `gcc` (or `gcc-14`), OpenMP, and [FFTW3](http://www.fftw.org/).

```bash
gcc-14 -g -O3 -ffast-math -DMAC \
  main.c basic_functions.c simulation.c config.c file.c movement.c \
  neighborlist.c potentials.c initial_structures.c mt19937ar.c \
  traj_binary.c msd.c forces.c vmd_export.c \
  -Iinclude -I$(brew --prefix fftw)/include -L$(brew --prefix fftw)/lib \
  -lfftw3 -lm -fopenmp -o main \
  -Wall -Wno-unused-variable -Wno-unused-but-set-variable \
  -Wno-return-type -Wno-maybe-uninitialized -Wno-unused-result -Wno-comment
```

On Linux, drop `-DMAC` and point `-I`/`-L` at your system's FFTW3 installation
(e.g. via `pkg-config --cflags --libs fftw3`).

## Usage

```bash
./main <seed> <Delta> <gamma_fric> <r_conf> <T>
```

| Argument     | Description                                   |
|--------------|------------------------------------------------|
| `seed`       | RNG seed                                       |
| `Delta`      | Integration timestep                           |
| `gamma_fric` | Friction coefficient                           |
| `r_conf`     | Radius of the confining sphere                 |
| `T`          | Total number of simulation steps               |

Example:

```bash
./main 1 0.0001 1 1 5000000
```

Additional parameters (number of monomers, bending stiffness, spring
constant, recording periods, initial structure type, etc.) are set via
defaults in `config.c` and can be edited there.

Output is written to `./Results/`, including the binary trajectory
(`trajectory.bin`), a LAMMPS-format trajectory (`brownian_LJ_equilibrium.lammpstrj`)
for visualization in VMD/OVITO, end-to-end distance and MSD text files, and
center-of-mass tracking data.

## Analysis

```bash
python3 plot_r2.py Results/trajectory.bin [bead_i] [bead_j]
```

Produces end-to-end distance, autocorrelation, and Gaussian-chain analysis
plots from the recorded trajectory.

## Citation

If you use this code, please cite:

> Rey-Millet M.†, Costes L.†, Le-Floch E., Ayoub H., Saccomani Q., Manghi M.\*,
> Bystricky K.\* (2026). **Cell Cycle-Dependent Chromatin Motion: A Role for
> DNA Content Doubling Over Cohesion.** *bioRxiv*.
> https://doi.org/10.64898/2026.03.19.712877

† Joint first authors · \* Corresponding authors

**Affiliations**
1. Univ of Toulouse, CNRS, Centre de Biologie Intégrative (CBI), Toulouse, France
2. Univ Toulouse, CNRS, Laboratoire de Physique Théorique (LPT), Toulouse, France

## License

This code is released under the [MIT License](LICENSE).
