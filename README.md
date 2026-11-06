# CUBE FoF (Friends-of-Friends) Halo Finder

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Python interface for the Fortran FOF implementation, providing a user-friendly API for identifying halos in cosmological simulations.

## Features

- **High-Performance FOF**: Core algorithm based on 2-pass spatial CRS + 13 half-neighborhood search + lock-free atomic Disjoint-Set Union (DSU) in Fortran/C with OpenMP multi-threading.
- **Adaptive Gravitational Unbinding**: Multi-resolution isolated PM (48~256 grids) + C2 short-range continuous softening residuals and analytic direct summation, iteratively removing unbound particles in 2% energy batches.
- **Dual Unit Modes & Variable Mass**: Seamless support for CUBE grid units (`units='cube'`) and astrophysical physical units (`units='physical'`: Mpc/h, km/s, \(M_\odot/h\)), supporting equal and per-particle variable masses.
- **Pythonic API**: Single or multi-halo inputs, periodic boundary handling, $O(N)$ linear catalog generation, and angular momentum calculation.
- **Automated Mixed Build**: On-demand timestamp checking and f2py/gcc/gfortran compilation into `py_fof/`.

## Installation

### Requirements

- Python >= 3.8
- NumPy >= 1.20
- Matplotlib >= 3.0
- C compiler (`gcc`)
- Fortran compiler (`gfortran` recommended, with OpenMP)
- FFTW3 single-precision library (`libfftw3f`)

### Install from Source

```bash
git clone https://github.com/yourusername/py_fof.git
cd py_fof
pip install .
```

The Fortran & C extensions will be automatically compiled during installation or upon first import.

### Manual Compilation (Optional)

See [INSTALL.md](INSTALL.md) for detailed manual compilation instructions. Alternatively:

```bash
python -c "import py_fof"
```

## Quick Start

### Basic Usage

```python
import numpy as np
from py_fof import fof

# Generate test particles (single halo)
np.random.seed(42)
xp = np.random.rand(3, 1000).astype(np.float64) * 512  # 1000 particles in 512^3 box

# Run FOF
in_halo = fof(xp, L_Box=512.0, linking_length=0.2, np_halo_min=10)

print(f"Particles in halo: {np.sum(in_halo)} / {len(in_halo)}")
```

### Multiple Halos with Binding Check

```python
# Generate multiple halos
xp1 = np.random.rand(3, 500).astype(np.float64) * 100 + 50
xp2 = np.random.rand(3, 500).astype(np.float64) * 100 + 200
xp3 = np.random.rand(3, 500).astype(np.float64) * 100 + 350

xp_list = [xp1, xp2, xp3]

# Generate velocities
vp1 = np.random.randn(3, 500).astype(np.float64) * 0.5
vp2 = np.random.randn(3, 500).astype(np.float64) * 0.5
vp3 = np.random.randn(3, 500).astype(np.float64) * 0.5

vp_list = [vp1, vp2, vp3]

# Run FOF with binding check
in_halo = fof(
    xp_list, 
    vp_list=vp_list,
    L_Box=512.0, 
    linking_length=0.2, 
    np_halo_min=10,
    bound_check=True,
    verbose=1
)

print(f"Total particles in halos: {np.sum(in_halo)}")
```

### Visualization

```python
from py_fof.test_fof import plot_particles

# Plot results
plot_particles(
    np.hstack(xp_list), 
    in_halo, 
    title="FOF Halo Finder Result",
    save_path="halo_result.png"
)
```

## API Reference

### `fof(...)`
```python
def fof(xp_list, vp_list=None, L_Box=512.0, linking_length=0.2,
        np_halo_min=10, bound_check=False, verbose=0, return_ang_mom=False,
        n_threads=None, mp_list=None, softening=None, z=0.0,
        box=300.0, h0=0.67, omega_m=0.3, omega_l=0.7, n_sim=None,
        units='cube', include_hubble_flow=False):
```

**Parameters:**

- `xp_list` (array_like or list of array_like): Particle coordinates. Shape (3, N) or list of (3, Ni) arrays. Units depend on `units`.
- `vp_list` (array_like or list of array_like, optional): Particle velocities. Required when `bound_check=True`.
- `L_Box` (float, optional): Simulation box size for periodic boundary conditions. Default: 512.0.
- `linking_length` (float, optional): FOF linking length. Default: 0.2.
- `np_halo_min` (int, optional): Minimum number of particles to form a valid halo. Default: 10.
- `bound_check` (bool, optional): Enable self-consistent gravitational binding check. Default: False.
- `verbose` (int, optional): Verbosity level (0=quiet, 1=progress, 2=detailed). Default: 0.
- `return_ang_mom` (bool, optional): Whether to return accumulated angular momentum. Default: False.
- `n_threads` (int, optional): Max OpenMP threads. Default None uses all available cores (up to 256).
- `mp_list` (scalar, array, or list of arrays, optional): Particle mass. Dimensionless grid weights in `units='cube'`, or \(M_\odot/h\) in `units='physical'`.
- `softening` (float, optional): Gravitational softening length. Default 0.1 grid units in `units='cube'`; must be explicitly provided in `units='physical'` (comoving Mpc/h).
- `z` (float, optional): Cosmological redshift. Default: 0.0.
- `units` (`'cube'` or `'physical'`, optional): Unit system. Default `'cube'`.
- `include_hubble_flow` (bool, optional): Whether to include relative Hubble expansion flow in binding kinetic energy. Default: False.

**Returns:**
- `in_halo` (numpy.ndarray): 1D boolean array indicating bound halo membership.
- `ang_mom` (numpy.ndarray, optional): Accumulated angular momentum shape (3, nhalo) if `return_ang_mom=True`.

### `fof_with_info(...)`
Same arguments as `fof(...)`, returning a dictionary via $O(N)$ linear extraction:
- `'in_halo'`: boolean mask array
- `'nhalo_found'`: total number of identified halos
- `'np_per_halo'`: array of particle counts per halo
- `'halo_indices'`: list of arrays with global particle indices for each halo
- `'ang_mom'`: accumulated angular momentum (if `return_ang_mom=True`)

## Algorithm

1. **Two-Pass Spatial Compressed Row Storage (CRS)**:
   Binning particles into 3D cells using two sequential passes, eliminating linked-list pointers and achieving high cache locality.
2. **13 Half-Neighborhood & Lock-Free Atomic DSU**:
   Symmetric half-neighborhood search avoiding redundant pair calculations; lock-free atomic DSU (`atomic_dsu.c`) for multi-core graph union.
3. **Adaptive Self-Consistent Gravitational Unbinding (`bound_test`)**:
   - **Direct Summation** (\(N \le 1000\)): Continuous softening kernel with incremental subtraction of unpinned particles;
   - **Multi-Resolution Isolated PM + C2 Short-Range Residuals** (\(N > 1000\)): 48/64/128/256 isolated grids with optimal Green functions and self-potential removal, coupled with cell-CSR short-range interactions;
   - **2% Energy Batch Removal**: Iterative unbinding removing the least bound 2% particles per step up to 20 outer potential refreshes.

## Testing

Run the test suite:

```bash
python test_fof.py
```

This will generate several test images:
- `test_predictable_halo.png`: Single halo test
- `test_predictable_halo_velocity.png`: Halo with velocity binding check
- `test_multiple_predictable_halos.png`: Multiple halos test
- `test_multiple_predictable_halos_velocity.png`: Multiple halos with velocity

## Performance

The Fortran implementation with OpenMP parallelization provides excellent performance:
- ~10^5 particles: < 1 second
- ~10^6 particles: ~10 seconds
- ~10^7 particles: ~2 minutes

(Performance depends on hardware and linking length)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this code in your research, please cite:

```bibtex
@software{py_fof,
  author = {Your Name},
  title = {CUBE FoF: Friends-of-Friends Halo Finder},
  url = {https://github.com/yourusername/py_fof},
  year = {2024}
}
```

## Acknowledgments

- Original Fortran implementation based on standard FOF algorithms
- Inspired by cosmological simulation codes like GADGET and RAMSES
