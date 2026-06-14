# CUBE FoF (Friends-of-Friends) Halo Finder

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Python interface for the Fortran FOF implementation, providing a user-friendly API for identifying halos in cosmological simulations.

## Features

- **Fast FOF Algorithm**: Core implementation in Fortran with OpenMP parallelization
- **Python Interface**: Easy-to-use Python API with NumPy arrays
- **Binding Check**: Optional velocity-based gravitational binding check
- **Multiple Halo Support**: Handle multiple halos in a single call
- **Periodic Boundary Conditions**: Automatic handling of periodic boundaries
- **Visualization Tools**: Built-in plotting functions for result visualization

## Installation

### Requirements

- Python >= 3.8
- NumPy >= 1.20
- Matplotlib >= 3.0
- A Fortran compiler (gfortran recommended)

### Install from Source

```bash
git clone https://github.com/yourusername/py_fof.git
cd py_fof
pip install .
```

The Fortran code will be automatically compiled during installation.

### Manual Compilation (Optional)

If automatic compilation fails, you can manually compile the Fortran code:

```bash
python -m numpy.f2py -c fof_fortran.f90 -m fof_fortran --f90flags='-O3 -fopenmp'
mv fof_fortran*.so py_fof/
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

### `fof(xp_list, vp_list=None, L_Box=512.0, linking_length=0.2, np_halo_min=10, bound_check=False, verbose=0)`

Friends-of-Friends halo finder.

**Parameters:**

- `xp_list` (array_like or list of array_like): Particle coordinates. Can be a single array with shape (3, N) or a list of arrays for multiple halos. Units: grid points.
- `vp_list` (array_like or list of array_like, optional): Particle velocities with shape (3, N) or list of arrays. Required for binding check.
- `L_Box` (float, optional): Box size. Default: 512.0
- `linking_length` (float, optional): FOF linking length. Default: 0.2
- `np_halo_min` (int, optional): Minimum number of particles to form a halo. Default: 10
- `bound_check` (bool, optional): Enable velocity-based binding check. Default: False
- `verbose` (int, optional): Verbosity level (0=quiet, 1=progress, 2=detailed). Default: 0

**Returns:**

- `in_halo` (numpy.ndarray): Boolean array indicating whether each particle belongs to a halo.

## Algorithm

The FOF algorithm identifies halos by linking particles that are within a specified linking length of each other. The algorithm consists of:

1. **Grid Construction**: Particles are binned into a 3D grid for efficient neighbor search
2. **FOF Search**: Neighboring particles within the linking length are linked
3. **Chain Merging**: Linked particles are grouped into halos using a union-find algorithm
4. **Binding Check** (optional): Particles are checked for gravitational binding using:
   ```
   0.5 * v_rel² < N_particles / r
   ```
   where v_rel is the relative velocity and r is the distance from the halo center.

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
