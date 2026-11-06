# Installation Guide

## Quick Install

### Method 1: Automatic Build & Install (Recommended)

```bash
# Clone the repository
git clone https://github.com/yourusername/py_fof.git
cd py_fof

# Install dependencies
pip install numpy matplotlib

# Install package (will automatically compile Fortran & C extensions)
pip install .
```

You can also install in editable mode for development:
```bash
pip install -e .
```

### Method 2: Import-time Auto Compilation

If the package directory is directly in your `PYTHONPATH`, importing `py_fof` will automatically detect and compile the extensions into `py_fof/`:

```bash
export PYTHONPATH="${PYTHONPATH}:$(pwd)"
python -c "import py_fof"
```

### Method 3: Manual Compilation

If you prefer to compile manually using `gcc`, `gfortran`, and `numpy.f2py`:

```bash
# Create temporary build directory
mkdir -p .test_ai/f2py_build

# 1. Compile C atomic DSU
gcc -O3 -fPIC -c atomic_dsu.c -o .test_ai/f2py_build/atomic_dsu.o

# 2. Compile Fortran adaptive bound solver
gfortran -O3 -fopenmp -fPIC -c bound_fortran.f90 -J.test_ai/f2py_build -o .test_ai/f2py_build/bound_fortran.o

# 3. Create static archive
ar rcs .test_ai/f2py_build/libbound_fof.a .test_ai/f2py_build/bound_fortran.o .test_ai/f2py_build/atomic_dsu.o

# 4. Compile FOF f2py module
cd .test_ai/f2py_build
python -m numpy.f2py -c ../../fof_fortran.f90 -m fof_fortran \
    --f90flags='-O3 -fopenmp' \
    -I. -L. -lbound_fof -lfftw3f -lgomp

# 5. Move compiled extension to package directory
mv fof_fortran*.so ../../py_fof/
cd ../..
```

## Requirements

### System Requirements
- Python >= 3.8
- C compiler (`gcc`)
- Fortran compiler (`gfortran` recommended, with OpenMP support)
- FFTW3 single-precision library (`libfftw3f`)

### Python Dependencies
```
numpy >= 1.20
matplotlib >= 3.0
```

Install Python dependencies:
```bash
pip install numpy matplotlib
```

### System Libraries & Compiler Installation

**Ubuntu / Debian:**
```bash
sudo apt-get update
sudo apt-get install gcc gfortran libfftw3-dev
```

**CentOS / RHEL / Fedora:**
```bash
sudo yum install gcc gcc-gfortran fftw-devel
```

**macOS:**
```bash
brew install gcc fftw
```

## Verification

Test the installation:

```bash
python -c "from py_fof import fof; print('Installation successful!')"
```

Run basic usage examples:

```bash
python examples/basic_usage.py
```

## Troubleshooting

### ImportError: libfftw3f.so.3: cannot open shared object file
Make sure FFTW3 single-precision library is installed (`libfftw3-dev` on Debian/Ubuntu or `fftw-devel` on RHEL/CentOS).

### Compilation Error: "f2py not found"
Ensure NumPy is installed in the active Python environment:
```bash
pip install numpy
```

## Uninstallation

```bash
pip uninstall py_fof
```
