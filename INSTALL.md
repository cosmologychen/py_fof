# Installation Guide

## Quick Install

### Method 1: Direct Installation (Recommended)

```bash
# Clone the repository
git clone https://github.com/yourusername/py_fof.git
cd py_fof

# Install dependencies
pip install numpy matplotlib

# The Fortran module is already compiled and included
# Just add the package to your Python path
export PYTHONPATH="${PYTHONPATH}:$(pwd)"
```

### Method 2: Build from Source

If you need to recompile the Fortran code:

```bash
# Install dependencies
pip install numpy matplotlib

# Compile Fortran code
python -m numpy.f2py -c fof_fortran.f90 -m fof_fortran --f90flags='-O3 -fopenmp'

# Move compiled module to package directory
mv fof_fortran*.so py_fof/

# Add to Python path
export PYTHONPATH="${PYTHONPATH}:$(pwd)"
```

### Method 3: Using setup.py (for distribution)

```bash
# Build source distribution
python setup.py sdist

# The package will be in dist/py_fof-1.0.0.tar.gz
# Users can install with:
pip install py_fof-1.0.0.tar.gz
```

## Requirements

### System Requirements
- Python >= 3.8
- Fortran compiler (gfortran recommended) - only needed if recompiling

### Python Dependencies
```
numpy >= 1.20
matplotlib >= 3.0
```

Install dependencies:
```bash
pip install numpy matplotlib
```

### Fortran Compiler Installation

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install gfortran
```

**macOS:**
```bash
brew install gcc
```

**CentOS/RHEL:**
```bash
sudo yum install gcc-gfortran
```

## Verification

Test the installation:

```bash
python -c "from py_fof import fof; print('Installation successful!')"
```

Run the test suite:

```bash
python test_fof.py
```

## Troubleshooting

### ImportError: No module named 'fof_fortran'

The Fortran module is not compiled or not in the correct location.

**Solution:**
```bash
cd /path/to/py_fof
python -m numpy.f2py -c fof_fortran.f90 -m fof_fortran --f90flags='-O3'
mv fof_fortran*.so py_fof/
```

### Compilation Error: "f2py not found"

NumPy is not installed or f2py is not in PATH.

**Solution:**
```bash
pip install numpy
# or
conda install numpy
```

### RuntimeError: "Fortran source file not found"

The working directory is incorrect.

**Solution:**
Make sure you're in the py_fof directory when running the code.

## Development Installation

For development, install in editable mode:

```bash
git clone https://github.com/yourusername/py_fof.git
cd py_fof
pip install -e .
```

## Uninstallation

```bash
pip uninstall py_fof
```
