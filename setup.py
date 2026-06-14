"""
CUBE FoF (Friends-of-Friends) Halo Finder - Python Package

Python interface for the Fortran FOF implementation.
Provides a user-friendly API for identifying halos in cosmological simulations.

Installation:
    pip install .

Requirements:
    - Python >= 3.8
    - NumPy >= 1.20
    - A Fortran compiler (gfortran or ifort)
    - numpy.f2py for building the Fortran extension
"""

from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext
import subprocess
import sys
import os
import glob


class F2pyBuildExt(build_ext):
    """
    自定义构建命令，使用f2py编译Fortran代码
    """
    def run(self):
        # 编译Fortran扩展
        self._compile_fortran()
        # 继续执行标准的构建过程
        super().run()
    
    def _compile_fortran(self):
        """使用f2py编译Fortran代码"""
        package_dir = os.path.dirname(os.path.abspath(__file__))
        fortran_source = os.path.join(package_dir, 'fof_fortran.f90')
        output_dir = os.path.join(package_dir, 'py_fof')
        
        if not os.path.exists(fortran_source):
            raise FileNotFoundError(
                f"Fortran source file not found: {fortran_source}"
            )
        
        # 检查是否已经存在编译好的模块
        so_pattern = os.path.join(output_dir, 'fof_fortran*.so')
        if glob.glob(so_pattern):
            print("Fortran module already compiled.")
            return
        
        print("="*60)
        print("Compiling Fortran extension with f2py...")
        print(f"  Source: {fortran_source}")
        print(f"  Output: {output_dir}")
        print("="*60)
        
        cmd = [
            sys.executable, '-m', 'numpy.f2py',
            '-c', fortran_source,
            '-m', 'fof_fortran',
            '--f90flags=-O3 -fopenmp'
        ]
        
        try:
            result = subprocess.run(
                cmd,
                cwd=output_dir,
                capture_output=True,
                text=True,
                check=True
            )
            print("Compilation successful!")
        except subprocess.CalledProcessError as e:
            print("="*60)
            print("ERROR: Compilation failed!")
            print("="*60)
            print(f"Command: {' '.join(cmd)}")
            print(f"\nError output:\n{e.stderr}")
            raise RuntimeError(
                "Failed to compile Fortran code. "
                "Please make sure you have a Fortran compiler installed.\n"
                "For example:\n"
                "  Ubuntu/Debian: sudo apt-get install gfortran\n"
                "  macOS: brew install gcc\n"
                "  CentOS/RHEL: sudo yum install gcc-gfortran"
            )


# 读取README文件
readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
long_description = ""
if os.path.exists(readme_path):
    with open(readme_path, 'r', encoding='utf-8') as f:
        long_description = f.read()


setup(
    name='py_fof',
    version='1.0.0',
    description='Friends-of-Friends (FOF) Halo Finder for Cosmological Simulations',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Your Name',
    author_email='your.email@example.com',
    url='https://github.com/yourusername/py_fof',
    packages=find_packages(),
    package_data={
        'py_fof': ['*.so', '*.f90'],
    },
    include_package_data=True,
    cmdclass={
        'build_ext': F2pyBuildExt,
    },
    install_requires=[
        'numpy>=1.20',
        'matplotlib>=3.0',
    ],
    extras_require={
        'dev': [
            'pytest>=6.0',
            'pytest-cov',
            'flake8',
            'black',
        ],
    },
    python_requires='>=3.8',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Fortran',
    ],
    keywords='cosmology halo fof friends-of-friends n-body simulation',
    project_urls={
        'Bug Reports': 'https://github.com/yourusername/py_fof/issues',
        'Source': 'https://github.com/yourusername/py_fof',
    },
)
