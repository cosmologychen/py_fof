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
from setuptools.dist import Distribution
import subprocess
import sys
import os
import glob
import shutil


class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        return True


def _find_system_fftw3f():
    candidate_paths = [
        "/usr/lib/x86_64-linux-gnu/libfftw3f.so",
        "/usr/lib/x86_64-linux-gnu/libfftw3f.so.3",
        "/lib/x86_64-linux-gnu/libfftw3f.so",
        "/lib/x86_64-linux-gnu/libfftw3f.so.3",
        "/usr/lib64/libfftw3f.so",
        "/usr/lib64/libfftw3f.so.3",
        "/usr/local/lib/libfftw3f.so",
        "/usr/local/lib/libfftw3f.so.3",
    ]
    for path in candidate_paths:
        if os.path.exists(path):
            return path
    try:
        out = subprocess.check_output(['ldconfig', '-p'], text=True, stderr=subprocess.DEVNULL)
        for line in out.splitlines():
            if 'libfftw3f.so' in line and '=>' in line:
                p = line.split('=>')[-1].strip()
                if os.path.exists(p):
                    return p
    except Exception:
        pass
    return None


class F2pyBuildExt(build_ext):
    """
    自定义构建命令，使用f2py混合编译Fortran/C代码
    """
    def run(self):
        # 编译Fortran扩展
        self._compile_fortran()
        # 继续执行标准的构建过程
        super().run()
        # 把扩展复制进 build_lib：bdist_wheel 先跑 build_py 再跑 build_ext，
        # 若不在此处补拷，生成的 wheel 内将不含 fof_fortran*.so。
        self._copy_extension_to_build_lib()

    def _copy_extension_to_build_lib(self):
        """将已编译的扩展复制到 build_lib/py_fof，确保 wheel 包含二进制模块。"""
        build_lib = getattr(self, 'build_lib', None)
        if not build_lib:
            return
        src_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'py_fof')
        dest_dir = os.path.join(build_lib, 'py_fof')
        os.makedirs(dest_dir, exist_ok=True)
        for so in glob.glob(os.path.join(src_dir, 'fof_fortran*.so')):
            dest = os.path.join(dest_dir, os.path.basename(so))
            if (not os.path.exists(dest)
                    or os.path.getmtime(so) > os.path.getmtime(dest)):
                shutil.copy2(so, dest)
        print(f"Copied extension(s) into {dest_dir}")
    
    def _compile_fortran(self):
        """使用f2py与gcc/gfortran混合编译Fortran/C代码"""
        package_dir = os.path.dirname(os.path.abspath(__file__))
        src_fof = os.path.join(package_dir, 'fof_fortran.f90')
        src_bound = os.path.join(package_dir, 'bound_fortran.f90')
        src_c = os.path.join(package_dir, 'atomic_dsu.c')
        output_dir = os.path.join(package_dir, 'py_fof')
        
        for src_path in (src_fof, src_bound, src_c):
            if not os.path.exists(src_path):
                raise FileNotFoundError(
                    f"Required source file not found: {src_path}"
                )
        
        max_src_mtime = max(os.path.getmtime(p) for p in (src_fof, src_bound, src_c))
        
        # 检查是否已经存在且比源文件更新的模块
        so_pattern = os.path.join(output_dir, 'fof_fortran*.so')
        so_files = glob.glob(so_pattern)
        if so_files and all(os.path.getmtime(so) >= max_src_mtime for so in so_files):
            print("Fortran module already compiled and up-to-date.")
            return
        
        # 清除旧的 .so
        for so in so_files:
            try:
                os.remove(so)
            except OSError:
                pass
        
        print("="*60)
        print("Compiling Fortran & C extension with f2py...")
        print(f"  Sources: {src_fof}, {src_bound}, {src_c}")
        print(f"  Output: {output_dir}")
        print("="*60)
        
        build_dir = os.path.join(package_dir, '.test_ai', 'f2py_build')
        os.makedirs(build_dir, exist_ok=True)
        
        build_fftw = os.path.join(build_dir, 'libfftw3f.so')
        if not os.path.exists(build_fftw):
            sys_fftw = _find_system_fftw3f()
            if sys_fftw:
                try:
                    os.symlink(sys_fftw, build_fftw)
                except OSError:
                    pass
        
        env = os.environ.copy()
        env['CC'] = 'gcc'
        env['FC'] = 'gfortran'
        
        try:
            # 1. 编译 atomic_dsu.c
            cmd_c = [
                'gcc', '-O3', '-fPIC', '-c',
                src_c, '-o', os.path.join(build_dir, 'atomic_dsu.o')
            ]
            subprocess.run(cmd_c, cwd=build_dir, capture_output=True, text=True, check=True, env=env)
            
            # 2. 编译 bound_fortran.f90
            cmd_f90 = [
                'gfortran', '-O3', '-fopenmp', '-ffree-line-length-none', '-fPIC', '-c',
                src_bound, f'-J{build_dir}', '-o', os.path.join(build_dir, 'bound_fortran.o')
            ]
            subprocess.run(cmd_f90, cwd=build_dir, capture_output=True, text=True, check=True, env=env)
            
            # 3. 打包静态库 libbound_fof.a
            lib_a = os.path.join(build_dir, 'libbound_fof.a')
            cmd_ar = [
                'ar', 'rcs', lib_a,
                os.path.join(build_dir, 'bound_fortran.o'),
                os.path.join(build_dir, 'atomic_dsu.o')
            ]
            subprocess.run(cmd_ar, cwd=build_dir, capture_output=True, text=True, check=True, env=env)
            
            # 4. 调用 f2py 在 build_dir 中编译
            cmd_f2py = [
                sys.executable, '-m', 'numpy.f2py',
                '-c', src_fof,
                '-m', 'fof_fortran',
                '--f90flags=-O3 -fopenmp -ffree-line-length-none',
                f'-I{build_dir}',
                f'-L{build_dir}',
                '-lbound_fof', '-lfftw3f', '-lgomp'
            ]
            subprocess.run(cmd_f2py, cwd=build_dir, capture_output=True, text=True, check=True, env=env)
            
            # 5. 移动生成的 .so 到 output_dir
            built_sos = glob.glob(os.path.join(build_dir, 'fof_fortran*.so'))
            if not built_sos:
                raise RuntimeError("No fof_fortran*.so produced in build directory.")
            for bso in built_sos:
                dest = os.path.join(output_dir, os.path.basename(bso))
                if os.path.exists(dest):
                    os.remove(dest)
                os.replace(bso, dest)
                
            print("Compilation successful!")
        except subprocess.CalledProcessError as e:
            print("="*60)
            print("ERROR: Compilation failed!")
            print("="*60)
            print(f"Command: {' '.join(e.cmd) if isinstance(e.cmd, list) else e.cmd}")
            print(f"\nError output:\n{e.stderr}")
            raise RuntimeError(
                "Failed to compile Fortran/C code. "
                "Please make sure you have gcc, gfortran and libfftw3f installed.\n"
                "For example:\n"
                "  Ubuntu/Debian: sudo apt-get install gcc gfortran libfftw3-dev\n"
                "  CentOS/RHEL: sudo yum install gcc gcc-gfortran fftw-devel"
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
    distclass=BinaryDistribution,
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
