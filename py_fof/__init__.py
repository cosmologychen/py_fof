"""
CUBE FoF (Friends-of-Friends) Halo Finder

Python interface for the Fortran FOF implementation.
Provides a user-friendly API for identifying halos in cosmological simulations.

Author: [Your Name]
Date: 2024

Usage:
    from py_fof import fof
    
    # Simple usage: single halo
    xp = np.random.rand(3, 1000) * 512  # 1000 particles in a 512^3 box
    in_halo = fof(xp)
    
    # Multiple halos with velocity-based binding check
    xp_list = [xp1, xp2, xp3]  # list of particle arrays
    vp_list = [vp1, vp2, vp3]  # corresponding velocities
    in_halo = fof(xp_list, vp_list=vp_list, bound_check=True)
"""

import numpy as np
import os
import sys
import subprocess
import glob


def _compile_fortran():
    """
    自动编译Fortran代码
    
    如果编译后的模块不存在，自动调用f2py编译
    """
    package_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(package_dir)
    
    fortran_source = os.path.join(parent_dir, 'fof_fortran.f90')
    
    so_pattern = os.path.join(package_dir, 'fof_fortran*.so')
    so_files = glob.glob(so_pattern)
    
    if so_files:
        return
    
    if not os.path.exists(fortran_source):
        raise FileNotFoundError(
            f"Fortran source file not found: {fortran_source}\n"
            "Please make sure fof_fortran.f90 is in the package directory."
        )
    
    print(f"Compiling Fortran code...")
    print(f"  Source: {fortran_source}")
    print(f"  Output: {package_dir}")
    
    cmd = [
        sys.executable, '-m', 'numpy.f2py',
        '-c', fortran_source,
        '-m', 'fof_fortran',
        '--f90flags=-O3'
    ]
    
    try:
        result = subprocess.run(
            cmd,
            cwd=package_dir,
            capture_output=True,
            text=True,
            check=True
        )
        print("Compilation successful!")
    except subprocess.CalledProcessError as e:
        print(f"Compilation failed!")
        print(f"Error output:\n{e.stderr}")
        raise RuntimeError(
            "Failed to compile Fortran code. "
            "Please compile manually with:\n"
            f"  cd {parent_dir}\n"
            f"  python -m numpy.f2py -c fof_fortran.f90 -m fof_fortran --f90flags='-O3'\n"
            f"  mv fof_fortran*.so py_fof/"
        )


_compile_fortran()

from . import fof_fortran


def fof(xp_list, vp_list=None, L_Box=512.0, linking_length=0.2, 
        np_halo_min=10, bound_check=False, verbose=0, return_ang_mom=False):
    """
    Friends-of-Friends (FOF) halo finder.
    
    Parameters
    ----------
    xp_list : array_like or list of array_like
        粒子坐标。可以是单个数组 shape=(3, N)，或者是多个halo的列表
        [xp_halo1, xp_halo2, ...]，每个元素 shape=(3, Ni)
        单位：格点
        
    vp_list : array_like or list of array_like, optional
        粒子速度。格式与xp_list相同。
        仅当 bound_check=True 时需要。
        单位：格点/时间步
        
    L_Box : float, default=512.0
        盒子边长，用于周期性边界条件。
        单位：格点
        
    linking_length : float, default=0.2
        FOF链接长度。两个粒子距离小于此值时被链接。
        单位：格点
        
    np_halo_min : int, default=10
        最小halo粒子数阈值。粒子数少于此值的FOF链不被识别为halo。
        
    bound_check : bool, default=False
        是否执行束缚检查。如果为True，需要提供vp_list。
        束缚检查会迭代移除不满足引力束缚条件的粒子（最多3轮）。
        
    verbose : int, default=0
        详细输出级别：
        - 0: 关闭详细输出
        - 1: 输出halo处理进度和bound_check信息
        - 2: 输出更详细的信息（包括网格大小）
        
    Returns
    -------
    in_halo : ndarray
        逻辑数组，shape=(total_np,)，标记每个粒子是否属于有效halo。
        
    Examples
    --------
    >>> import numpy as np
    >>> from py_fof import fof
    >>> 
    >>> # 单个halo
    >>> xp = np.random.rand(3, 1000).astype(np.float64) * 512
    >>> in_halo = fof(xp, L_Box=512, linking_length=0.2)
    >>> print(f"Found {in_halo.sum()} particles in halos")
    >>> 
    >>> # 多个halo
    >>> xp1 = np.random.rand(3, 500).astype(np.float64) * 100 + 200
    >>> xp2 = np.random.rand(3, 800).astype(np.float64) * 100 + 350
    >>> in_halo = fof([xp1, xp2], L_Box=512)
    >>> 
    >>> # 带束缚检查
    >>> vp1 = np.random.randn(3, 500).astype(np.float64) * 0.1
    >>> vp2 = np.random.randn(3, 800).astype(np.float64) * 0.1
    >>> in_halo = fof([xp1, xp2], vp_list=[vp1, vp2], bound_check=True)
    >>> 
    >>> # 带详细输出
    >>> in_halo = fof([xp1, xp2], vp_list=[vp1, vp2], bound_check=True, verbose=1)
    
    Notes
    -----
    1. 所有坐标和长度参数均以格点为单位。
    2. 束缚检查条件: 0.5 * v_rel^2 < np_fof / r
       这是一个简化的引力束缚判据。
    3. 束缚检查最多迭代3轮，或直到所有粒子都束缚。
    """
    
    if bound_check and vp_list is None:
        raise ValueError("vp_list is required when bound_check=True")
    
    if isinstance(xp_list, np.ndarray) and xp_list.ndim == 2:
        xp_list = [xp_list]
        single_halo = True
    else:
        single_halo = False
    
    nhalo = len(xp_list)
    np_halo = np.array([xp.shape[1] for xp in xp_list], dtype=np.int32)
    total_np = int(np.sum(np_halo))
    
    xp_halo = np.hstack(xp_list).astype(np.float64, order='F')
    
    if xp_halo.shape[0] != 3:
        raise ValueError(f"xp_halo must have shape (3, N), got {xp_halo.shape}")
    
    if vp_list is not None:
        if isinstance(vp_list, np.ndarray) and vp_list.ndim == 2:
            vp_list = [vp_list]
        if len(vp_list) != nhalo:
            raise ValueError(f"vp_list length ({len(vp_list)}) must match xp_list length ({nhalo})")
        vp_halo = np.hstack(vp_list).astype(np.float64, order='F')
        if vp_halo.shape != xp_halo.shape:
            raise ValueError(f"vp_halo shape {vp_halo.shape} must match xp_halo shape {xp_halo.shape}")
    else:
        vp_halo = np.zeros((3, total_np), dtype=np.float64, order='F')
    
    in_halo, ang_mom = fof_fortran.fof_fortran(
        xp_halo, vp_halo, np_halo, nhalo,
        float(L_Box), float(linking_length), 
        int(np_halo_min), int(bound_check), int(verbose), total_np
    )
    
    if return_ang_mom:
        return in_halo.astype(bool), ang_mom
    else:
        return in_halo.astype(bool)


def fof_with_info(xp_list, vp_list=None, L_Box=512.0, linking_length=0.2, 
                  np_halo_min=10, bound_check=False, verbose=0, return_ang_mom=False):
    """
    FOF halo finder with detailed output.
    
    Same as fof(), but returns additional information about each halo.
    
    Returns
    -------
    result : dict
        Dictionary containing:
        - 'in_halo': logical array marking particles in halos
        - 'nhalo_found': number of halos found
        - 'np_per_halo': array of particle counts per halo
        - 'halo_indices': list of arrays, each containing particle indices for a halo
    """
    if return_ang_mom:
        in_halo, ang_mom = fof(xp_list, vp_list, L_Box, linking_length, np_halo_min, bound_check, verbose, return_ang_mom=True)
    else:
        in_halo = fof(xp_list, vp_list, L_Box, linking_length, np_halo_min, bound_check, verbose, return_ang_mom=False)
    
    if isinstance(xp_list, np.ndarray) and xp_list.ndim == 2:
        xp_concat = xp_list
    else:
        xp_concat = np.hstack(xp_list)
    
    halo_indices = []
    np_per_halo = []
    
    visited = np.zeros(len(in_halo), dtype=bool)
    
    for i in range(len(in_halo)):
        if in_halo[i] and not visited[i]:
            halo_mask = _find_connected_particles(i, in_halo, xp_concat, linking_length, L_Box)
            halo_idx = np.where(halo_mask)[0]
            halo_indices.append(halo_idx)
            np_per_halo.append(len(halo_idx))
            visited[halo_idx] = True
    
    res = {
        'in_halo': in_halo,
        'nhalo_found': len(halo_indices),
        'np_per_halo': np.array(np_per_halo),
        'halo_indices': halo_indices
    }
    if return_ang_mom:
        res['ang_mom'] = ang_mom
    return res


def _find_connected_particles(start_idx, in_halo, xp, linking_length, L_Box):
    """Find all particles connected to start_idx using BFS."""
    from collections import deque
    
    n_particles = xp.shape[1]
    visited = np.zeros(n_particles, dtype=bool)
    halo_mask = np.zeros(n_particles, dtype=bool)
    
    queue = deque([start_idx])
    visited[start_idx] = True
    L_Box2 = L_Box / 2.0
    
    while queue:
        i = queue.popleft()
        halo_mask[i] = True
        
        for j in range(n_particles):
            if not visited[j] and in_halo[j]:
                dx = xp[:, i] - xp[:, j]
                dx = np.mod(dx + L_Box2, L_Box) - L_Box2
                if np.sum(dx**2) <= linking_length**2:
                    visited[j] = True
                    queue.append(j)
    
    return halo_mask
