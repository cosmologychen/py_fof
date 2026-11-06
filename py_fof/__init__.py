"""
CUBE FoF (Friends-of-Friends) Halo Finder

Python interface for the Fortran FOF implementation.
Provides a user-friendly API for identifying halos in cosmological simulations.
"""

import numpy as np
import os
import sys
import subprocess
import glob


def _find_system_fftw3f():
    """查找系统中存在的 libfftw3f 动态库路径"""
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

    # 尝试 ldconfig 查找
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


def _compile_fortran():
    """
    自动多文件混合编译 Fortran + C 扩展模块。
    检查源文件修改时间，在不存在或源文件更新时调用 gcc/gfortran/f2py 重新编译。
    所有中间编译产物隔离在 .test_ai/f2py_build 目录下。
    """
    package_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(package_dir)

    src_fof = os.path.join(parent_dir, 'fof_fortran.f90')
    src_bound = os.path.join(parent_dir, 'bound_fortran.f90')
    src_c = os.path.join(parent_dir, 'atomic_dsu.c')

    so_pattern = os.path.join(package_dir, 'fof_fortran*.so')
    so_files = glob.glob(so_pattern)

    srcs_exist = all(os.path.exists(p) for p in (src_fof, src_bound, src_c))
    if so_files and not srcs_exist:
        # 安装态（如 wheel 安装后无源码），直接复用已有 .so
        return

    if not srcs_exist:
        missing = [p for p in (src_fof, src_bound, src_c) if not os.path.exists(p)]
        raise FileNotFoundError(
            f"Required source files not found: {missing}\n"
            "Please make sure fof_fortran.f90, bound_fortran.f90, and atomic_dsu.c exist."
        )

    max_src_mtime = max(os.path.getmtime(p) for p in (src_fof, src_bound, src_c))
    if so_files and all(os.path.getmtime(so) >= max_src_mtime for so in so_files):
        return

    # 清理旧的 .so 文件
    for so in so_files:
        try:
            os.remove(so)
        except OSError:
            pass

    print("Compiling Fortran & C extensions with f2py...")
    print(f"  Sources: {src_fof}, {src_bound}, {src_c}")
    print(f"  Target: {package_dir}")

    build_dir = os.path.join(parent_dir, '.test_ai', 'f2py_build')
    os.makedirs(build_dir, exist_ok=True)

    # 准备 libfftw3f 软链接（若需要）
    build_fftw = os.path.join(build_dir, 'libfftw3f.so')
    if not os.path.exists(build_fftw):
        sys_fftw = _find_system_fftw3f()
        if sys_fftw:
            try:
                os.symlink(sys_fftw, build_fftw)
            except OSError:
                pass

    # 确保编译器统一为 gcc / gfortran
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

        # 4. 调用 f2py 在 build_dir 中编译 fof_fortran.f90
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

        # 5. 将生成的 .so 移动到 package_dir
        built_sos = glob.glob(os.path.join(build_dir, 'fof_fortran*.so'))
        if not built_sos:
            raise RuntimeError("No fof_fortran*.so produced in build directory.")
        for bso in built_sos:
            dest = os.path.join(package_dir, os.path.basename(bso))
            if os.path.exists(dest):
                os.remove(dest)
            os.replace(bso, dest)

        print("Compilation successful!")
    except subprocess.CalledProcessError as e:
        print("Compilation failed!")
        print(f"Command: {' '.join(e.cmd) if isinstance(e.cmd, list) else e.cmd}")
        print(f"Error output:\n{e.stderr}")
        raise RuntimeError(
            "Failed to compile Fortran/C code.\n"
            f"Error details:\n{e.stderr}"
        )


_compile_fortran()

from . import fof_fortran


def _fof_core(xp_list, vp_list=None, L_Box=512.0, linking_length=0.2,
              np_halo_min=10, bound_check=False, verbose=0, return_ang_mom=False,
              n_threads=None, mp_list=None, softening=None, z=0.0,
              box=300.0, h0=0.67, omega_m=0.3, omega_l=0.7, n_sim=None,
              units='cube', include_hubble_flow=False,
              ng=None, ng_global=None):
    """
    底层核心参数校验、单位转换与 Fortran 调用逻辑。
    返回: (in_halo, halo_id, ang_mom)
    """
    if units not in ('cube', 'physical'):
        raise ValueError(f"units must be 'cube' or 'physical', got '{units}'")

    if n_threads is None or n_threads <= 0:
        n_threads_val = 0
    else:
        if n_threads > 256:
            raise ValueError("n_threads cannot exceed 256")
        n_threads_val = int(n_threads)

    # 规范化 xp_list
    if isinstance(xp_list, np.ndarray) and xp_list.ndim == 2:
        xp_list_norm = [xp_list]
    elif isinstance(xp_list, (list, tuple)):
        xp_list_norm = list(xp_list)
    else:
        raise ValueError("xp_list must be a 2D numpy array (3, N) or a list of 2D arrays")

    nhalo = len(xp_list_norm)
    for i, xp in enumerate(xp_list_norm):
        if not isinstance(xp, np.ndarray) or xp.ndim != 2 or xp.shape[0] != 3:
            raise ValueError(f"xp_list[{i}] must have shape (3, N), got {getattr(xp, 'shape', type(xp))}")

    np_halo = np.array([xp.shape[1] for xp in xp_list_norm], dtype=np.int32)
    total_np = int(np.sum(np_halo))

    if total_np == 0:
        in_halo = np.zeros(0, dtype=bool)
        halo_id = np.zeros(0, dtype=np.int32)
        ang_mom = np.zeros((3, nhalo), dtype=np.float64)
        return in_halo, halo_id, ang_mom

    xp_halo = np.asfortranarray(np.hstack(xp_list_norm), dtype=np.float64)

    # 规范化 vp_list
    if bound_check and vp_list is None:
        raise ValueError("vp_list is required when bound_check=True")

    if vp_list is not None:
        if isinstance(vp_list, np.ndarray) and vp_list.ndim == 2:
            vp_list_norm = [vp_list]
        elif isinstance(vp_list, (list, tuple)):
            vp_list_norm = list(vp_list)
        else:
            raise ValueError("vp_list must be a 2D numpy array (3, N) or a list of 2D arrays")

        if len(vp_list_norm) != nhalo:
            raise ValueError(f"vp_list length ({len(vp_list_norm)}) must match xp_list length ({nhalo})")

        for i, vp in enumerate(vp_list_norm):
            if not isinstance(vp, np.ndarray) or vp.ndim != 2 or vp.shape != xp_list_norm[i].shape:
                raise ValueError(f"vp_list[{i}] shape {getattr(vp, 'shape', None)} must match xp_list[{i}] shape {xp_list_norm[i].shape}")

        vp_halo = np.asfortranarray(np.hstack(vp_list_norm), dtype=np.float64)
    else:
        vp_halo = np.zeros((3, total_np), dtype=np.float64, order='F')

    # 规范化 mp_list
    mp_raw = None
    if mp_list is not None:
        if np.isscalar(mp_list) or (isinstance(mp_list, np.ndarray) and mp_list.ndim == 0):
            val = float(mp_list)
            if val <= 0 or not np.isfinite(val):
                raise ValueError("mp_list scalar must be positive and finite")
            mp_raw = np.full(total_np, val, dtype=np.float64)
        elif isinstance(mp_list, np.ndarray):
            if mp_list.ndim == 1 and len(mp_list) == total_np:
                mp_raw = np.ascontiguousarray(mp_list, dtype=np.float64)
            elif nhalo > 1 and mp_list.ndim == 1 and len(mp_list) == nhalo:
                mp_raw = np.repeat(np.ascontiguousarray(mp_list, dtype=np.float64), np_halo)
            else:
                raise ValueError(f"mp_list 1D array length ({len(mp_list)}) does not match total particles ({total_np})")
        elif isinstance(mp_list, (list, tuple)):
            if len(mp_list) == nhalo and all(isinstance(m, np.ndarray) and m.ndim == 1 and len(m) == np_halo[i] for i, m in enumerate(mp_list)):
                mp_raw = np.concatenate(mp_list).astype(np.float64)
            elif len(mp_list) == nhalo and all(np.isscalar(m) for m in mp_list):
                scalars = np.array([float(m) for m in mp_list], dtype=np.float64)
                mp_raw = np.repeat(scalars, np_halo)
            elif len(mp_list) == total_np and all(np.isscalar(m) for m in mp_list):
                mp_raw = np.array(mp_list, dtype=np.float64)
            else:
                raise ValueError("mp_list structure does not match xp_list structure or particle counts")
        else:
            raise ValueError("Unsupported mp_list type")

        if not np.all(np.isfinite(mp_raw)) or not np.all(mp_raw > 0):
            raise ValueError("All particle masses in mp_list must be positive and finite")

    # 尺度因子与哈勃流参数
    z_val = float(z)
    if 1.0 + z_val <= 0.0:
        raise ValueError(f"1 + z must be positive, got {1.0 + z_val}")
    a_val = 1.0 / (1.0 + z_val)

    if include_hubble_flow:
        t0 = 1.0 / (150.0 * float(h0) * np.sqrt(float(omega_m)))
        hubble = 100.0 * float(h0) * np.sqrt(float(omega_m) / (a_val ** 3) + float(omega_l))
        b_val = hubble * (a_val ** 2) * t0
    else:
        b_val = 0.0

    # 双单位转换
    if units == 'cube':
        softening_grid = 0.1 if softening is None else float(softening)
        if bound_check and softening_grid <= 0:
            raise ValueError("softening must be positive")
        ng_eff_global = float(ng_global) if ng_global is not None else float(L_Box)
        s_15_val = 15.0 * float(h0) * ng_eff_global / float(box)
        xp_code = xp_halo
        L_Box_code = float(L_Box)
        link_code = float(linking_length)
        vp_code = vp_halo
        if mp_raw is None:
            mp_code = np.ones(total_np, dtype=np.float64)
        else:
            mp_code = mp_raw
        ang_mom_scale = 1.0
    else:  # units == 'physical'
        if bound_check:
            if softening is None or float(softening) <= 0:
                raise ValueError("softening must be explicitly provided and > 0 when units='physical' and bound_check=True")
            if mp_raw is None and n_sim is None:
                raise ValueError("Either mp_list or n_sim must be provided when units='physical' and bound_check=True")

        ng_eff = float(round(float(n_sim) ** (1.0 / 3.0))) if (n_sim is not None and float(n_sim) > 0) else (float(ng_global) if ng_global is not None else 512.0)
        box_mpc_h = float(L_Box)
        dx = box_mpc_h / ng_eff

        xp_code = np.asfortranarray(xp_halo / dx, dtype=np.float64)
        L_Box_code = ng_eff
        link_code = float(linking_length) / dx
        softening_grid = 0.1 if softening is None else (float(softening) / dx)
        s_15_val = 15.0 * float(h0) * ng_eff / box_mpc_h

        v_scale = (a_val * ng_eff) / (150.0 * np.sqrt(float(omega_m)) * box_mpc_h)
        vp_code = np.asfortranarray(vp_halo * v_scale, dtype=np.float64)

        rho_crit_0 = 2.77536627e11
        m_grid = float(omega_m) * rho_crit_0 * (box_mpc_h ** 3) / (ng_eff ** 3)

        if mp_raw is not None:
            mp_code = np.ascontiguousarray(mp_raw / m_grid, dtype=np.float64)
        elif n_sim is not None:
            mp_code = np.full(total_np, (ng_eff ** 3) / float(n_sim), dtype=np.float64)
        else:
            mp_code = np.ones(total_np, dtype=np.float64)
            m_grid = 1.0

        ang_mom_scale = m_grid * dx / v_scale

    ng_val = int(ng) if ng is not None else 256
    ng_global_val = float(ng_global) if ng_global is not None else 512.0

    in_halo_out, halo_id_out, ang_mom_out, ierr = fof_fortran.fof_fortran(
        xp_code, vp_code, mp_code, np_halo, nhalo, L_Box_code,
        link_code, int(np_halo_min), int(bound_check), int(verbose), softening_grid,
        a_val, b_val, s_15_val, n_threads_val, total_np,
        ng_val=ng_val, ng_global_val=ng_global_val
    )

    if ierr != 0:
        if ierr == 1:
            raise ValueError("Invalid input parameters or configuration in FOF/bound_test")
        elif ierr == 2:
            raise ValueError("Candidate halo span exceeds 15 Mpc limit (s_15)")
        elif ierr == 4:
            raise ValueError("Non-finite values (NaN/Inf) encountered during calculation")
        elif ierr == 10:
            raise ValueError("n_threads exceeds maximum supported slots (256)")
        else:
            raise RuntimeError(f"FOF Fortran computation failed with error code ierr={ierr}")

    ang_mom_scaled = ang_mom_out * ang_mom_scale
    return in_halo_out.astype(bool), halo_id_out, ang_mom_scaled


def fof(xp_list, vp_list=None, L_Box=512.0, linking_length=0.2, 
        np_halo_min=10, bound_check=False, verbose=0, return_ang_mom=False,
        n_threads=None, mp_list=None, softening=None, z=0.0,
        box=300.0, h0=0.67, omega_m=0.3, omega_l=0.7, n_sim=None,
        units='cube', include_hubble_flow=False,
        ng=None, ng_global=None):
    """
    Friends-of-Friends (FOF) halo finder.
    
    Parameters
    ----------
    xp_list : array_like or list of array_like
        粒子坐标。可以是单个数组 shape=(3, N)，或者是多个halo的列表
        [xp_halo1, xp_halo2, ...]，每个元素 shape=(3, Ni)。
        单位：格点（units='cube'）或 Mpc/h（units='physical'）。
        
    vp_list : array_like or list of array_like, optional
        粒子速度。格式与xp_list相同。
        仅当 bound_check=True 时需要。
        单位：格点超共形速度（units='cube'）或 km/s（units='physical'）。
        
    L_Box : float, default=512.0
        盒子边长，用于周期性边界条件。
        单位：格点（units='cube'）或 Mpc/h（units='physical'）。
        
    linking_length : float, default=0.2
        FOF链接长度。两个粒子距离小于此值时被链接。
        单位：格点（units='cube'）或 Mpc/h（units='physical'）。
        
    np_halo_min : int, default=10
        最小halo粒子数阈值。粒子数少于此值的FOF链不被识别为halo。
        
    bound_check : bool, default=False
        是否执行自洽去束缚检查。如果为True，需要提供vp_list。
        
    verbose : int, default=0
        详细输出级别（0=静默，1=阶段进度，2=详细信息）。
        
    return_ang_mom : bool, default=False
        是否返回各个输入区域累加的自转角动量。
        
    n_threads : int, optional
        OpenMP 最大线程数。默认 None 自动使用系统最大核心数。最大支持 256。
        
    mp_list : scalar, array_like, or list of array_like, optional
        粒子质量。可传入标量、单数组或列表。
        units='cube' 时为相对格点质量权重，units='physical' 时单位为 M_sun/h。
        
    softening : float, optional
        引力软化长度。units='cube' 默认 0.1 格点；units='physical' 必须显式指定（单位 Mpc/h）。
        
    z : float, default=0.0
        宇宙学红移。用于计算尺度因子 a = 1 / (1 + z)。
        
    box : float, default=300.0
        用于计算 15 Mpc 对应的格点跨度上限的物理盒子大小（Mpc/h）。
        
    h0 : float, default=0.67
        无量纲哈勃常数 h = H0 / (100 km/s/Mpc)。
        
    omega_m : float, default=0.3
        物质密度参数。
        
    omega_l : float, default=0.7
        暗能量密度参数。
        
    n_sim : float or int, optional
        全模拟总粒子数，用于 units='physical' 且未指定 mp_list 时推导等效单粒子质量。
        
    units : {'cube', 'physical'}, default='cube'
        单位模式。'cube' 为标准格点单位，'physical' 为宇宙学物理单位。
        
    include_hubble_flow : bool, default=False
        束缚检查时是否计入哈勃哈勃膨胀相对退行速度。
        
    Returns
    -------
    in_halo : ndarray of bool
        一维布尔数组，shape=(total_np,)，标记每个粒子是否属于有效束缚 halo。
        
    ang_mom : ndarray of float64, optional
        仅当 return_ang_mom=True 时返回，shape=(3, nhalo)，每个输入区域累加的角动量。
    """
    in_halo, _, ang_mom = _fof_core(
        xp_list=xp_list, vp_list=vp_list, L_Box=L_Box, linking_length=linking_length,
        np_halo_min=np_halo_min, bound_check=bound_check, verbose=verbose,
        return_ang_mom=return_ang_mom, n_threads=n_threads, mp_list=mp_list,
        softening=softening, z=z, box=box, h0=h0, omega_m=omega_m, omega_l=omega_l,
        n_sim=n_sim, units=units, include_hubble_flow=include_hubble_flow,
        ng=ng, ng_global=ng_global
    )
    if return_ang_mom:
        return in_halo, ang_mom
    return in_halo


def fof_with_info(xp_list, vp_list=None, L_Box=512.0, linking_length=0.2, 
                  np_halo_min=10, bound_check=False, verbose=0, return_ang_mom=False,
                  n_threads=None, mp_list=None, softening=None, z=0.0,
                  box=300.0, h0=0.67, omega_m=0.3, omega_l=0.7, n_sim=None,
                  units='cube', include_hubble_flow=False,
                  ng=None, ng_global=None):
    """
    FOF halo finder with detailed halo catalog output.
    
    采用 O(N) 线性分组提取每个有效 Halo 的粒子下标列表与成员粒子数。
    
    Returns
    -------
    result : dict
        Dictionary containing:
        - 'in_halo': 布尔数组，shape=(total_np,)
        - 'nhalo_found': 识别出的有效 Halo 总数
        - 'np_per_halo': 一维数组，各 Halo 的成员粒子数
        - 'halo_indices': list of arrays，各 Halo 包含的粒子下标
        - 'bound_status': shape=(nhalo_found,) 数组，各 Halo 的束缚求解状态（有效收敛均为 0）
        - 'ang_mom': (3, nhalo) 数组（仅当 return_ang_mom=True 时包含）
    """
    in_halo, halo_id, ang_mom = _fof_core(
        xp_list=xp_list, vp_list=vp_list, L_Box=L_Box, linking_length=linking_length,
        np_halo_min=np_halo_min, bound_check=bound_check, verbose=verbose,
        return_ang_mom=return_ang_mom, n_threads=n_threads, mp_list=mp_list,
        softening=softening, z=z, box=box, h0=h0, omega_m=omega_m, omega_l=omega_l,
        n_sim=n_sim, units=units, include_hubble_flow=include_hubble_flow,
        ng=ng, ng_global=ng_global
    )

    total_np = len(halo_id)
    max_id = int(np.max(halo_id)) if total_np > 0 else 0

    if max_id == 0:
        nhalo_found = 0
        np_per_halo = np.array([], dtype=int)
        halo_indices = []
    else:
        valid_mask = (halo_id > 0)
        valid_indices = np.where(valid_mask)[0]
        valid_hids = halo_id[valid_indices]

        # 统计每个有效 halo 的粒子数（halo_id 编号为 1..max_id）
        counts = np.bincount(valid_hids, minlength=max_id + 1)[1:max_id + 1]

        # 按 halo_id 稳定排序粒子下标
        order = np.argsort(valid_hids, kind='stable')
        sorted_indices = valid_indices[order]

        # 线性切分
        split_points = np.cumsum(counts)[:-1]
        raw_splits = np.split(sorted_indices, split_points) if len(split_points) > 0 else [sorted_indices]

        # 仅保留非空 halo（正常情况下 counts 全部 > 0）
        halo_indices = [idx for idx, c in zip(raw_splits, counts) if c > 0]
        np_per_halo = counts[counts > 0]
        nhalo_found = len(halo_indices)

    res = {
        'in_halo': in_halo,
        'nhalo_found': nhalo_found,
        'np_per_halo': np_per_halo,
        'halo_indices': halo_indices,
        'bound_status': np.zeros(nhalo_found, dtype=np.int32)
    }
    if return_ang_mom:
        res['ang_mom'] = ang_mom
    return res

