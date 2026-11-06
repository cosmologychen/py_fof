# CUBE FoF (Friends-of-Friends) Halo Finder

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

用于宇宙学模拟的 Friends-of-Friends (FOF) 光环识别算法的 Python 接口，核心算法使用 Fortran 实现。

## 功能特性

- **高效 FOF 算法**：基于空间两遍 CRS + 13 半邻域 + 无锁原子并查集 (DSU) 的 Fortran/C 核心实现，支持 OpenMP 多核并行
- **自适应引力去束缚**：内嵌多分辨率孤立 PM (48~256 网格) + C2 短程连续软化残差与解析直接求和，按比能量分批迭代自洽剔除未束缚粒子
- **双单位模式与变质量支持**：原生支持 CUBE 模拟格点单位 (`cube`) 与天体物理宇宙学单位 (`physical`：Mpc/h, km/s, \(M_\odot/h\))，支持等质量与逐粒子变质量权重
- **Python 易用 API**：支持单/多 Halo、周期性边界展开、自动 $O(N)$ 提取 Halo 成员目录与角动量
- **自动化构建**：导入时智能按源文件时间戳混合编译 C 与 Fortran 扩展，无缝集成


## 安装

### 系统要求

- Python >= 3.8
- NumPy >= 1.20
- Matplotlib >= 3.0
- C 编译器（gcc）
- Fortran 编译器（推荐 gfortran，需支持 OpenMP）
- 单精度 FFTW3 库（Ubuntu: `libfftw3-dev`，CentOS: `fftw-devel`）

### 从源码安装

```bash
git clone https://github.com/yourusername/py_fof.git
cd py_fof
pip install .
```

Fortran 和 C 代码会在安装或首次导入时自动调用混合编译生成扩展。

### 手动编译（可选）

如需手动编译，详见 [INSTALL.md](INSTALL.md)。也可直接使用：

```bash
python -c "import py_fof"
```
系统将自动检测源码时间戳并在 `.test_ai/f2py_build` 中构建扩展放入 `py_fof/`。

## 快速开始

### 基本用法

```python
import numpy as np
from py_fof import fof

# 生成测试粒子（单个光环）
np.random.seed(42)
xp = np.random.rand(3, 1000).astype(np.float64) * 512  # 1000个粒子在512^3盒子中

# 运行 FOF
in_halo = fof(xp, L_Box=512.0, linking_length=0.2, np_halo_min=10)

print(f"光环内粒子数: {np.sum(in_halo)} / {len(in_halo)}")
```

### 多个光环与束缚检查

```python
# 生成多个光环
xp1 = np.random.rand(3, 500).astype(np.float64) * 100 + 50
xp2 = np.random.rand(3, 500).astype(np.float64) * 100 + 200
xp3 = np.random.rand(3, 500).astype(np.float64) * 100 + 350

xp_list = [xp1, xp2, xp3]

# 生成速度
vp1 = np.random.randn(3, 500).astype(np.float64) * 0.5
vp2 = np.random.randn(3, 500).astype(np.float64) * 0.5
vp3 = np.random.randn(3, 500).astype(np.float64) * 0.5

vp_list = [vp1, vp2, vp3]

# 运行带束缚检查的 FOF
in_halo = fof(
    xp_list, 
    vp_list=vp_list,
    L_Box=512.0, 
    linking_length=0.2, 
    np_halo_min=10,
    bound_check=True,
    verbose=1
)

print(f"光环内总粒子数: {np.sum(in_halo)}")
```

### 可视化

```python
from py_fof.test_fof import plot_particles

# 绘制结果
plot_particles(
    np.hstack(xp_list), 
    in_halo, 
    title="FOF 光环识别结果",
    save_path="halo_result.png"
)
```

## API 参考

### `fof(...)`
```python
def fof(xp_list, vp_list=None, L_Box=512.0, linking_length=0.2,
        np_halo_min=10, bound_check=False, verbose=0, return_ang_mom=False,
        n_threads=None, mp_list=None, softening=None, z=0.0,
        box=300.0, h0=0.67, omega_m=0.3, omega_l=0.7, n_sim=None,
        units='cube', include_hubble_flow=False):
```

**核心参数说明：**

- `xp_list` (array_like 或 list of array_like): 粒子坐标。单个数组 shape=(3, N) 或多个 Halo 的列表 [xp1, xp2, ...]，每个元素 shape=(3, Ni)。单位由 `units` 决定。
- `vp_list` (array_like 或 list of array_like, 可选): 粒子速度。shape=(3, N) 或列表。当 `bound_check=True` 时必填。
- `L_Box` (float, 可选): 模拟盒边长，用于周期性边界条件。默认：512.0。
- `linking_length` (float, 可选): FOF 链接长度。默认：0.2。
- `np_halo_min` (int, 可选): 构成有效 Halo 的最小成员粒子数。少于此数的候选将被过滤。默认：10。
- `bound_check` (bool, 可选): 是否开启自洽引力去束缚检查。默认：False。
- `verbose` (int, 可选): 详细输出级别（0=静默，1=阶段进度，2=详细信息）。默认：0。
- `return_ang_mom` (bool, 可选): 是否返回输入区域累加的自转角动量。默认：False。
- `n_threads` (int, 可选): OpenMP 最大线程数。默认 None 自动利用系统全部核心，最大支持 256。
- `mp_list` (标量、数组或数组列表, 可选): 粒子质量。`units='cube'` 下为无量纲格点质量权重（默认 1.0）；`units='physical'` 下单位为 \(M_\odot/h\)。
- `softening` (float, 可选): 引力软化长度。`units='cube'` 默认为 0.1 格点；`units='physical'` 下必须显式传入（共动 Mpc/h）。
- `z` (float, 可选): 宇宙学红移，默认 0.0。
- `units` (`'cube'` 或 `'physical'`, 可选): 单位系统。默认 `'cube'`。当为 `'physical'` 时，坐标/链接长度单位为 Mpc/h，速度单位为 km/s，质量单位为 \(M_\odot/h\)。
- `include_hubble_flow` (bool, 可选): 去束缚计算时是否将相对哈勃膨胀退行速度计入动能。默认：False。

**返回：**
- `in_halo` (numpy.ndarray): 一维布尔数组，指示每个粒子是否属于最终有效束缚 Halo。
- `ang_mom` (numpy.ndarray, 可选): 当 `return_ang_mom=True` 时返回，shape=(3, nhalo)。

### `fof_with_info(...)`
参数与 `fof(...)` 完全相同，额外通过 $O(N)$ 线性提取并返回详细字典：
- `'in_halo'`: 成员布尔标记数组；
- `'nhalo_found'`: 识别出的有效 Halo 总数；
- `'np_per_halo'`: 各 Halo 成员粒子数一维数组；
- `'halo_indices'`: list of ndarray，每个 Halo 包含的原输入粒子全局下标；
- `'ang_mom'`: 各区域角动量（仅当 `return_ang_mom=True` 时）。

## 算法原理

1. **两遍空间 CRS 构建与周期展开**：
   通过第一遍统计与第二遍填充，快速构建三维网格的压缩行存储（CRS）粒子链表，消除动态指针并保证内存局部性。
2. **13 半邻域与并行原子并查集 (DSU)**：
   采用 13 半邻域（Half-neighborhood）对称搜索避免重复距离计算；在多核并行下利用 C 语言无锁原子并查集 (`atomic_dsu.c`) 执行高效连通图合并。
3. **自适应自洽引力去束缚 (`bound_test`)**：
   - **直接求和**（\(N \le 1000\)）：基于连续软化解析势直接求和，迭代中增量扣除已剔除粒子的贡献；
   - **多分辨率孤立 PM + C2 短程连续软化残差**（\(N > 1000\)）：采用 48/64/128/256 自适应孤立网格、最优 Green 函数与自势扣除，结合 cell-CSR 连续软化短程相互作用；
   - **2% 比能量分批收敛迭代**：按比束缚能每次剔除 2% 束缚最弱的粒子，最多迭代 20 轮外层势场刷新，严格收敛判定。

## 测试

运行测试套件：

```bash
python test_fof.py
```

这将生成多个测试图像：
- `test_predictable_halo.png`: 单光环测试
- `test_predictable_halo_velocity.png`: 带速度束缚检查的光环
- `test_multiple_predictable_halos.png`: 多光环测试
- `test_multiple_predictable_halos_velocity.png`: 带速度的多光环

## 性能

Fortran 实现配合 OpenMP 并行化提供了出色的性能：
- ~10^5 粒子：< 1 秒
- ~10^6 粒子：~10 秒
- ~10^7 粒子：~2 分钟

（性能取决于硬件和链接长度）

### 性能测试结果

测试配置：10个光环，每个10,000个粒子，总共100,000个粒子

| 测试项              | 时间       | 速度                    |
| ------------------- | ---------- | ----------------------- |
| FOF (无bound_check) | **0.44秒** | **225,117 particles/s** |
| FOF + bound_check   | **0.43秒** | **230,554 particles/s** |

## 贡献

欢迎贡献！请随时提交 Pull Request。

## 许可证

本项目采用 MIT 许可证 - 详情见 LICENSE 文件。

## 引用

如果您在研究中使用此代码，请引用：

```bibtex
@software{py_fof,
  author = {Your Name},
  title = {CUBE FoF: Friends-of-Friends Halo Finder},
  url = {https://github.com/yourusername/py_fof},
  year = {2024}
}
```

## 常见问题

### Q: 为什么有些粒子被排除在光环之外？

A: 当启用 `bound_check=True` 时，算法会检查粒子的引力束缚性。高速粒子可能无法被光环引力束缚，因此会被排除。这是物理上合理的行为。

### Q: 如何选择合适的链接长度？

A: 链接长度通常设置为平均粒子间距的 0.2 倍左右。对于周期性盒子，可以使用：
```python
linking_length = 0.2 * L_Box / N_particles**(1/3)
```

### Q: 支持哪些数据格式？

A: 输入应为 NumPy 数组，shape 为 (3, N) 或列表形式。支持 float32 和 float64 类型。

## 联系方式

如有问题或建议，请提交 Issue 或联系：your.email@example.com
