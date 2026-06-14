# CUBE FoF (Friends-of-Friends) Halo Finder

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

用于宇宙学模拟的 Friends-of-Friends (FOF) 光环识别算法的 Python 接口，核心算法使用 Fortran 实现。

## 功能特性

- **快速 FOF 算法**：Fortran 核心实现，支持 OpenMP 并行
- **Python 接口**：简洁易用的 Python API，支持 NumPy 数组
- **束缚检查**：可选的基于速度的引力束缚检查
- **多光环支持**：单次调用可处理多个光环
- **周期性边界条件**：自动处理周期性边界
- **可视化工具**：内置绘图函数用于结果可视化

## 安装

### 系统要求

- Python >= 3.8
- NumPy >= 1.20
- Matplotlib >= 3.0
- Fortran 编译器（推荐 gfortran）

### 从源码安装

```bash
git clone https://github.com/yourusername/py_fof.git
cd py_fof
pip install .
```

Fortran 代码会在安装过程中自动编译。

### 手动编译（可选）

如果自动编译失败，可以手动编译 Fortran 代码：

```bash
python -m numpy.f2py -c fof_fortran.f90 -m fof_fortran --f90flags='-O3 -fopenmp'
mv fof_fortran*.so py_fof/
```

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

### `fof(xp_list, vp_list=None, L_Box=512.0, linking_length=0.2, np_halo_min=10, bound_check=False, verbose=0)`

Friends-of-Friends 光环识别器。

**参数：**

- `xp_list` (array_like 或 list of array_like): 粒子坐标。可以是单个数组 shape=(3, N)，或者是多个光环的列表。单位：格点。
- `vp_list` (array_like 或 list of array_like, 可选): 粒子速度，shape=(3, N) 或列表。束缚检查时需要。
- `L_Box` (float, 可选): 盒子边长。默认：512.0
- `linking_length` (float, 可选): FOF 链接长度。默认：0.2
- `np_halo_min` (int, 可选): 形成光环的最小粒子数。默认：10
- `bound_check` (bool, 可选): 启用基于速度的束缚检查。默认：False
- `verbose` (int, 可选): 详细输出级别（0=静默，1=进度，2=详细）。默认：0

**返回：**

- `in_halo` (numpy.ndarray): 布尔数组，指示每个粒子是否属于光环。

## 算法原理

FOF 算法通过链接彼此距离在指定链接长度内的粒子来识别光环。算法包括：

1. **网格构建**：将粒子分箱到 3D 网格中，用于高效的邻居搜索
2. **FOF 搜索**：链接链接长度内的相邻粒子
3. **链合并**：使用并查算法将链接的粒子分组为光环
4. **束缚检查**（可选）：使用以下条件检查粒子的引力束缚：
   ```
   0.5 * v_rel² < N_particles / r
   ```
   其中 v_rel 是相对速度，r 是到光环中心的距离。

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

## 致谢

- 原始 Fortran 实现基于标准 FOF 算法
- 受 GADGET 和 RAMSES 等宇宙学模拟代码启发

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
