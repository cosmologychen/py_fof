"""
Test script for py_fof package

Tests the FOF halo finder with random particle distributions.
"""

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# =============================================================================
# 全局测试参数 - 可根据需要调整
# =============================================================================

N_PARTICLES = 10000           # 每个halo的粒子数
BOX_SIZE = 40.0            # 粒子分布范围（0到BOX_SIZE）
L_BOX = 512.0               # 周期性盒子边长
LINKING_LENGTH = 2.0       # FOF链接长度
NP_HALO_MIN = 10            # 最小halo粒子数阈值

VELOCITY_SCALE = 0.1        # 速度分布的标准差（适配 a/(6*pi) 势场深度）

# =============================================================================

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

OS_TEST_DIR = ".test_ai"
if not os.path.exists(OS_TEST_DIR):
    os.makedirs(OS_TEST_DIR)

try:
    from py_fof import fof, fof_with_info
    print("Successfully imported py_fof module")
except ImportError as e:
    print(f"Failed to import py_fof: {e}")
    print("\nPlease build the Fortran extension first:")
    print("  cd /home/cbh/py_fof")
    print("  python setup.py build_ext --inplace")
    sys.exit(1)


def plot_particles(xp, in_halo, unbound_in_halo=None, title="FOF Halo Finder Result", save_path=None):
    """
    绘制粒子分布图，区分halo内外的粒子
    
    参数:
        xp: 粒子坐标数组，形状为 (3, N)
        in_halo: 布尔数组，标记粒子是否在halo中
        unbound_in_halo: 布尔数组，标记在halo中但不束缚的粒子（绿色显示）
        title: 图标题
        save_path: 保存路径，如果为None则显示图像
    """
    fig = plt.figure(figsize=(14, 12))
    
    # 2x2布局：3D图 + XY/XZ/YZ三个投影
    
    # 左上：3D图
    ax1 = fig.add_subplot(221, projection='3d')
    
    not_in_halo = ~in_halo
    if np.any(not_in_halo):
        ax1.scatter(xp[0, not_in_halo], xp[1, not_in_halo], xp[2, not_in_halo], 
                   c='blue', s=20, alpha=0.5, edgecolors='darkblue', linewidths=0.5,
                   label=f'Not in halo ({np.sum(not_in_halo)})')
    
    # 在halo中但不束缚的粒子（绿色）
    if unbound_in_halo is not None and np.any(unbound_in_halo):
        ax1.scatter(xp[0, unbound_in_halo], xp[1, unbound_in_halo], xp[2, unbound_in_halo], 
                   c='green', s=25, alpha=0.7, edgecolors='darkgreen', linewidths=0.5,
                   label=f'Unbound in halo ({np.sum(unbound_in_halo)})')
    
    # 在halo中且束缚的粒子（红色）
    bound_in_halo = in_halo.copy()
    if unbound_in_halo is not None:
        bound_in_halo = in_halo & ~unbound_in_halo
    if np.any(bound_in_halo):
        ax1.scatter(xp[0, bound_in_halo], xp[1, bound_in_halo], xp[2, bound_in_halo], 
                   c='red', s=25, alpha=0.7, edgecolors='darkred', linewidths=0.5,
                   label=f'Bound in halo ({np.sum(bound_in_halo)})')
    
    ax1.set_xlabel('X', fontsize=10)
    ax1.set_ylabel('Y', fontsize=10)
    ax1.set_zlabel('Z', fontsize=10)
    ax1.set_title(f'{title}\n(3D View)', fontsize=11)
    ax1.legend(fontsize=9)
    
    # 右上：XY投影
    ax2 = fig.add_subplot(222)
    
    if np.any(not_in_halo):
        ax2.scatter(xp[0, not_in_halo], xp[1, not_in_halo], 
                   c='blue', s=15, alpha=0.5, edgecolors='darkblue', linewidths=0.5,
                   label=f'Not in halo ({np.sum(not_in_halo)})')
    
    if unbound_in_halo is not None and np.any(unbound_in_halo):
        ax2.scatter(xp[0, unbound_in_halo], xp[1, unbound_in_halo], 
                   c='green', s=25, alpha=0.7, edgecolors='darkgreen', linewidths=0.5,
                   label=f'Unbound in halo ({np.sum(unbound_in_halo)})')
    
    if np.any(bound_in_halo):
        ax2.scatter(xp[0, bound_in_halo], xp[1, bound_in_halo], 
                   c='red', s=25, alpha=0.7, edgecolors='darkred', linewidths=0.5,
                   label=f'Bound in halo ({np.sum(bound_in_halo)})')
    
    ax2.set_xlabel('X', fontsize=10)
    ax2.set_ylabel('Y', fontsize=10)
    ax2.set_title('X-Y Projection', fontsize=11)
    ax2.legend(fontsize=9)
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)
    
    # 左下：XZ投影
    ax3 = fig.add_subplot(223)
    
    if np.any(not_in_halo):
        ax3.scatter(xp[0, not_in_halo], xp[2, not_in_halo], 
                   c='blue', s=15, alpha=0.5, edgecolors='darkblue', linewidths=0.5,
                   label=f'Not in halo ({np.sum(not_in_halo)})')
    
    if unbound_in_halo is not None and np.any(unbound_in_halo):
        ax3.scatter(xp[0, unbound_in_halo], xp[2, unbound_in_halo], 
                   c='green', s=25, alpha=0.7, edgecolors='darkgreen', linewidths=0.5,
                   label=f'Unbound in halo ({np.sum(unbound_in_halo)})')
    
    if np.any(bound_in_halo):
        ax3.scatter(xp[0, bound_in_halo], xp[2, bound_in_halo], 
                   c='red', s=25, alpha=0.7, edgecolors='darkred', linewidths=0.5,
                   label=f'Bound in halo ({np.sum(bound_in_halo)})')
    
    ax3.set_xlabel('X', fontsize=10)
    ax3.set_ylabel('Z', fontsize=10)
    ax3.set_title('X-Z Projection', fontsize=11)
    ax3.legend(fontsize=9)
    ax3.set_aspect('equal')
    ax3.grid(True, alpha=0.3)
    
    # 右下：YZ投影
    ax4 = fig.add_subplot(224)
    
    if np.any(not_in_halo):
        ax4.scatter(xp[1, not_in_halo], xp[2, not_in_halo], 
                   c='blue', s=15, alpha=0.5, edgecolors='darkblue', linewidths=0.5,
                   label=f'Not in halo ({np.sum(not_in_halo)})')
    
    if unbound_in_halo is not None and np.any(unbound_in_halo):
        ax4.scatter(xp[1, unbound_in_halo], xp[2, unbound_in_halo], 
                   c='green', s=25, alpha=0.7, edgecolors='darkgreen', linewidths=0.5,
                   label=f'Unbound in halo ({np.sum(unbound_in_halo)})')
    
    if np.any(bound_in_halo):
        ax4.scatter(xp[1, bound_in_halo], xp[2, bound_in_halo], 
                   c='red', s=25, alpha=0.7, edgecolors='darkred', linewidths=0.5,
                   label=f'Bound in halo ({np.sum(bound_in_halo)})')
    
    ax4.set_xlabel('Y', fontsize=10)
    ax4.set_ylabel('Z', fontsize=10)
    ax4.set_title('Y-Z Projection', fontsize=11)
    ax4.legend(fontsize=9)
    ax4.set_aspect('equal')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"  Figure saved to: {save_path}")
    else:
        plt.show()
    
    plt.close()


# %%
def test_performance():
    """
    性能测试：测试大规模FOF算法的速度
    
    测试配置：
    - 10个halo
    - 每个halo 10000个粒子
    - 总共100000个粒子
    """
    import time
    
    print("="*60)
    print("Test: Performance Test (Large Scale)")
    print("="*60)
    
    n_halos = 10
    n_particles_per_halo = 10000
    
    print(f"\nConfiguration:")
    print(f"  Number of halos: {n_halos}")
    print(f"  Particles per halo: {n_particles_per_halo}")
    print(f"  Total particles: {n_halos * n_particles_per_halo}")
    print(f"  Linking length: {LINKING_LENGTH}")
    print(f"  L_Box: {L_BOX}")
    
    # 生成halo数据
    print("\nGenerating halo data...")
    start_time = time.time()
    
    xp_list = []
    vp_list = []
    
    # 定义每个halo的偏移量，确保它们在空间上分开
    halo_shifts = []
    for i in range(n_halos):
        # 在盒子中均匀分布
        shift_x = (i % 3) * 150 + 50
        shift_y = (i // 3) * 150 + 50
        shift_z = (i // 9) * 150 + 50
        halo_shifts.append(np.array([shift_x, shift_y, shift_z]))
    
    for ihalo in range(n_halos):
        # 生成随机粒子分布（跨度 12.0 格点 < 17.152，满足 15 Mpc 物理上限并形成致密大 Halo）
        xp = np.random.rand(3, n_particles_per_halo).astype(np.float64) * 12.0 + halo_shifts[ihalo].reshape(3, 1)
        
        # 生成速度
        vp = np.random.randn(3, n_particles_per_halo).astype(np.float64) * VELOCITY_SCALE
        
        xp_list.append(xp)
        vp_list.append(vp)
    
    gen_time = time.time() - start_time
    print(f"  Data generation time: {gen_time:.2f} seconds")
    
    # 测试FOF（不带bound_check）
    print("\nRunning FOF (without bound_check)...")
    start_time = time.time()
    in_halo = fof(xp_list, L_Box=L_BOX, linking_length=LINKING_LENGTH, np_halo_min=NP_HALO_MIN)
    fof_time = time.time() - start_time
    
    print(f"  FOF time: {fof_time:.2f} seconds")
    print(f"  Particles in halos: {np.sum(in_halo)} / {n_halos * n_particles_per_halo}")
    print(f"  Speed: {n_halos * n_particles_per_halo / fof_time:.0f} particles/second")
    
    # 测试FOF（带bound_check）
    print("\nRunning FOF (with bound_check)...")
    start_time = time.time()
    in_halo_bound = fof(xp_list, vp_list=vp_list, L_Box=L_BOX, linking_length=LINKING_LENGTH, 
                        np_halo_min=NP_HALO_MIN, bound_check=True, verbose=0)
    bound_time = time.time() - start_time
    
    print(f"  FOF + bound_check time: {bound_time:.2f} seconds")
    print(f"  Particles in halos (bound): {np.sum(in_halo_bound)} / {n_halos * n_particles_per_halo}")
    print(f"  Speed: {n_halos * n_particles_per_halo / bound_time:.0f} particles/second")
    
    print("\n" + "="*60)
    print("Performance Summary:")
    print(f"  Total particles: {n_halos * n_particles_per_halo}")
    print(f"  FOF time: {fof_time:.2f} s ({n_halos * n_particles_per_halo / fof_time:.0f} particles/s)")
    print(f"  FOF + bound_check time: {bound_time:.2f} s ({n_halos * n_particles_per_halo / bound_time:.0f} particles/s)")
    print("="*60)
    
    assert np.sum(in_halo) > 0, "Expected non-zero particles in halos"
    assert np.sum(in_halo_bound) > 0, "Expected non-zero particles in bound halos"
    
    return xp_list, in_halo, in_halo_bound


# %%
# def generate_random_halo(n_particles, box_size, offset=0.0, seed=None):
#     """
#     生成完全随机的halo粒子分布
#     
#     参数:
#         n_particles: 粒子数量
#         box_size: 分布范围大小
#         offset: 坐标偏移量
#         seed: 随机种子
#     
#     返回:
#         xp: 粒子坐标数组 (3, n_particles)
#     """
#     if seed is not None:
#         np.random.seed(seed)
#     
#     xp = np.random.rand(3, n_particles).astype(np.float64) * box_size + offset
#     return xp


# %%
# def generate_multiple_halos_with_velocity(n_halos, n_halo_particles_per_halo, linking_length,
#                                           velocity_scale=0.5, high_velocity_fraction=0.1,
#                                           seed=None):
#     """
#     生成多个可预测的halo，用于验证FOF算法
#     
#     参数:
#         n_halos: halo数量
#         n_halo_particles_per_halo: 每个halo的粒子数量（可以是整数或列表）
#         linking_length: FOF链接长度
#         velocity_scale: 速度的标准差
#         high_velocity_fraction: 高速粒子的比例（这些粒子会被bound_check排除）
#         seed: 随机种子
#     
#     返回:
#         xp_list: halo坐标列表 [xp1, xp2, ...]
#         vp_list: 速度列表 [vp1, vp2, ...]
#         expected_in_halo_list: 每个halo期望的成员掩码列表
#     """
#     if seed is not None:
#         np.random.seed(seed)
#     
#     # 处理n_halo_particles_per_halo参数
#     if isinstance(n_halo_particles_per_halo, int):
#         n_halo_particles_per_halo = [n_halo_particles_per_halo] * n_halos
#     
#     xp_list = []
#     vp_list = []
#     expected_in_halo_list = []
#     
#     # 生成多个halo
#     for ihalo in range(n_halos):
#         n_halo = n_halo_particles_per_halo[ihalo]
#         
#         # 生成halo粒子
#         xp_halo = np.zeros((3, n_halo), dtype=np.float64)
#         
#         # 随机偏移，让不同halo在不同位置
#         offset = np.random.uniform(0, 100, 3)
#         xp_halo[:, 0] = offset
#         
#         for i in range(1, n_halo):
#             r = np.random.uniform(linking_length * 0.5, linking_length * 0.95)
#             theta = np.random.uniform(0, 2 * np.pi)
#             phi = np.arccos(np.random.uniform(-1, 1))
#             
#             dx = r * np.sin(phi) * np.cos(theta)
#             dy = r * np.sin(phi) * np.sin(theta)
#             dz = r * np.cos(phi)
#             
#             xp_halo[0, i] = xp_halo[0, i-1] + dx
#             xp_halo[1, i] = xp_halo[1, i-1] + dy
#             xp_halo[2, i] = xp_halo[2, i-1] + dz
#         
#         # 生成速度
#         vp_halo = np.random.randn(3, n_halo).astype(np.float64) * velocity_scale
#         
#         # 给一部分粒子极高的速度（这些粒子会被bound_check排除）
#         n_high_v = max(1, int(n_halo * high_velocity_fraction))
#         high_v_indices = np.random.choice(n_halo, n_high_v, replace=False)
#         vp_halo[:, high_v_indices] *= 10.0  # 10倍速度
#         
#         # 期望的halo成员（排除高速粒子）
#         expected_in_halo = np.ones(n_halo, dtype=bool)
#         expected_in_halo[high_v_indices] = False
#         
#         xp_list.append(xp_halo)
#         vp_list.append(vp_halo)
#         expected_in_halo_list.append(expected_in_halo)
#     
#     return xp_list, vp_list, expected_in_halo_list


# %%
# 生成可预测的单个halo
def generate_predictable_halo(n_halo_particles, linking_length, n_noise_particles, seed=None):
    """
    生成可预测的halo，用于验证FOF算法
    
    生成策略:
        1. 第一个粒子在原点(0,0,0)
        2. 后续粒子在前一个粒子的linking_length距离内
        3. 这样所有halo粒子一定被FOF识别为一个halo
        4. 在halo范围外5倍linking_length距离处生成噪声粒子
        5. 噪声粒子之间距离大于linking_length
    
    参数:
        n_halo_particles: halo内粒子数量
        linking_length: FOF链接长度
        n_noise_particles: 噪声粒子数量
        seed: 随机种子
    
    返回:
        xp: 粒子坐标数组 (3, n_halo_particles + n_noise_particles)
        expected_in_halo: 期望的halo成员掩码（前n_halo_particles个为True）
        halo_center: halo的质心位置
        halo_range: halo的坐标范围
    """
    if seed is not None:
        np.random.seed(seed)
    
    # 生成halo粒子：每个粒子在前一个粒子的linking_length距离内
    # 这样形成一个链状结构，保证所有粒子都能被FOF连接
    xp_halo = np.zeros((3, n_halo_particles), dtype=np.float64)
    
    for i in range(1, n_halo_particles):
        # 在前一个粒子的linking_length距离内随机生成新粒子
        # 使用球坐标生成均匀分布的点
        # 使用较大的距离系数（0.95）让halo扩散更大
        r = np.random.uniform(linking_length * 0.5, linking_length * 0.95)
        theta = np.random.uniform(0, 2 * np.pi)
        phi = np.arccos(np.random.uniform(-1, 1))
        
        dx = r * np.sin(phi) * np.cos(theta)
        dy = r * np.sin(phi) * np.sin(theta)
        dz = r * np.cos(phi)
        
        xp_halo[0, i] = xp_halo[0, i-1] + dx
        xp_halo[1, i] = xp_halo[1, i-1] + dy
        xp_halo[2, i] = xp_halo[2, i-1] + dz
    
    # 计算halo的范围
    halo_min = xp_halo.min(axis=1)
    halo_max = xp_halo.max(axis=1)
    halo_center = (halo_min + halo_max) / 2
    halo_range = halo_max - halo_min
    
    # 生成噪声粒子：在halo范围外5倍linking_length距离处
    xp_noise = np.zeros((3, n_noise_particles), dtype=np.float64)
    
    # 计算噪声粒子的最小距离（halo边界 + 5倍linking_length）
    noise_min_dist = halo_range.max() + 5 * linking_length
    
    for i in range(n_noise_particles):
        # 在halo中心周围生成噪声粒子，距离至少为noise_min_dist
        while True:
            # 随机方向
            theta = np.random.uniform(0, 2 * np.pi)
            phi = np.arccos(np.random.uniform(-1, 1))
            r = np.random.uniform(noise_min_dist, noise_min_dist * 2)
            
            x = halo_center[0] + r * np.sin(phi) * np.cos(theta)
            y = halo_center[1] + r * np.sin(phi) * np.sin(theta)
            z = halo_center[2] + r * np.cos(phi)
            
            # 检查与已有噪声粒子的距离
            if i == 0:
                xp_noise[0, i] = x
                xp_noise[1, i] = y
                xp_noise[2, i] = z
                break
            else:
                distances = np.sqrt(np.sum((xp_noise[:, :i] - np.array([x, y, z])[:, np.newaxis])**2, axis=0))
                if np.all(distances > linking_length):
                    xp_noise[0, i] = x
                    xp_noise[1, i] = y
                    xp_noise[2, i] = z
                    break
    
    # 合并halo粒子和噪声粒子
    xp = np.hstack([xp_halo, xp_noise])
    
    # 期望的halo成员掩码
    expected_in_halo = np.zeros(n_halo_particles + n_noise_particles, dtype=bool)
    expected_in_halo[:n_halo_particles] = True
    
    # 打乱粒子顺序
    n_total = n_halo_particles + n_noise_particles
    shuffle_indices = np.random.permutation(n_total)
    xp = xp[:, shuffle_indices]
    expected_in_halo = expected_in_halo[shuffle_indices]
    
    return xp, expected_in_halo, halo_center, halo_range


# %%
# def test_single_halo():
#     """测试单个halo的情况"""
#     print("\n" + "="*60)
#     print("Test 1: Single halo")
#     print("="*60)
#     
#     np.random.seed(42)
#     
#     xp = np.random.rand(3, N_PARTICLES).astype(np.float64) * BOX_SIZE
#     
#     print(f"  N_particles: {N_PARTICLES}")
#     print(f"  Box size: {BOX_SIZE}")
#     print(f"  Linking length: {LINKING_LENGTH}")
#     print(f"  Particle coordinates shape: {xp.shape}")
#     print(f"  Coordinate range: [{xp.min():.2f}, {xp.max():.2f}]")
#     
#     in_halo = fof(xp, L_Box=L_BOX, linking_length=LINKING_LENGTH, np_halo_min=NP_HALO_MIN)
#     
#     n_in_halo = np.sum(in_halo)
#     print(f"\n  Particles in halos: {n_in_halo} / {N_PARTICLES}")
#     print(f"  Fraction: {n_in_halo / N_PARTICLES * 100:.1f}%")
#     
#     plot_particles(xp, in_halo, title=f"Single Halo Test (N={N_PARTICLES}, L={LINKING_LENGTH})",
#                    save_path="test_single_halo.png")
#     
#     return xp, in_halo


# %%
# def test_multiple_halos():
#     """测试多个halo的情况"""
#     print("\n" + "="*60)
#     print("Test 2: Multiple halos")
#     print("="*60)
#     
#     np.random.seed(123)
#     
#     # 三个halo，分别放在不同位置
#     offset1 = BOX_SIZE * 0.5
#     offset2 = BOX_SIZE * 2.0
#     offset3 = BOX_SIZE * 3.5
#     
#     xp1 = np.random.rand(3, N_PARTICLES).astype(np.float64) * BOX_SIZE + offset1
#     xp2 = np.random.rand(3, N_PARTICLES).astype(np.float64) * BOX_SIZE + offset2
#     xp3 = np.random.rand(3, N_PARTICLES).astype(np.float64) * BOX_SIZE + offset3
#     
#     xp_list = [xp1, xp2, xp3]
#     
#     print(f"  N_particles per halo: {N_PARTICLES}")
#     print(f"  Box size: {BOX_SIZE}")
#     print(f"  Linking length: {LINKING_LENGTH}")
#     print(f"  Number of halos: {len(xp_list)}")
#     for i, xp in enumerate(xp_list):
#         print(f"  Halo {i+1}: {xp.shape[1]} particles, range [{xp.min():.1f}, {xp.max():.1f}]")
#     
#     in_halo = fof(xp_list, L_Box=L_BOX, linking_length=LINKING_LENGTH, np_halo_min=NP_HALO_MIN)
#     
#     total_particles = sum(xp.shape[1] for xp in xp_list)
#     n_in_halo = np.sum(in_halo)
#     print(f"\n  Total particles: {total_particles}")
#     print(f"  Particles in halos: {n_in_halo}")
#     print(f"  Fraction: {n_in_halo / total_particles * 100:.1f}%")
#     
#     xp_all = np.hstack(xp_list)
#     plot_particles(xp_all, in_halo, title=f"Multiple Halos Test (3 halos, N={N_PARTICLES} each)",
#                    save_path="test_multiple_halos.png")
#     
#     return xp_all, in_halo


# %%
# def test_with_binding_check():
#     """测试带束缚检查的情况"""
#     print("\n" + "="*60)
#     print("Test 3: Multiple halos with binding check")
#     print("="*60)
#     
#     np.random.seed(456)
#     
#     offset1 = BOX_SIZE * 0.5
#     offset2 = BOX_SIZE * 2.0
#     offset3 = BOX_SIZE * 3.5
#     
#     xp1 = np.random.rand(3, N_PARTICLES).astype(np.float64) * BOX_SIZE + offset1
#     xp2 = np.random.rand(3, N_PARTICLES).astype(np.float64) * BOX_SIZE + offset2
#     xp3 = np.random.rand(3, N_PARTICLES).astype(np.float64) * BOX_SIZE + offset3
#     
#     vp1 = np.random.randn(3, N_PARTICLES).astype(np.float64) * VELOCITY_SCALE
#     vp2 = np.random.randn(3, N_PARTICLES).astype(np.float64) * VELOCITY_SCALE
#     vp3 = np.random.randn(3, N_PARTICLES).astype(np.float64) * VELOCITY_SCALE
#     
#     xp_list = [xp1, xp2, xp3]
#     vp_list = [vp1, vp2, vp3]
#     
#     print(f"  N_particles per halo: {N_PARTICLES}")
#     print(f"  Box size: {BOX_SIZE}")
#     print(f"  Linking length: {LINKING_LENGTH}")
#     print(f"  Number of halos: {len(xp_list)}")
#     print(f"  Binding check: enabled")
#     
#     in_halo = fof(xp_list, vp_list=vp_list, L_Box=L_BOX, linking_length=LINKING_LENGTH, 
#                   np_halo_min=NP_HALO_MIN, bound_check=True)
#     
#     total_particles = sum(xp.shape[1] for xp in xp_list)
#     n_in_halo = np.sum(in_halo)
#     print(f"\n  Total particles: {total_particles}")
#     print(f"  Particles in halos (after binding check): {n_in_halo}")
#     print(f"  Fraction: {n_in_halo / total_particles * 100:.1f}%")
#     
#     xp_all = np.hstack(xp_list)
#     plot_particles(xp_all, in_halo, title=f"With Binding Check (v_scale={VELOCITY_SCALE})",
#                    save_path="test_binding_check.png")
#     
#     return xp_all, in_halo


# %%
# def test_detailed_output():
#     """测试详细输出"""
#     print("\n" + "="*60)
#     print("Test 4: Detailed output with halo information")
#     print("="*60)
#     
#     np.random.seed(789)
#     
#     offset1 = BOX_SIZE * 0.5
#     offset2 = BOX_SIZE * 2.0
#     offset3 = BOX_SIZE * 3.5
#     
#     xp1 = np.random.rand(3, N_PARTICLES).astype(np.float64) * BOX_SIZE + offset1
#     xp2 = np.random.rand(3, N_PARTICLES).astype(np.float64) * BOX_SIZE + offset2
#     xp3 = np.random.rand(3, N_PARTICLES).astype(np.float64) * BOX_SIZE + offset3
#     
#     xp_list = [xp1, xp2, xp3]
#     
#     print(f"  N_particles per halo: {N_PARTICLES}")
#     print(f"  Box size: {BOX_SIZE}")
#     print(f"  Linking length: {LINKING_LENGTH}")
#     
#     result = fof_with_info(xp_list, L_Box=L_BOX, linking_length=LINKING_LENGTH, np_halo_min=NP_HALO_MIN)
#     
#     print(f"\n  Number of halos found: {result['nhalo_found']}")
#     print(f"  Particles per halo: {result['np_per_halo']}")
#     print(f"  Total particles in halos: {np.sum(result['np_per_halo'])}")
#     
#     xp_all = np.hstack(xp_list)
#     plot_particles(xp_all, result['in_halo'], 
#                    title=f"Detailed Output (found {result['nhalo_found']} halos)",
#                    save_path="test_detailed_output.png")
#     
#     return result


# %%
# def test_detailed_output():
#     """测试详细输出"""
#     print("\n" + "="*60)
#     print("Test 4: Detailed output with halo information")
#     print("="*60)
#     
#     np.random.seed(789)
#     
#     offset1 = BOX_SIZE * 0.5
#     offset2 = BOX_SIZE * 2.0
#     offset3 = BOX_SIZE * 3.5
#     
#     xp1 = np.random.rand(3, N_PARTICLES).astype(np.float64) * BOX_SIZE + offset1
#     xp2 = np.random.rand(3, N_PARTICLES).astype(np.float64) * BOX_SIZE + offset2
#     xp3 = np.random.rand(3, N_PARTICLES).astype(np.float64) * BOX_SIZE + offset3
#     
#     xp_list = [xp1, xp2, xp3]
#     
#     print(f"  N_particles per halo: {N_PARTICLES}")
#     print(f"  Box size: {BOX_SIZE}")
#     print(f"  Linking length: {LINKING_LENGTH}")
#     
#     result = fof_with_info(xp_list, L_Box=L_BOX, linking_length=LINKING_LENGTH, np_halo_min=NP_HALO_MIN)
#     
#     print(f"\n  Number of halos found: {result['nhalo_found']}")
#     print(f"  Particles per halo: {result['np_per_halo']}")
#     print(f"  Total particles in halos: {np.sum(result['np_per_halo'])}")
#     
#     xp_all = np.hstack(xp_list)
#     plot_particles(xp_all, result['in_halo'], title=f"Detailed Output ({result['nhalo_found']} halos found)",
#                    save_path="test_detailed_output.png")
#     
#     return result


# %%
# 测试可预测的halo生成器
def test_predictable_halo():
    """测试可预测的halo生成器，验证FOF算法正确性"""
    print("\n" + "="*60)
    print("Test 5: Predictable halo validation")
    print("="*60)
    
    np.random.seed(999)
    
    # 生成可预测的halo
    n_halo = 200
    n_noise = 400
    
    xp, expected_in_halo, halo_center, halo_range = generate_predictable_halo(
        n_halo_particles=n_halo,
        linking_length=LINKING_LENGTH,
        n_noise_particles=n_noise,
        seed=999
    )
    
    print(f"  Generated halo particles: {n_halo}")
    print(f"  Generated noise particles: {n_noise}")
    print(f"  Halo center: [{halo_center[0]:.2f}, {halo_center[1]:.2f}, {halo_center[2]:.2f}]")
    print(f"  Halo range: [{halo_range[0]:.2f}, {halo_range[1]:.2f}, {halo_range[2]:.2f}]")
    print(f"  Linking length: {LINKING_LENGTH}")
    
    # 运行FOF
    in_halo, ang_mom_zero = fof(xp, L_Box=L_BOX, linking_length=LINKING_LENGTH, np_halo_min=NP_HALO_MIN, return_ang_mom=True)
    assert np.allclose(ang_mom_zero, 0.0), f"Expected zero angular momentum, got {ang_mom_zero}"
    
    # 对比结果
    expected_count = np.sum(expected_in_halo)
    found_count = np.sum(in_halo)
    
    # 计算匹配度
    true_positive = np.sum(in_halo & expected_in_halo)
    false_positive = np.sum(in_halo & ~expected_in_halo)
    false_negative = np.sum(~in_halo & expected_in_halo)
    
    print(f"\n  Expected particles in halo: {expected_count}")
    print(f"  Found particles in halo: {found_count}")
    print(f"\n  Confusion Matrix:")
    print(f"    True Positive:  {true_positive} (correctly identified)")
    print(f"    False Positive: {false_positive} (incorrectly included)")
    print(f"    False Negative: {false_negative} (incorrectly excluded)")
    
    if expected_count > 0:
        precision = true_positive / found_count if found_count > 0 else 0
        recall = true_positive / expected_count
        f1_score = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
        
        print(f"\n  Metrics:")
        print(f"    Precision: {precision:.3f}")
        print(f"    Recall:    {recall:.3f}")
        print(f"    F1 Score:  {f1_score:.3f}")
    
    # 显示每个粒子的判断结果
    correct_mask = in_halo == expected_in_halo
    n_correct = np.sum(correct_mask)
    n_total = len(correct_mask)
    
    print(f"\n  Per-particle accuracy: {n_correct}/{n_total} ({100*n_correct/n_total:.1f}%)")
    assert n_correct == n_total, f"Expected 100% accuracy, got {n_correct}/{n_total}"
    
    # 显示错误判断的粒子
    if false_positive > 0:
        fp_indices = np.where(in_halo & ~expected_in_halo)[0]
        print(f"\n  False Positive particle indices: {fp_indices[:10]}{'...' if len(fp_indices) > 10 else ''}")
    
    if false_negative > 0:
        fn_indices = np.where(~in_halo & expected_in_halo)[0]
        print(f"  False Negative particle indices: {fn_indices[:10]}{'...' if len(fn_indices) > 10 else ''}")
    
    # 绘图
    plot_particles(xp, in_halo, title=f"Predictable Halo Test (expected {n_halo}, found {found_count})",
                   save_path=os.path.join(OS_TEST_DIR, "test_predictable_halo.png"))
    
    return xp, in_halo, expected_in_halo


# %%
# 测试可预测的halo带速度（部分粒子高速）
def test_predictable_halo_with_velocity():
    """测试可预测的halo带速度，部分粒子有高速会被bound_check排除"""
    print("\n" + "="*60)
    print("Test: Predictable halo with velocity (some high velocity)")
    print("="*60)
    
    np.random.seed(888)
    
    # 生成可预测的halo
    n_halo = 200
    n_noise = 400
    
    xp, expected_in_halo, halo_center, halo_range = generate_predictable_halo(
        n_halo_particles=n_halo,
        linking_length=LINKING_LENGTH,
        n_noise_particles=n_noise,
        seed=888
    )
    
    print(f"  Generated halo particles: {n_halo}")
    print(f"  Generated noise particles: {n_noise}")
    print(f"  Halo center: [{halo_center[0]:.2f}, {halo_center[1]:.2f}, {halo_center[2]:.2f}]")
    print(f"  Halo range: [{halo_range[0]:.2f}, {halo_range[1]:.2f}, {halo_range[2]:.2f}]")
    
    # 生成速度（基础热速度标准差 0.1，适配 a/(6*pi) 势阱深度）
    vp = np.random.randn(3, n_halo + n_noise).astype(np.float64) * 0.1
    
    # 给halo粒子中的10%设置极高的速度
    n_high_v = max(1, int(n_halo * 0.1))
    halo_indices = np.where(expected_in_halo)[0]
    high_v_indices = np.random.choice(halo_indices, n_high_v, replace=False)
    vp[:, high_v_indices] *= 100.0  # 100倍速度，保证高速粒子不束缚但不影响其他粒子
    
    print(f"  High velocity particles: {n_high_v}")
    print(f"  High velocity particle indices: {high_v_indices[:10]}{'...' if len(high_v_indices) > 10 else ''}")
    
    # 运行FOF（不带bound_check）
    in_halo_no_check = fof(xp, L_Box=L_BOX, linking_length=LINKING_LENGTH, np_halo_min=NP_HALO_MIN)
    
    # 运行FOF（带bound_check）
    in_halo_with_check, ang_mom = fof(xp, vp_list=vp, L_Box=L_BOX, linking_length=LINKING_LENGTH, 
                                     np_halo_min=NP_HALO_MIN, bound_check=True, verbose=1, return_ang_mom=True)
    print(f"    Calculated Angular Momentum: {ang_mom[:, 0]}")
    assert ang_mom.shape == (3, 1), f"Expected shape (3, 1), got {ang_mom.shape}"
    
    # 找出在FOF中但不束缚的粒子
    unbound_in_halo = in_halo_no_check & ~in_halo_with_check
    
    # 统计结果
    print(f"\n  Results:")
    print(f"    Particles in FOF (no bound check): {np.sum(in_halo_no_check)}")
    print(f"    Particles in FOF (with bound check): {np.sum(in_halo_with_check)}")
    print(f"    Unbound particles (removed by bound check): {np.sum(unbound_in_halo)}")
    
    # 验证高速粒子是否被正确排除
    high_v_removed = np.sum(unbound_in_halo[high_v_indices])
    print(f"    High velocity particles removed: {high_v_removed}/{n_high_v}")
    assert high_v_removed == n_high_v, f"Expected all {n_high_v} high-v particles removed, got {high_v_removed}"
    assert np.sum(in_halo_with_check) == (n_halo - n_high_v), f"Expected {n_halo - n_high_v} bound particles, got {np.sum(in_halo_with_check)}"
    
    # 绘图
    plot_particles(xp, in_halo_with_check, unbound_in_halo=unbound_in_halo,
                   title=f"Predictable Halo with Velocity (bound_check=True)",
                   save_path=os.path.join(OS_TEST_DIR, "test_predictable_halo_velocity.png"))
    
    return xp, in_halo_with_check, unbound_in_halo


# %%
# 测试多个可预测的halo
def test_multiple_predictable_halos():
    """测试多个可预测的halo"""
    print("\n" + "="*60)
    print("Test: Multiple predictable halos")
    print("="*60)
    
    np.random.seed(777)
    
    # 生成3个可预测的halo
    n_halos = 3
    n_halo_per_halo = 100
    n_noise_per_halo = 200  # 比例1:2
    
    xp_list = []
    expected_in_halo_list = []
    
    # 定义每个halo的偏移量，确保它们在空间上分开
    halo_shifts = [
        np.array([-50.0, -50.0, 0.0]),   # Halo 1: 左下
        np.array([50.0, -50.0, 0.0]),    # Halo 2: 右下
        np.array([0.0, 50.0, 30.0])      # Halo 3: 上方偏后
    ]
    
    for ihalo in range(n_halos):
        xp, expected_in_halo, halo_center, halo_range = generate_predictable_halo(
            n_halo_particles=n_halo_per_halo,
            linking_length=LINKING_LENGTH,
            n_noise_particles=n_noise_per_halo,
            seed=777 + ihalo
        )
        
        # 给整个数据集添加大偏移
        xp = xp + halo_shifts[ihalo].reshape(3, 1)
        
        xp_list.append(xp)
        expected_in_halo_list.append(expected_in_halo)
        print(f"  Halo {ihalo+1}: {n_halo_per_halo} halo particles, {n_noise_per_halo} noise particles")
    
    # 运行FOF
    in_halo = fof(xp_list, L_Box=L_BOX, linking_length=LINKING_LENGTH, np_halo_min=NP_HALO_MIN)
    
    # 统计结果
    total_particles = sum(xp.shape[1] for xp in xp_list)
    total_expected = sum(np.sum(expected) for expected in expected_in_halo_list)
    total_found = np.sum(in_halo)
    
    print(f"\n  Total particles: {total_particles}")
    print(f"  Expected in halos: {total_expected}")
    print(f"  Found in halos: {total_found}")
    
    # 计算准确率
    expected_all = np.hstack(expected_in_halo_list)
    correct = np.sum(in_halo == expected_all)
    print(f"  Per-particle accuracy: {correct}/{total_particles} ({100*correct/total_particles:.1f}%)")
    assert correct == total_particles, f"Expected 100% accuracy, got {correct}/{total_particles}"
    
    # 绘图
    xp_all = np.hstack(xp_list)
    plot_particles(xp_all, in_halo, 
                   title=f"Multiple Predictable Halos (expected {total_expected}, found {total_found})",
                   save_path=os.path.join(OS_TEST_DIR, "test_multiple_predictable_halos.png"))
    
    return xp_list, in_halo, expected_in_halo_list


# %%
# 测试多个可预测的halo带速度（部分粒子高速）
def test_multiple_predictable_halos_with_velocity():
    """测试多个可预测的halo带速度，部分粒子有高速会被bound_check排除"""
    print("\n" + "="*60)
    print("Test: Multiple predictable halos with velocity (some high velocity)")
    print("="*60)
    
    np.random.seed(666)
    
    # 生成3个可预测的halo
    n_halos = 3
    n_halo_per_halo = 100
    n_noise_per_halo = 200  # 比例1:2
    
    xp_list = []
    vp_list = []
    expected_in_halo_list = []
    high_v_indices_list = []
    
    # 定义每个halo的偏移量，确保它们在空间上分开
    halo_shifts = [
        np.array([-50.0, -50.0, 0.0]),   # Halo 1: 左下
        np.array([50.0, -50.0, 0.0]),    # Halo 2: 右下
        np.array([0.0, 50.0, 30.0])      # Halo 3: 上方偏后
    ]
    
    for ihalo in range(n_halos):
        # 生成可预测的halo
        xp, expected_in_halo, halo_center, halo_range = generate_predictable_halo(
            n_halo_particles=n_halo_per_halo,
            linking_length=LINKING_LENGTH,
            n_noise_particles=n_noise_per_halo,
            seed=666 + ihalo
        )
        
        # 给整个数据集添加大偏移，防止不同halo连在一起
        xp = xp + halo_shifts[ihalo].reshape(3, 1)
        
        # 生成速度（基础热速度标准差 0.1，适配 a/(6*pi) 势阱深度）
        vp = np.random.randn(3, n_halo_per_halo + n_noise_per_halo).astype(np.float64) * 0.1
        
        # 给halo粒子中的10%设置极高的速度
        n_high_v = max(1, int(n_halo_per_halo * 0.1))
        halo_indices = np.where(expected_in_halo)[0]
        high_v_indices = np.random.choice(halo_indices, n_high_v, replace=False)
        vp[:, high_v_indices] *= 100.0  # 100倍速度，保证高速粒子不束缚
        
        xp_list.append(xp)
        vp_list.append(vp)
        expected_in_halo_list.append(expected_in_halo)
        high_v_indices_list.append(high_v_indices)
        
        print(f"  Halo {ihalo+1}: {n_halo_per_halo} halo particles, {n_noise_per_halo} noise particles, {n_high_v} high velocity")
    
    # 运行FOF（不带bound_check）
    in_halo_no_check = fof(xp_list, L_Box=L_BOX, linking_length=LINKING_LENGTH, np_halo_min=NP_HALO_MIN)
    
    # 运行FOF（带bound_check）
    in_halo_with_check = fof(xp_list, vp_list=vp_list, L_Box=L_BOX, linking_length=LINKING_LENGTH, 
                             np_halo_min=NP_HALO_MIN, bound_check=True, verbose=1)
    
    # 找出在FOF中但不束缚的粒子
    unbound_in_halo = in_halo_no_check & ~in_halo_with_check
    
    # 统计结果
    total_particles = sum(xp.shape[1] for xp in xp_list)
    total_expected = sum(np.sum(expected) for expected in expected_in_halo_list)
    
    print(f"\n  Results:")
    print(f"    Total particles: {total_particles}")
    print(f"    Expected in halos: {total_expected}")
    print(f"    Particles in FOF (no bound check): {np.sum(in_halo_no_check)}")
    print(f"    Particles in FOF (with bound check): {np.sum(in_halo_with_check)}")
    print(f"    Unbound particles (removed by bound check): {np.sum(unbound_in_halo)}")
    
    # 验证高速粒子是否被正确排除
    total_high_v_removed = 0
    total_high_v = 0
    for ihalo in range(n_halos):
        high_v_indices = high_v_indices_list[ihalo]
        # 计算该halo在总粒子数组中的起始索引
        start_idx = sum(xp.shape[1] for xp in xp_list[:ihalo])
        # 调整high_v_indices为全局索引
        global_high_v_indices = high_v_indices + start_idx
        high_v_removed = np.sum(unbound_in_halo[global_high_v_indices])
        total_high_v_removed += high_v_removed
        total_high_v += len(high_v_indices)
    
    print(f"    High velocity particles removed: {total_high_v_removed}/{total_high_v}")
    assert total_high_v_removed == total_high_v, f"Expected all {total_high_v} high-v particles removed, got {total_high_v_removed}"
    expected_bound = total_expected - total_high_v
    assert np.sum(in_halo_with_check) == expected_bound, f"Expected {expected_bound} bound particles, got {np.sum(in_halo_with_check)}"
    
    # 绘图
    xp_all = np.hstack(xp_list)
    plot_particles(xp_all, in_halo_with_check, unbound_in_halo=unbound_in_halo,
                   title=f"Multiple Predictable Halos with Velocity (bound_check=True)",
                   save_path=os.path.join(OS_TEST_DIR, "test_multiple_predictable_halos_velocity.png"))
    
    return xp_list, in_halo_with_check, unbound_in_halo


# %%
# def test_multiple_halos_with_velocity():
#     """测试多个halo带速度，部分粒子有高速会被bound_check排除"""
#     print("\n" + "="*60)
#     print("Test 6: Multiple halos with velocity (some high velocity)")
#     print("="*60)
#     
#     # 生成多个halo，每个halo有部分高速粒子
#     xp_list, vp_list, expected_in_halo_list = generate_multiple_halos_with_velocity(
#         n_halos=3,
#         n_halo_particles_per_halo=100,
#         linking_length=LINKING_LENGTH,
#         velocity_scale=0.5,
#         high_velocity_fraction=0.1,  # 10%的粒子有高速
#         seed=1234
#     )
#     
#     print(f"  Number of halos: {len(xp_list)}")
#     for i, (xp, vp, expected) in enumerate(zip(xp_list, vp_list, expected_in_halo_list)):
#         n_expected = np.sum(expected)
#         n_high_v = len(expected) - n_expected
#         print(f"  Halo {i+1}: {xp.shape[1]} particles, {n_expected} expected in halo, {n_high_v} high velocity")
#     
#     # 运行FOF（不带bound_check）
#     in_halo_no_check = fof(xp_list, L_Box=L_BOX, linking_length=LINKING_LENGTH, np_halo_min=NP_HALO_MIN)
#     
#     # 运行FOF（带bound_check）
#     in_halo_with_check = fof(xp_list, vp_list=vp_list, L_Box=L_BOX, linking_length=LINKING_LENGTH, 
#                              np_halo_min=NP_HALO_MIN, bound_check=True)
#     
#     # 统计结果
#     total_particles = sum(xp.shape[1] for xp in xp_list)
#     total_expected = sum(np.sum(expected) for expected in expected_in_halo_list)
#     
#     print(f"\n  Total particles: {total_particles}")
#     print(f"  Expected in halos (excluding high velocity): {total_expected}")
#     print(f"  Found in halos (no bound check): {np.sum(in_halo_no_check)}")
#     print(f"  Found in halos (with bound check): {np.sum(in_halo_with_check)}")
#     
#     # 绘图
#     xp_all = np.hstack(xp_list)
#     plot_particles(xp_all, in_halo_with_check, 
#                    title=f"Multiple Halos with Velocity (bound_check=True)",
#                    save_path="test_multiple_halos_velocity.png")
#     
#     return xp_list, vp_list, in_halo_with_check


# %% 新增功能专项测试
# 新增功能专项测试
def test_new_features():
    """
    新增功能专项验证测试：
    1. fof_with_info 的 O(N) 准确性与一致性
    2. 不同 n_threads (1, 2, 4, None) 并行一致性
    3. 标量 vs 数组等质量一致性及变质量束缚测试
    4. 双单位模式 (units='cube' vs 'physical') 及哈勃流测试
    5. N > 1000 候选跨度 > 15 Mpc 时的安全 ValueError 异常捕获
    """
    print("\n" + "="*60)
    print("Test: New Features Comprehensive Verification")
    print("="*60)

    # 1. fof_with_info 的 O(N) 准确性
    print("\n[1/5] Verifying fof_with_info O(N) accuracy...")
    np.random.seed(101)
    n_h = 3
    xp_list = []
    vp_list = []
    for i in range(n_h):
        xp, exp, _, _ = generate_predictable_halo(100, LINKING_LENGTH, 50, seed=101 + i)
        xp = xp + np.array([i * 100.0, 0.0, 0.0]).reshape(3, 1)
        vp = np.random.randn(3, 150).astype(np.float64) * 0.05
        xp_list.append(xp)
        vp_list.append(vp)

    res_info = fof_with_info(xp_list, vp_list=vp_list, L_Box=L_BOX, linking_length=LINKING_LENGTH,
                             np_halo_min=NP_HALO_MIN, bound_check=True, return_ang_mom=True)
    res_fof, ang_mom = fof(xp_list, vp_list=vp_list, L_Box=L_BOX, linking_length=LINKING_LENGTH,
                           np_halo_min=NP_HALO_MIN, bound_check=True, return_ang_mom=True)

    assert np.array_equal(res_info['in_halo'], res_fof), "fof_with_info in_halo must match fof in_halo"
    assert res_info['nhalo_found'] == n_h, f"Expected {n_h} halos, got {res_info['nhalo_found']}"
    assert len(res_info['halo_indices']) == n_h, "halo_indices length mismatch"
    assert np.sum(res_info['np_per_halo']) == np.sum(res_info['in_halo']), "Sum of np_per_halo must equal total bound particles"
    all_indices = np.concatenate(res_info['halo_indices']) if len(res_info['halo_indices']) > 0 else np.array([], dtype=int)
    assert len(all_indices) == len(np.unique(all_indices)), "halo_indices must be disjoint"
    assert np.all(res_info['in_halo'][all_indices]), "All indices in halo_indices must have in_halo=True"
    assert np.allclose(res_info['ang_mom'], ang_mom), "ang_mom mismatch in fof_with_info"
    print("  ✓ fof_with_info passed with exact O(N) halo membership matching fof().")

    # 2. n_threads (1, 2, 4, None) 并行一致性
    print("\n[2/5] Verifying n_threads consistency (1, 2, 4, None)...")
    res_threads = {}
    ang_threads = {}
    for nt in [1, 2, 4, None]:
        inh, am = fof(xp_list, vp_list=vp_list, L_Box=L_BOX, linking_length=LINKING_LENGTH,
                      np_halo_min=NP_HALO_MIN, bound_check=True, n_threads=nt, return_ang_mom=True)
        res_threads[nt] = inh
        ang_threads[nt] = am

    base_inh = res_threads[1]
    base_am = ang_threads[1]
    for nt in [2, 4, None]:
        assert np.array_equal(res_threads[nt], base_inh), f"n_threads={nt} in_halo differs from n_threads=1"
        assert np.allclose(ang_threads[nt], base_am, atol=1e-12), f"n_threads={nt} ang_mom differs from n_threads=1"

    # 测试超限线程数异常
    try:
        fof(xp_list, n_threads=300)
        assert False, "Should have raised ValueError for n_threads > 256"
    except ValueError as e:
        print(f"  ✓ Caught expected exception for n_threads > 256: {e}")
    print("  ✓ Multi-threading (1, 2, 4, None) produces identical bound halo results.")

    # 3. 质量参数 mp_list: 标量 vs 数组等质量一致性及变质量束缚测试
    print("\n[3/5] Verifying mp_list (scalar vs array and variable mass effects)...")
    total_p = sum(x.shape[1] for x in xp_list)
    # 标量 1.0 vs 全 1 数组
    inh_scalar = fof(xp_list, vp_list=vp_list, L_Box=L_BOX, linking_length=LINKING_LENGTH,
                     np_halo_min=NP_HALO_MIN, bound_check=True, mp_list=1.0)
    inh_array = fof(xp_list, vp_list=vp_list, L_Box=L_BOX, linking_length=LINKING_LENGTH,
                    np_halo_min=NP_HALO_MIN, bound_check=True, mp_list=np.ones(total_p, dtype=np.float64))
    assert np.array_equal(inh_scalar, inh_array), "Scalar mp_list and equivalent array mp_list must match"

    # 变质量测试：构造单 Halo，给定较高的速度，如果等质量时部分粒子逃逸，但中心增加超大质量后，势阱加深粒子被束缚
    np.random.seed(202)
    xp_single, exp_s, _, _ = generate_predictable_halo(100, LINKING_LENGTH, 0, seed=202)
    # 给定临界速度
    vp_single = np.random.randn(3, 100).astype(np.float64) * 0.25
    inh_light = fof(xp_single, vp_list=vp_single, L_Box=L_BOX, linking_length=LINKING_LENGTH,
                    np_halo_min=NP_HALO_MIN, bound_check=True, mp_list=1.0)
    # 增加中心粒子的质量到 200 倍
    mp_heavy = np.ones(100, dtype=np.float64)
    mp_heavy[0] = 200.0
    inh_heavy = fof(xp_single, vp_list=vp_single, L_Box=L_BOX, linking_length=LINKING_LENGTH,
                    np_halo_min=NP_HALO_MIN, bound_check=True, mp_list=mp_heavy)
    assert np.sum(inh_heavy) >= np.sum(inh_light), "Deeper potential well from heavy center should bind equal or more particles"
    print(f"  ✓ Variable mass test passed: light bound = {np.sum(inh_light)}, heavy center bound = {np.sum(inh_heavy)}.")

    # 4. 双单位模式 (units='cube' vs 'physical') 及哈勃流测试
    print("\n[4/5] Verifying units='cube' vs 'physical' and hubble flow...")
    # units='physical' 缺少 softening 或质量时应抛出异常
    try:
        fof(xp_single, vp_list=vp_single, units='physical', bound_check=True)
        assert False, "Should raise ValueError when softening is missing in units='physical'"
    except ValueError as e:
        print(f"  ✓ Caught expected exception for missing softening in units='physical': {e}")

    # 构造等价物理单位测试
    ng = 512.0
    box_mpc = 300.0
    dx = box_mpc / ng
    xp_phys = xp_single * dx
    soft_phys = 0.1 * dx
    om_val = 0.3
    a_val = 1.0  # z = 0.0
    v_factor = 150.0 * np.sqrt(om_val) * box_mpc / (a_val * ng)
    vp_phys = vp_single * v_factor
    m_grid = om_val * 2.77536627e11 * (box_mpc ** 3) / (ng ** 3)
    mp_phys = mp_heavy * m_grid

    inh_cube, am_cube = fof(xp_single, vp_list=vp_single, L_Box=ng, linking_length=LINKING_LENGTH,
                            np_halo_min=NP_HALO_MIN, bound_check=True, mp_list=mp_heavy,
                            softening=0.1, units='cube', return_ang_mom=True)
    inh_phys, am_phys = fof(xp_phys, vp_list=vp_phys, L_Box=box_mpc, linking_length=LINKING_LENGTH * dx,
                            np_halo_min=NP_HALO_MIN, bound_check=True, mp_list=mp_phys,
                            softening=soft_phys, units='physical', return_ang_mom=True)
    assert np.array_equal(inh_cube, inh_phys), "units='cube' and equivalent units='physical' must produce identical in_halo"
    assert np.allclose(am_cube * (m_grid * dx * v_factor), am_phys, rtol=1e-5), "Converted physical angular momentum mismatch"
    print("  ✓ units='cube' and units='physical' equivalence verified.")

    # 验证 include_hubble_flow
    inh_hf = fof(xp_phys, vp_list=vp_phys, L_Box=box_mpc, linking_length=LINKING_LENGTH * dx,
                 np_halo_min=NP_HALO_MIN, bound_check=True, mp_list=mp_phys,
                 softening=soft_phys, units='physical', include_hubble_flow=True)
    assert isinstance(inh_hf, np.ndarray), "include_hubble_flow=True should return valid boolean array"
    print("  ✓ include_hubble_flow=True verified.")

    # 5. N > 1000 候选跨度 > 15 Mpc 时的安全 ValueError 异常捕获
    print("\n[5/5] Verifying 15 Mpc span limit ValueError for N > 1000...")
    # 构造 N = 1050 粒子，沿 X 轴分布，总跨度为 25.0 格点（> s_15 ≈ 17.152 格点）
    # 相邻间隔 25.0 / 1050 ≈ 0.0238 格点 << linking_length=2.0，必定连成一个 Halo 候选
    np_large = 1050
    x_coords = np.linspace(10.0, 35.0, np_large, endpoint=True)  # span = 25.0 > 17.152
    y_coords = np.zeros(np_large)
    z_coords = np.zeros(np_large)
    xp_overspan = np.vstack([x_coords, y_coords, z_coords])
    vp_overspan = np.zeros((3, np_large), dtype=np.float64)

    try:
        fof(xp_overspan, vp_list=vp_overspan, L_Box=L_BOX, linking_length=LINKING_LENGTH,
            np_halo_min=NP_HALO_MIN, bound_check=True, box=300.0, h0=0.67)
        assert False, "Should raise ValueError when candidate halo span exceeds 15 Mpc"
    except ValueError as e:
        print(f"  ✓ Successfully caught expected ValueError without crashing Python: {e}")
        assert "Candidate halo span exceeds 15 Mpc limit" in str(e) or "exceeds" in str(e)
    print("  ✓ 15 Mpc span safety limit verified.")


# %%
# 主函数
def main():
    """运行所有测试"""
    print("\n" + "#"*60)
    print("#  py_fof Test Suite")
    print("#"*60)
    print(f"\nGlobal Parameters:")
    print(f"  N_PARTICLES    = {N_PARTICLES}")
    print(f"  BOX_SIZE       = {BOX_SIZE}")
    print(f"  L_BOX          = {L_BOX}")
    print(f"  LINKING_LENGTH = {LINKING_LENGTH}")
    print(f"  NP_HALO_MIN    = {NP_HALO_MIN}")
    print(f"  VELOCITY_SCALE = {VELOCITY_SCALE}")
    
    try:
        # test_single_halo()
        # test_multiple_halos()
        # test_with_binding_check()
        # test_detailed_output()
        test_predictable_halo()
        test_predictable_halo_with_velocity()
        test_multiple_predictable_halos()
        test_multiple_predictable_halos_with_velocity()
        # test_multiple_halos_with_velocity()
        
        # 专项新功能验证
        test_new_features()
        
        # 性能测试
        test_performance()
        
        print("\n" + "#"*60)
        print("#  All tests completed successfully!")
        print("#  Figures saved:")
        print("#    .test_ai/test_predictable_halo.png")
        print("#    .test_ai/test_predictable_halo_velocity.png")
        print("#    .test_ai/test_multiple_predictable_halos.png")
        print("#    .test_ai/test_multiple_predictable_halos_velocity.png")
        print("#"*60 + "\n")
        
    except Exception as e:
        print(f"\nError during testing: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
