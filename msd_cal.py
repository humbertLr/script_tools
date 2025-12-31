#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VASP XDATCAR MSD计算脚本
从AIMD的XDATCAR文件中读取轨迹数据并计算均方位移(MSD)
fasfa dsaf
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# 设置中文字体
# rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
# rcParams['axes.unicode_minus'] = False

def read_xdatcar(filename='XDATCAR'):
    """
    读取VASP的XDATCAR文件
    
    参数:
        filename: XDATCAR文件路径
    
    返回:
        lattice: 晶格矩阵 (3x3)
        atom_types: 原子种类列表
        atom_numbers: 各种原子的数量
        positions: 所有时间步的原子位置 (n_steps, n_atoms, 3)
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # 读取晶格常数
    scale = float(lines[1].strip())
    lattice = np.array([[float(x) for x in lines[i].split()] for i in range(2, 5)])
    lattice *= scale
    
    # 读取原子种类和数量
    atom_types = lines[5].split()
    atom_numbers = [int(x) for x in lines[6].split()]
    n_atoms = sum(atom_numbers)
    
    # 读取所有时间步的坐标
    positions = []
    i = 7
    while i < len(lines):
        if 'Direct configuration=' in lines[i] or 'configuration=' in lines[i]:
            coords = []
            for j in range(i+1, i+1+n_atoms):
                coord = [float(x) for x in lines[j].split()[:3]]
                coords.append(coord)
            positions.append(coords)
            i += n_atoms + 1
        else:
            i += 1
    
    positions = np.array(positions)
    
    print(f"成功读取 {len(positions)} 个时间步")
    print(f"原子总数: {n_atoms}")
    print(f"原子种类: {atom_types}")
    print(f"各种原子数量: {atom_numbers}")
    
    return lattice, atom_types, atom_numbers, positions

def fractional_to_cartesian(positions, lattice):
    """
    将分数坐标转换为笛卡尔坐标
    
    参数:
        positions: 分数坐标 (n_steps, n_atoms, 3)
        lattice: 晶格矩阵 (3x3)
    
    返回:
        cartesian: 笛卡尔坐标 (n_steps, n_atoms, 3)
    """
    return np.dot(positions, lattice)

def unwrap_positions(positions, lattice):
    """
    处理周期性边界条件,展开原子轨迹
    
    参数:
        positions: 分数坐标 (n_steps, n_atoms, 3)
        lattice: 晶格矩阵 (3x3)
    
    返回:
        unwrapped: 展开后的笛卡尔坐标 (n_steps, n_atoms, 3)
    """
    n_steps, n_atoms, _ = positions.shape
    unwrapped = np.zeros_like(positions)
    unwrapped[0] = positions[0]
    
    for i in range(1, n_steps):
        delta = positions[i] - positions[i-1]
        # 处理跨越周期性边界的情况
        delta = delta - np.round(delta)
        unwrapped[i] = unwrapped[i-1] + delta
    
    # 转换为笛卡尔坐标
    unwrapped_cart = fractional_to_cartesian(unwrapped, lattice)
    
    return unwrapped_cart

def calculate_msd(positions):
    """
    计算均方位移(MSD)
    
    参数:
        positions: 笛卡尔坐标 (n_steps, n_atoms, 3)
    
    返回:
        msd: 各时间间隔的MSD值
    """
    n_steps, n_atoms, _ = positions.shape
    msd = np.zeros(n_steps)
    
    for dt in range(n_steps):
        displacements = []
        for t0 in range(n_steps - dt):
            displacement = positions[t0 + dt] - positions[t0]
            displacements.append(np.sum(displacement**2, axis=1))
        msd[dt] = np.mean(displacements)
    
    return msd

def calculate_msd_by_species(positions, atom_numbers):
    """
    分别计算每种原子的MSD
    
    参数:
        positions: 笛卡尔坐标 (n_steps, n_atoms, 3)
        atom_numbers: 各种原子的数量列表
    
    返回:
        msd_list: 各种原子的MSD列表
    """
    msd_list = []
    start_idx = 0
    
    for num in atom_numbers:
        end_idx = start_idx + num
        species_positions = positions[:, start_idx:end_idx, :]
        msd = calculate_msd(species_positions)
        msd_list.append(msd)
        start_idx = end_idx
    
    return msd_list

def calculate_diffusion_coefficient(msd, time, fit_range=None):
    """
    从MSD计算扩散系数 D = MSD / (6t)
    
    参数:
        msd: MSD数组
        time: 时间数组
        fit_range: 拟合范围 (start_idx, end_idx)
    
    返回:
        D: 扩散系数 (Å²/ps)
    """
    if fit_range is None:
        # 使用后半段数据进行拟合
        fit_range = (len(msd)//2, len(msd))
    
    start, end = fit_range
    coeffs = np.polyfit(time[start:end], msd[start:end], 1)
    D = coeffs[0] / 6.0  # MSD = 6Dt
    
    return D, coeffs

def plot_msd(time, msd_list, atom_types, timestep=20, save_fig=True):
    """
    绘制MSD曲线
    
    参数:
        time: 时间数组
        msd_list: MSD列表
        atom_types: 原子种类列表
        timestep: 时间步长 (fs)
        save_fig: 是否保存图片
    """
    plt.figure(figsize=(10, 6))
    
    colors = plt.cm.tab10(np.linspace(0, 1, len(msd_list)))
    
    for i, (msd, atom_type) in enumerate(zip(msd_list, atom_types)):
        # 计算扩散系数
        D, coeffs = calculate_diffusion_coefficient(msd, time)
        
        plt.plot(time, msd, 'o-', markersize=3, color=colors[i], 
                label=f'{atom_type} (D={D:.4e} Å²/ps)', linewidth=2)
        
        # 绘制拟合直线
        fit_start = len(msd)//2
        fit_line = coeffs[0] * time + coeffs[1]
        plt.plot(time[fit_start:], fit_line[fit_start:], '--', 
                color=colors[i], alpha=0.5, linewidth=1.5)
    
    plt.xlabel('time (ps)', fontsize=12)
    plt.ylabel('MSD (Ų)', fontsize=12)
    plt.title('MSD vs time', fontsize=14, fontweight='bold')
    plt.legend(fontsize=10, frameon=True)
    plt.grid(True, alpha=0.3, linestyle='--')
    plt.tight_layout()
    
    if save_fig:
        plt.savefig('msd_plot.png', dpi=300, bbox_inches='tight')
        print("图片已保存为 msd_plot.png")
    
    plt.show()

def main():
    """主函数"""
    # 参数设置
    xdatcar_file = 'XDATCAR'  # XDATCAR文件路径
    timestep = 1.0  # 时间步长 (fs)
    
    print("="*50)
    print("VASP AIMD MSD计算程序")
    print("="*50)
    
    # 读取XDATCAR文件
    print("\n正在读取XDATCAR文件...")
    lattice, atom_types, atom_numbers, positions = read_xdatcar(xdatcar_file)
    
    # 展开周期性边界
    print("\n正在处理周期性边界条件...")
    unwrapped_positions = unwrap_positions(positions, lattice)
    
    # 计算MSD
    print("\n正在计算MSD...")
    msd_list = calculate_msd_by_species(unwrapped_positions, atom_numbers)
    
    # 计算总体MSD
    total_msd = calculate_msd(unwrapped_positions)
    msd_list.append(total_msd)
    atom_types_plot = atom_types + ['Total']
    
    # 生成时间数组 (转换为ps)
    n_steps = len(positions)
    time = np.arange(n_steps) * timestep / 1000.0  # fs -> ps
    
    # 输出扩散系数
    print("\n扩散系数:")
    print("-"*50)
    for i, atom_type in enumerate(atom_types_plot):
        D, _ = calculate_diffusion_coefficient(msd_list[i], time)
        print(f"{atom_type:10s}: D = {D:.4e} Ų/ps = {D*1e-4:.4e} cm²/s")
    print("-"*50)
    
    # 绘制MSD图
    print("\n正在绘制MSD曲线...")
    plot_msd(time, msd_list, atom_types_plot, timestep)
    
    # 保存数据
    print("\n正在保存数据...")
    data = np.column_stack([time] + msd_list)
    header = 'Time(ps) ' + ' '.join([f'MSD_{at}(A^2)' for at in atom_types_plot])
    np.savetxt('msd_data.txt', data, header=header, fmt='%.6e')
    print("数据已保存为 msd_data.txt")
    
    print("\n计算完成!")
    print("="*50)

if __name__ == '__main__':
    main()