#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
将VASP的XDATCAR文件转换为包含晶格信息的XYZ格式轨迹文件

XYZ文件的注释行将包含完整的晶格信息,可以被OVITO等软件正确读取

使用方法:
    python xdatcar_to_xyz_with_lattice.py XDATCAR output.xyz
"""

import sys
import numpy as np
from ase.io import read


def write_xyz_with_lattice(filename, trajectory):
    """
    写入包含晶格信息的XYZ文件
    
    格式遵循extended XYZ格式标准,在注释行中包含Lattice信息
    """
    with open(filename, 'w') as f:
        for i, atoms in enumerate(trajectory):
            # 写入原子数
            f.write(f"{len(atoms)}\n")
            
            # 写入注释行,包含晶格信息
            cell = atoms.get_cell()
            # 将晶格向量转换为字符串格式
            lattice_str = ' '.join([f'{x:.10f}' for row in cell for x in row])
            
            # 获取其他信息
            formula = atoms.get_chemical_formula()
            
            # 写入extended XYZ格式的注释行
            # 格式: Lattice="a1 a2 a3 b1 b2 b3 c1 c2 c3" Properties=species:S:1:pos:R:3
            comment = f'Lattice="{lattice_str}" Properties=species:S:1:pos:R:3'
            
            # 如果有能量信息也可以添加
            if atoms.calc is not None and hasattr(atoms.calc, 'results'):
                if 'energy' in atoms.calc.results:
                    energy = atoms.calc.results['energy']
                    comment += f' energy={energy:.10f}'
            
            f.write(f"{comment}\n")
            
            # 写入原子坐标
            symbols = atoms.get_chemical_symbols()
            positions = atoms.get_positions()
            
            for symbol, pos in zip(symbols, positions):
                f.write(f"{symbol:2s}  {pos[0]:18.10f}  {pos[1]:18.10f}  {pos[2]:18.10f}\n")


def xdatcar_to_xyz_with_lattice(xdatcar_file, xyz_file='trajectory.xyz'):
    """
    将XDATCAR文件转换为包含晶格信息的XYZ格式
    
    参数:
        xdatcar_file: XDATCAR文件路径
        xyz_file: 输出的XYZ文件路径
    """
    try:
        # 读取XDATCAR文件
        print(f"正在读取 {xdatcar_file}...")
        trajectory = read(xdatcar_file, index=':', format='vasp-xdatcar')
        
        # 确保trajectory是列表
        if not isinstance(trajectory, list):
            trajectory = [trajectory]
        
        n_frames = len(trajectory)
        print(f"读取到 {n_frames} 帧结构")
        
        # 写入包含晶格信息的XYZ文件
        print(f"正在写入 {xyz_file} (包含晶格信息)...")
        write_xyz_with_lattice(xyz_file, trajectory)
        
        print(f"\n✓ 转换完成!")
        print(f"  输出文件: {xyz_file}")
        print(f"  帧数: {n_frames}")
        
        # 显示第一帧的信息
        first_frame = trajectory[0]
        print(f"\n第一帧信息:")
        print(f"  原子数: {len(first_frame)}")
        print(f"  化学式: {first_frame.get_chemical_formula()}")
        
        # 显示晶胞信息
        cell = first_frame.get_cell()
        print(f"\n晶胞信息:")
        print(f"  a向量: {cell[0]}")
        print(f"  b向量: {cell[1]}")
        print(f"  c向量: {cell[2]}")
        print(f"  晶胞参数: a={cell.lengths()[0]:.4f} Å, "
              f"b={cell.lengths()[1]:.4f} Å, "
              f"c={cell.lengths()[2]:.4f} Å")
        angles = cell.angles()
        print(f"  晶胞角度: α={angles[0]:.2f}°, "
              f"β={angles[1]:.2f}°, "
              f"γ={angles[2]:.2f}°")
        print(f"  体积: {cell.volume:.4f} Å³")
        
        # 显示元素组成
        from collections import Counter
        symbols = first_frame.get_chemical_symbols()
        element_counts = Counter(symbols)
        print(f"\n元素组成:")
        for element, count in sorted(element_counts.items()):
            print(f"  {element}: {count}")
        
        print(f"\n注意: XYZ文件的注释行包含完整的晶格信息")
        print(f"      格式: Lattice=\"a1 a2 a3 b1 b2 b3 c1 c2 c3\"")
        print(f"      可被OVITO, ASE等软件正确读取")
        
    except FileNotFoundError:
        print(f"错误: 找不到文件 {xdatcar_file}")
        sys.exit(1)
    except Exception as e:
        print(f"错误: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def main():
    """主函数"""
    if len(sys.argv) < 2:
        print("用法: python xdatcar_to_xyz_with_lattice.py XDATCAR [output.xyz]")
        print("示例: python xdatcar_to_xyz_with_lattice.py XDATCAR trajectory.xyz")
        sys.exit(1)
    
    xdatcar_file = sys.argv[1]
    xyz_file = sys.argv[2] if len(sys.argv) > 2 else 'trajectory_with_lattice.xyz'
    
    xdatcar_to_xyz_with_lattice(xdatcar_file, xyz_file)


if __name__ == '__main__':
    main()