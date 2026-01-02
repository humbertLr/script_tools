#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
将VASP的XDATCAR文件或LAMMPS的data文件转换为包含晶格信息的XYZ格式轨迹文件

XYZ文件的注释行将包含完整的晶格信息,可以被OVITO等软件正确读取

使用方法:
    python trajectory_converter.py XDATCAR output.xyz
    python trajectory_converter.py input.data output.xyz
    python trajectory_converter.py input.data output.xyz --format lammps-data
"""

import sys
import numpy as np
from ase.io import read
import argparse


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


def detect_file_format(input_file):
    """
    自动检测输入文件格式
    
    返回:
        'vasp-xdatcar' 或 'lammps-data'
    """
    with open(input_file, 'r') as f:
        first_lines = [f.readline().strip() for _ in range(10)]
    
    # 检查是否是LAMMPS data文件
    lammps_keywords = ['atoms', 'atom types', 'xlo xhi', 'ylo yhi', 'zlo zhi', 
                       'xy xz yz', 'Atoms', 'Masses', 'Velocities']
    lammps_score = sum(1 for line in first_lines for kw in lammps_keywords if kw in line)
    
    # 检查是否是XDATCAR文件
    # XDATCAR通常第6行是元素名称,第7行是原子数
    xdatcar_pattern = False
    if len(first_lines) >= 7:
        try:
            # 尝试解析第7行是否为整数(原子数)
            numbers = first_lines[6].split()
            if all(num.isdigit() for num in numbers):
                xdatcar_pattern = True
        except:
            pass
    
    if lammps_score >= 3:
        return 'lammps-data'
    elif xdatcar_pattern or 'Direct configuration' in ' '.join(first_lines):
        return 'vasp-xdatcar'
    else:
        # 默认尝试VASP格式
        return 'vasp-xdatcar'


def convert_to_xyz(input_file, xyz_file='trajectory.xyz', file_format=None):
    """
    将输入文件转换为包含晶格信息的XYZ格式
    
    参数:
        input_file: 输入文件路径 (XDATCAR或LAMMPS data)
        xyz_file: 输出的XYZ文件路径
        file_format: 文件格式 ('vasp-xdatcar', 'lammps-data', 或 None自动检测)
    """
    try:
        # 自动检测文件格式
        if file_format is None:
            detected_format = detect_file_format(input_file)
            print(f"自动检测到文件格式: {detected_format}")
            file_format = detected_format
        
        # 读取输入文件
        print(f"正在读取 {input_file} (格式: {file_format})...")
        
        if file_format == 'vasp-xdatcar':
            trajectory = read(input_file, index=':', format='vasp-xdatcar')
        elif file_format == 'lammps-data':
            # LAMMPS data文件通常只包含一个结构
            trajectory = read(input_file, format='lammps-data', style='atomic')
        else:
            raise ValueError(f"不支持的文件格式: {file_format}")
        
        # 确保trajectory是列表
        if not isinstance(trajectory, list):
            trajectory = [trajectory]
        
        n_frames = len(trajectory)
        print(f"读取到 {n_frames} 帧结构")
        
        # 写入包含晶格信息的XYZ文件
        print(f"正在写入 {xyz_file} (包含晶格信息)...")
        write_xyz_with_lattice(xyz_file, trajectory)
        
        print(f"\n✓ 转换完成!")
        print(f"  输入文件: {input_file}")
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
        print(f"  体积: {cell.volume:.4f} Ų")
        
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
        print(f"错误: 找不到文件 {input_file}")
        sys.exit(1)
    except Exception as e:
        print(f"错误: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='将VASP XDATCAR或LAMMPS data文件转换为包含晶格信息的XYZ格式',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 转换VASP XDATCAR文件 (自动检测)
  python trajectory_converter.py XDATCAR output.xyz
  
  # 转换LAMMPS data文件 (自动检测)
  python trajectory_converter.py input.data output.xyz
  
  # 手动指定格式
  python trajectory_converter.py input.data output.xyz --format lammps-data
  python trajectory_converter.py XDATCAR output.xyz --format vasp-xdatcar
        """
    )
    
    parser.add_argument('input_file', help='输入文件 (XDATCAR或LAMMPS data)')
    parser.add_argument('output_file', nargs='?', default=None,
                       help='输出XYZ文件 (默认: trajectory.xyz)')
    parser.add_argument('--format', '-f', choices=['vasp-xdatcar', 'lammps-data'],
                       help='指定输入文件格式 (默认: 自动检测)')
    
    args = parser.parse_args()
    
    # 设置默认输出文件名
    if args.output_file is None:
        if args.input_file.endswith('.data'):
            args.output_file = args.input_file.replace('.data', '.xyz')
        else:
            args.output_file = 'trajectory.xyz'
    
    convert_to_xyz(args.input_file, args.output_file, args.format)


if __name__ == '__main__':
    main()