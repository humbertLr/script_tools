# -*- coding: utf-8 -*-
"""
Convert material density (g/cm³) to atomic density (atoms/Å³)
Author: Xinyu
"""

# ---------- 输入参数 ----------
density_min = 10.022  # 下限密度 g/cm³  source: https://www.sciencedirect.com/science/article/pii/002219026380233X
density_max = 10.05   # 上限密度 g/cm³  source: wikipedia; https://moltensalt.org/references/static/downloads/pdf/element-salt-densities.pdf
M = 208.980            # 摩尔质量 g/mol (例如 Bi)
NA = 6.022e23         # 阿伏伽德罗常数 atoms/mol

# ---------- 转换公式 ----------
def density_to_atom_density(rho, M, NA):
    """
    Convert density in g/cm³ to atomic density in atoms/Å³
    """
    return (rho * NA) / (M * 1e24)

# ---------- 计算 ----------
atom_density_min = density_to_atom_density(density_min, M, NA)
atom_density_max = density_to_atom_density(density_max, M, NA)

# ---------- 输出 ----------
print(f"Density range: {density_min} - {density_max} g/cm³")
print(f"Atomic density range: {atom_density_min:.7f} - {atom_density_max:.7f} atoms/Å³")

'''
Density range: 10.022 - 10.05 g/cm³
Atomic density range: 0.0288796 - 0.0289602 atoms/Å³
'''
