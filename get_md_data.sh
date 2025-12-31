#!/bin/bash

# 检查是否存在 OSZICAR 文件
if [ ! -f OSZICAR ]; then
    echo "错误: 当前目录下未找到 OSZICAR 文件。"
    exit 1
fi

echo "正在处理 OSZICAR..."

file1="temp.txt"
file2="energy.txt"
timestep=0.002

# 提取数据
# grep "T=" 筛选出包含温度和能量信息的行
# awk '{print $1, $3}' 提取第1列(步数)和第3列(温度)
# awk '{print $1, $5}' 提取第1列(步数)和第5列(能量E)

# 输出步数和温度到 temperature.dat
grep "T=" OSZICAR | awk -v dt=$timestep '{print $1*dt, $3}' > $file1

# 输出步数和能量到 energy.dat
grep "T=" OSZICAR | awk -v dt=$timestep '{print $1*dt, $5}' > $file2

echo "处理完成！"
echo "生成文件: $file1 (time vs temp)"
echo "生成文件: $file2 (time vs energy)"
