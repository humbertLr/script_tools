#!/bin/bash

# 步长 dt
dt=0.002

# 输出文件名
file1="temp.txt"
file2="energy.txt"
file3="msd.txt"

# 清空文件（如果文件已存在）
> "$file1"
> "$file2"
> "$file3"

# 逐行提取数据并输出到不同文件
awk '/^Step CPU/ {flag=1} /^Loop time/ {flag=0} flag' out > data.txt
awk -v dt="$dt" 'NR>1 {print $1 * dt, $4} ' data.txt > "$file1"
awk -v dt="$dt" 'NR>1 {print $1 * dt, $7} ' data.txt > "$file2"
awk -v dt="$dt" 'NR>1 {print $1 * dt, $10} ' data.txt > "$file3"

echo "Data extraction completed."
echo "Output files created: $file1, $file2, $file3"

