#!/bin/bash

# 设置文件路径和输出文件
CONFIG_FILE="/home/ln7/OctoSketch/CPU/CM/config.h"
OUTPUT_FILE="results.txt"

# 清空输出文件
> "$OUTPUT_FILE"

# 循环修改 LENGTH 值
for LENGTH in $(seq 65536 10240 1048576); do
    # 使用 sed 修改 config.h 中的 LENGTH
    sed -i "s/#define LENGTH .*/#define LENGTH $LENGTH/" "$CONFIG_FILE"
    
    # 编译程序
    make
    
    # 运行程序并抓取输出
    OUTPUT=$(./CM /home/ln7/dataset/output_0.99.bin)
    
    # 从输出中提取所需的值
    SKETCH_SIZE=$(echo "$OUTPUT" | grep "local sketch size" | awk '{print $4}')
    THROUGHPUT=$(echo "$OUTPUT" | grep "throughput" | awk '{print $2}')
    
    # 将结果写入文件
    echo "$SKETCH_SIZE $THROUGHPUT" >> "$OUTPUT_FILE"
done

echo "所有结果已写入 $OUTPUT_FILE"
