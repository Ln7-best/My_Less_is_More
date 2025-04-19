#!/bin/bash

# 设置文件路径
CONFIG_FILE="/home/ln7/Octosketch/CPU/CM/config.h"
OUTPUT_FILE="throughput_results.txt"

# 清空输出文件
> "$OUTPUT_FILE"

# 循环修改 THREAD_NUM 的值
for THREAD_NUM in $(seq 4 4 32); do
    # 使用 sed 修改 config.h 中的 THREAD_NUM
    sed -i "s/#define THREAD_NUM .*/#define THREAD_NUM $THREAD_NUM/" "$CONFIG_FILE"
    
    # 编译项目
    make

    # 存储当前 THREAD_NUM 的 throughput 值
    total_throughput=0

    # 运行程序 10 次并抓取 throughput
    for i in $(seq 1 10); do
        OUTPUT=$(./CM /home/ln7/Octosketch/dataset/output_0.99_super_large.bin)
        
        # 从输出中提取 throughput
        throughput=$(echo "$OUTPUT" | grep "throughput" | awk '{print $2}')
        echo $throughput
        # 累加 throughput 用于计算平均值
        total_throughput=$(echo "$total_throughput + $throughput" | bc -l)
    done

    # 计算并输出平均 throughput
    avg_throughput=$(echo "$total_throughput / 10" | bc -l)
    echo "THREAD_NUM: $THREAD_NUM, avg throughput: $avg_throughput" >> "$OUTPUT_FILE"
done

echo "所有结果已写入 $OUTPUT_FILE"
