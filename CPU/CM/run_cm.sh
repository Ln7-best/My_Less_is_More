#!/bin/bash

CMD="./CM /mydata/data_generator_for_octo/fb_0.99_3.2e9.bin 1"
RUNS=10
SUM=0
COUNT=0

echo "Running command: $CMD for $RUNS times..."

for ((i=1; i<=RUNS; i++)); do
    OUTPUT=$($CMD 2>&1)  # 运行命令并捕获输出
    THROUGHPUT=$(echo "$OUTPUT" | grep -oP 'average throughput:\s*\K[0-9.]+')  # 提取吞吐量

    if [[ -n "$THROUGHPUT" ]]; then
        echo "Run $i: Throughput = $THROUGHPUT"
        SUM=$(echo "$SUM + $THROUGHPUT" | bc)
        ((COUNT++))
    else
        echo "Run $i: Failed to extract throughput!"
    fi
done

if [[ $COUNT -gt 0 ]]; then
    AVG=$(echo "scale=2; $SUM / $COUNT" | bc)
    echo "Average Throughput: $AVG"
else
    echo "Error: No valid throughput values extracted!"
fi