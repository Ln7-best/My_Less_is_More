#!/bin/bash

set -e

echo "Compiling generator.cpp..."
g++ -O3 generator.cpp -o gen

if [ ! -f gen ]; then
    echo "Compilation failed: 'gen' not found."
    exit 1
fi

echo "Running gen..."
./gen simple_data.bin

echo "Done. Output written to simple_data.bin."