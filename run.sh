#!/usr/bin/env bash
for i in {100..500}; do
    echo "====="
    ratio=$(bc -l <<< "scale=2; ${i}/100")
    cmake-build-debug/bvh_incremental model.ply ray_100000.bin 2 "${ratio}"
    python plot.py
done