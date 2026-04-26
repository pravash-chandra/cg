#!/bin/bash

# Comprehensive benchmark comparison script

echo "=========================================================================="
echo "CG SOLVER - COMPREHENSIVE BENCHMARK COMPARISON"
echo "=========================================================================="

# Compile
echo ""
echo "[1/3] Compiling Serial Version..."
gfortran -O3 cg_bench.f90 -o cg_bench
if [ $? -ne 0 ]; then
    echo "❌ Serial compilation failed"
    exit 1
fi
echo "✓ Serial compiled"

echo "[2/3] Compiling OpenMP Version..."
gfortran -O3 -fopenmp cg_parallel_bench.f90 -o cg_parallel_bench
if [ $? -ne 0 ]; then
    echo "❌ OpenMP compilation failed"
    exit 1
fi
echo "✓ OpenMP compiled"

# Run benchmarks
echo ""
echo "[3/3] Running Benchmarks..."
echo ""
echo "=========================================================================="
echo "SERIAL VERSION"
echo "=========================================================================="

SERIAL_OUTPUT=$(./cg_bench)
echo "$SERIAL_OUTPUT"

echo ""
echo "=========================================================================="
echo "OPENMP VERSION (Running with different thread counts)"
echo "=========================================================================="

echo ""
echo "--- 1 Thread ---"
export OMP_NUM_THREADS=1
OPENMP_1=$(./cg_parallel_bench)
echo "$OPENMP_1"

echo ""
echo "--- 2 Threads ---"
export OMP_NUM_THREADS=2
OPENMP_2=$(./cg_parallel_bench)
echo "$OPENMP_2"

echo ""
echo "--- 4 Threads ---"
export OMP_NUM_THREADS=4
OPENMP_4=$(./cg_parallel_bench)
echo "$OPENMP_4"

# Extract times and calculate speedup
echo ""
echo "=========================================================================="
echo "PERFORMANCE ANALYSIS"
echo "=========================================================================="
echo ""

# Parse serial times
SERIAL_TIMES=$(echo "$SERIAL_OUTPUT" | grep -E "^[[:space:]]*[0-9]+" | awk '{print $2}')

# Parse OpenMP 4-thread times
OPENMP_TIMES=$(echo "$OPENMP_4" | grep -E "^[[:space:]]*[0-9]+" | awk '{print $2}')

# Calculate speedup
echo "Speedup Comparison (Serial vs OpenMP 4 threads):"
echo "=========================================================================="
echo ""
echo "Matrix Size | Serial Time | OpenMP Time | Speedup | Efficiency"
echo "------------|-------------|-------------|---------|----------"

serial_array=($SERIAL_TIMES)
openmp_array=($OPENMP_TIMES)

for i in "${!serial_array[@]}"; do
    serial_time=${serial_array[$i]}
    openmp_time=${openmp_array[$i]}
    
    if [ ! -z "$serial_time" ] && [ ! -z "$openmp_time" ]; then
        # Get matrix size
        matrix_size=$(echo "$SERIAL_OUTPUT" | grep -E "^[[:space:]]*[0-9]+" | sed -n "$((i+1))p" | awk '{print $1}')
        
        # Calculate speedup
        speedup=$(echo "scale=3; $serial_time / $openmp_time" | bc)
        
        # Calculate efficiency (speedup / 4 threads * 100)
        efficiency=$(echo "scale=1; ($speedup / 4) * 100" | bc)
        
        printf "%10s | %11s | %11s | %7s | %8s\n" "$matrix_size" "$serial_time" "$openmp_time" "${speedup}x" "${efficiency}%"
    fi
done

# echo ""
# echo "=========================================================================="
# echo "SUMMARY"
# echo "=========================================================================="
# echo ""
# echo "✓ Serial implementation: Tests 1000 to 15000 matrix sizes"
# echo "✓ OpenMP implementation: Tests with 1, 2, and 4 threads"
# echo "✓ Matrix type: Dense symmetric positive definite"
# echo "✓ Convergence tolerance: 1.0e-6"
# echo ""
# echo "To run with custom thread count:"
# echo "  export OMP_NUM_THREADS=8"
# echo "  ./cg_openmp_benchmark"
# echo ""
# echo "=========================================================================="
