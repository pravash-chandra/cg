= Verification of Algorithm with standard implementation
using Julia programming language created matrix  A of size 5 by 5 
with diagonal entries 2,4,6,8,10
and b vector of size 5 with entries 1 for verifications of algorithm 
and solution is obtained by operator A\b

using LinearAlgebra
A = diagm(2:2:10)
b = ones(5)
x = A\b

== Serial implementations cg.f90 
file name : cg_openmmp.f90 

output :
 Converged at iteration:           5
 Solution vector x:
    0.500000    0.250000    0.166667    0.125000    0.100000
 
 Matrix A:
    2.000000    0.000000    0.000000    0.000000    0.000000
    0.000000    4.000000    0.000000    0.000000    0.000000
    0.000000    0.000000    6.000000    0.000000    0.000000
    0.000000    0.000000    0.000000    8.000000    0.000000
    0.000000    0.000000    0.000000    0.000000   10.000000


== openmp implementations 
file name : cg_openmmp.f90 
 Converged at iteration:           5
 Solution vector x:
    0.500000    0.250000    0.166667    0.125000    0.100000
 
 Matrix A:
    2.000000    0.000000    0.000000    0.000000    0.000000
    0.000000    4.000000    0.000000    0.000000    0.000000
    0.000000    0.000000    6.000000    0.000000    0.000000
    0.000000    0.000000    0.000000    8.000000    0.000000
    0.000000    0.000000    0.000000    0.000000   10.000000

= Speed Comparison serial and openmp

==========================================
CG Solver Benchmarking Suite
==========================================

Step 1: Compiling sequential version...
✓ Sequential compilation successful

Step 2: Compiling OpenMP parallel version...
✓ Parallel compilation successful

==========================================
Running Benchmarks
==========================================

Running sequential benchmark (this may take several minutes)...

 ====== CG Solver Benchmarking (Sequential) ======
 
    Matrix Size     Time (sec)
 -------------------------------------------
 Converged at iteration:           3
           5000       0.162067
 Converged at iteration:           3
          10000       0.534784
 Converged at iteration:           3
          15000       1.281737
 Converged at iteration:           3
          20000       8.224290
 Converged at iteration:           3
          25000       6.121648
 
 Benchmark completed!


Running parallel benchmark (this may take several minutes)...

 ====== CG Solver Benchmarking (OpenMP Parallel) ======
 
    Matrix Size     Time (sec)
 -------------------------------------------
 Converged at iteration:           3
           5000       0.154804
 Converged at iteration:           3
          10000       0.595266
 Converged at iteration:           3
          15000       1.516473
 Converged at iteration:           3
          20000       5.110471
 Converged at iteration:           3
          25000       5.073851
 
 Benchmark completed!

==========================================
Benchmark Results Comparison
==========================================

Matrix Size | Sequential (s) | Parallel (s) | Speedup
---------------------------------------------------

==========================================
Performance Analysis
==========================================

Note: Speedup > 1.0 means parallel is faster
Speedup = Sequential Time / Parallel Time

Test Environment:
- Compiler: GNU Fortran (GCC) 15.2.1 20260209
- OpenMP Threads: 4
- Optimization Level: -O3

==========================================
Benchmark Complete!
==========================================
