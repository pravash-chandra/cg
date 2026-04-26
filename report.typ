= Verification of Algorithm with standard implementation
using Julia programming language created matrix  A of size 5 by 5 
with diagonal entries 2,4,6,8,10
and b vector of size 5 with entries 1 for verifications of algorithm 
and solution is obtained by operator A\b

using LinearAlgebra
A = diagm(2:2:10)
b = ones(5)
x = A\b

 Converged at iteration:           5
 Solution vector x:
    0.500000    0.250000    0.166667    0.125000    0.100000
 
 Matrix A:
    2.000000    0.000000    0.000000    0.000000    0.000000
    0.000000    4.000000    0.000000    0.000000    0.000000
    0.000000    0.000000    6.000000    0.000000    0.000000
    0.000000    0.000000    0.000000    8.000000    0.000000
    0.000000    0.000000    0.000000    0.000000   10.000000
