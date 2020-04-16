# Least Favorable Densities for Multiple Hypotheses

C implementation of the algorithm in 

[2] M. Fauß and A. M. Zoubir, "On the Minimization of Convex Functionals of Probability Distributions Under Band Constraints," in IEEE Transactions on Signal Processing, vol. 66, no. 6, pp. 1425-1437, March, 2018.

See comments and the example in '/example/cfpd_example.c' to get started.

### Build Requirements

C Math Library - \<math.h\>

GNU Scietific library (GSL), including development packages

CBLAS compatible BLAS library, such as OpenBLAS or ATLAS, including development packages

### Build

```
mkdir build
cd build
cmake ..
make
```
Run '/build/tests/cfpd_tests' to check if everything works correctly.