# Least Favorable Densities for Multiple Hypotheses

C implementation of the algorithm in 

[2] M. Fauß and A. M. Zoubir, "On the Minimization of Convex Functionals of Probability Distributions Under Band Constraints," in IEEE Transactions on Signal Processing, vol. 66, no. 6, pp. 1425-1437, March, 2018.

See comments and the example in '/example/lfds_example.c' to get started.

### Build Requirements

- C Math Library - \<math.h\>

- GNU Scietific library (GSL), including development packages

- OpenBLAS, including development packages

Other CBLAS compatible BLAS libraries can be used, but the make files have to be modified accordingly. 

### Build

```
mkdir build
cd build
cmake ..
make
sudo make install
```
This builds a library 'liblfds.so' and installs it into '/usr/local/lib'. The corresponding header is installed in '/usr/local/include'. Run '/build/tests/lfds_tests' to check if everything works correctly or solve the example problem by running '/build/example/lfds_example'.