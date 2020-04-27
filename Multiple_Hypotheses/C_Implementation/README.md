# Least Favorable Densities for Multiple Hypotheses

C implementation of the algorithm in 

[2] M. Fauß and A. M. Zoubir, "On the Minimization of Convex Functionals of Probability Distributions Under Band Constraints," in IEEE Transactions on Signal Processing, vol. 66, no. 6, pp. 1425-1437, March, 2018.



### Build Requirements

- C Math Library \<math.h\>

- GNU Scietific Library (GSL), including development packages

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

This builds a library 'liblfds.so' and installs it into '/usr/local/lib'. The corresponding header is installed in '/usr/local/include'. 

### Tests

Run '/build/tests/lfds_tests' to check if everything works correctly.

### Solve an Example Problem

Please see '/example/lfds_example.c' and the comments within this file to get started. It shows how to construct a problem, how to set the parameters of the algorithm, if desired, and how to solve it. The example can be used as a template to solve custom problems.

Running '/build/example/lfds_example' solves a small problem using both the regular and the proximal version of the algorithm. It also stores a text file 'lfds_example.dat' with the LFDs. In order to plot those to a PDF file using *gnuplot*, run

```
gnuplot lfds_example_plot.p
```