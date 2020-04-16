# Least Favorable Densities for multiple Hypotheses

Calculate LFDs for a multiple-hypotheses test under density band and epsilon contamination uncertainty. Python, MATLAB and C implementations of the algorithm in 

[2] M. Fauß and A. M. Zoubir, "On the Minimization of Convex Functionals of Probability Distributions Under Band Constraints," in IEEE Transactions on Signal Processing, vol. 66, no. 6, pp. 1425-1437, March, 2018.

The algorithm solves a discrete version of the optimization problem

$min_{p_0, ..., p_1} \int f(p_0(x), ..., p_N(x)) dx$,

where $f \colon \mathbb{R}^{K+1} \to \mathbb{R}$ is a convex and homogeneous function, subject to the density band constraints

$ p_k'(x) \leq p_k(x) \leq p_k''(x)$.

See comments and the 'examples' files to get started.

The algorithm comes in two varieties, a *standard* version and a *proximal* version. The former is faster and should work well in most cases. If $f$ is not strictly convex, however, it might fail to converge. In this case the slower, but more robust proximal version should be used. 

While the Python and MATLAB implementations are useful to get a feeling for the problem formulation and for solving small-scale problems, the C implementation is substantially faster (~100x) and uses a less efficient, but more reliable method to calculate the residual errors. Hence, we highly recommend using the C implementation for large-scale problems or problems where $f$ can only be evaluated/approximated with numerical noise. 