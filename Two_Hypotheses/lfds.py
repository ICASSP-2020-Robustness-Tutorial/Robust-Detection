import sys
import os

import numpy as np
from numpy import linalg as la

sys.path.append(os.path.abspath("../../Helper_Functions"))
import robust_detection_helpers as hlp


def density_band(
    P_min,
    P_max,
    dx,
    verbose=False,
    alpha=0.01,
    Q_init=np.nan,
    order=np.inf,
    tol=1e-6,
    itmax=100
):
    """
    Get least favourable densities for two hypotheses under density band uncertainty
    For details see:

    M. Fauß and A. M. Zoubir, "Old Bands, New Tracks—Revisiting the Band Model for
    Robust Hypothesis Testing," in IEEE Transactions on Signal Processing,
    vol. 64, no. 22, pp. 5875-5886, 15 Nov.15, 2016.

    INPUT
        P_min:          2xK matrix specifying the lower bounds
        P_max:          2xK matrix specifying the upper bounds
        dx:             grid size for numerical integraton

    OPTIONAL INPUT
        verbose:        display progress, defaults to false
        alpha:          regularization parameter, defaults to 0.0
        Q_init:         2xK matrix, initial guess for q0, defaults to uniform density
        order:          vector norm used for convergence criterion, defaults to np.inf
        tol:            vector-norm tolerance of fixed-point, defaults to 1e-6
        itmax:          maximum number of iterations, defaults to 100

    OUTPUT
        q0, q1:         least favorable densities
        llr:            log-likelihood ratio of q1 and q0, log(q1/q0)
        c:              clipping constants c0, c1
        nit:            number of iterations
    """

    # sanity checks
    if P_min.shape == P_min.shape:
        N, K = np.shape(P_min)
    else:
        raise ValueError("'P_min' and 'P_max' must be of the same shape")

    for n in range(N):
        if not hlp.is_valid_density_band(P_min[n, :], P_max[n, :], dx):
            raise ValueError("Invalid density band.")

    if not alpha >= 0:
        raise ValueError("The parameter 'alpha' must be a nonegative scalar")

    # initialize lfds
    Q = hlp.set_densities(Q_init, P_min, P_max, dx)

    # rename for easier reference
    p0_min, p0_max = P_min[0, :], P_max[0, :]
    p1_min, p1_max = P_min[1, :], P_max[1, :]
    q0_new, q1_new = Q[0, :], Q[1, :]

    # initialize counters
    res = np.inf
    nit = 0
    c0_max, c1_max = 100.0, 100.0

    # disply prgoress, if verbose flag is set
    if verbose:
        print("\n")
        print("Iteration | Residual q0 | Residual q1 |    c0   |    c1   ")
        print("----------|-------------|-------------|---------|---------")

    # solve fixed-point equation iteratively
    while res > tol and nit < itmax:

        # assigne updated lfds
        q0 = q0_new
        q1 = q1_new

        # update q0
        def func0(c1):
            return (
                np.sum(np.minimum(p0_max, np.maximum(c1 * (alpha * q0 + q1), p0_min)))
                * dx
                - 1
            )

        c1, c1_max = hlp.get_pos_root(func0, c1_max)
        q0_new = np.minimum(p0_max, np.maximum(c1 * (alpha * q0 + q1), p0_min))

        # update q1 using q0_new (!)
        def func1(c0):
            return (
                np.sum(
                    np.minimum(p1_max, np.maximum(c0 * (alpha * q1 + q0_new), p1_min))
                )
                * dx
                - 1
            )

        c0, c0_max = hlp.get_pos_root(func1, c0_max)
        q1_new = np.minimum(p1_max, np.maximum(c0 * (alpha * q1 + q0_new), p1_min))

        # calculate residuals
        res0 = la.norm(q0_new - q0, order)
        res1 = la.norm(q1_new - q1, order)
        res = np.maximum(res0, res1)

        # count iterations
        nit += 1
        if verbose:
            print("%9d |  %.4e |  %.4e |  %.4f |  %.4f" % (nit, res0, res1, c0, c1))

    # check results
    if nit == itmax:
        print(f"{nit:d} iterations exeeded, possible numerical problem!")
    elif la.norm(q0 - q1, order) < tol:
        print("Overlapping Densities")

    # clipping constants
    c = c0, c1

    # log-likelihood ratio
    llr = np.log(q1 / q0)

    return np.vstack((q0, q1)), llr, c, nit


def outliers(P, dx, eps, verbose=False):
    """
    Get least favourable densities for two hypotheses under epsilon contamination
    uncertainty as a specail case of band uncertainty.

    INPUT
        P:              nominal densities, 2xK vector
        dx:             grid size for numerical integraton
        eps:            outlier ration, can be a scalar or a 2-tuple (eps0, eps1)

    OPTIONAL INPUT
        verbose:        display progress, defaults to false

    OUTPUT
        Q:              least favorable densities, 2xK vector
        llr:            log-likelihood ratio of q1 and q0, log(q1/q0)
        c:              clipping constants c0, c1
        nit:            number of iterations
    """

    # sanity checks
    if not np.all(P >= 0.0):
        raise ValueError("'P' must be a nonnegative matrix.")

    N, K = P.shape

    # get outlier ratios
    if np.isscalar(eps):
        eps = eps, eps

    for n in range(2):
        if not 0.0 <= eps[n] <= 1.0:
            raise ValueError("outlier ratio must be between 0 and 1.")

    # initialize bands corresponding to outlier model
    P_min = np.vstack(((1 - eps[0]) * P[0, :], (1 - eps[1]) * P[1, :]))
    P_max = np.ones_like(P) / dx

    # solve via density band algorithm
    Q, llr, c, _ = density_band(
        P_min, P_max, dx, verbose=verbose, Q_init=P
    )

    return Q, llr, c
