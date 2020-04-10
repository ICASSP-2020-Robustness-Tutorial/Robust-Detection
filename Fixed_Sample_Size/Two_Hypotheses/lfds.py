import sys
import os

import numpy as np
from numpy import linalg as la

sys.path.append(os.path.abspath("../../Helper_Functions"))
import robust_detection_helpers as hlp


def density_band(
    p_min,
    p_max,
    dx,
    alpha=0.01,
    q0_init=np.nan,
    q1_init=np.nan,
    order=np.inf,
    tol=1e-6,
    itmax=100,
):
    """
    Get least favourable densities for two hypotheses under density band uncertainty
    For details see:

    M. Fauß and A. M. Zoubir, "Old Bands, New Tracks—Revisiting the Band Model for
    Robust Hypothesis Testing," in IEEE Transactions on Signal Processing,
    vol. 64, no. 22, pp. 5875-5886, 15 Nov.15, 2016.

    INPUT
        p_min:          2xK vector specifying the lower bounds
        p_max:          2xK vector specifying the upper bounds
        dx:             grid size for numerical integraton

    OPTIONAL INPUT
        alpha           regularization parameter, defaults to 0.0
        q0_init         initial guess for q0, defaults to uniform density
        q1_init         initial guess for q1, defaults to uniform density
        order           vector norm used for convergence criterion, defaults to np.inf
        tol             vector-norm tolerance of fixed-point, defaults to 1e-6
        itmax:          maximum number of iterations, defaults to 100

    OUTPUT
        q0, q1:         least favorable densities
        llr:            log-likelihood ratio of q1 and q0, log(q1/q0)
        c:              clipping constants c0, c1
        nit:            number of iterations
    """

    # sanity checks
    if p_min.shape[0] == p_min.shape[0] == 2 and p_min.shape[1] == p_max.shape[1]:
        K = p_min.shape[1]
        p0_min = p_min[0, :]
        p0_max = p_max[0, :]
        p1_min = p_min[1, :]
        p1_max = p_max[1, :]
    else:
        raise ValueError("'p_min' and 'p_max' must be of size 2xK")

    if not hlp.is_nonnegative_scalar(alpha):
        raise ValueError("The parameter 'alpha' must be a nonegative scalar")

    if not hlp.is_valid_density_band(p0_min, p0_max, dx):
        raise ValueError("Invalid density band under H0.")

    if not hlp.is_valid_density_band(p1_min, p1_max, dx):
        raise ValueError("Invalid density band under H1.")

    # user defined initialization for q0
    if np.any(np.isnan(q0_init)):
        q0_new = np.ones(K)
    else:
        if np.all(q0_init >= 0.0) and q0_init.ndim == 1 and q0_init.size == K:
            q0_new = q0_init
        else:
            raise ValueError("User supplied initialization for q0 is invalid.")

    # user defined initialization for q1
    if np.any(np.isnan(q1_init)):
        q1_new = np.ones(K)
    else:
        if np.all(q1_init >= 0.0) and q1_init.ndim == 1 and q1_init.size == K:
            q1_new = q1_init
        else:
            raise ValueError("User supplied initialization for q1 is invalid.")

    # initialize counters
    dist = np.inf
    nit = 0
    c0_max, c1_max = 100.0, 100.0

    # solve fixed-point equation iteratively
    while dist > tol and nit < itmax:

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

        # calculate norm of difference
        dist = np.maximum(la.norm(q0_new - q0, order), la.norm(q1_new - q1, order))

        # count iterations
        nit += 1

    # check results
    if nit == itmax:
        print(f"{nit:d} iterations exeeded, possible numerical problem!")
    elif la.norm(q0 - q1, order) < tol:
        print("Overlapping Densities")

    # clipping constants
    c = c0, c1

    # log-likelihood ratio
    llr = np.log(q1 / q0)

    return q0, q1, llr, c, nit


def outliers(p0, p1, dx, eps):
    """
    Get least favourable densities for two hypotheses under epsilon contamination
    uncertainty as a specail case of band uncertainty.

    INPUT
        p0:             nominal density under H0, 1xK vector
        p1:             nominal density under H1, 1xK vector
        dx:             grid size for numerical integraton
        eps:            outlier ration, can be a scalar or a 2-tuple (eps0, eps1)

    OUTPUT
        q0, q1:         least favorable densities
        llr:            log-likelihood ratio of q1 and q0, log(q1/q0)
        c:              clipping constants c0, c1
        nit:            number of iterations
    """

    # sanity checks
    if not np.all(p0 >= 0.0) or p0.ndim != 1:
        raise ValueError("'p0' must be a nonnegative array.")

    if not np.all(p1 >= 0.0) or p1.ndim != 1:
        raise ValueError("'p0' must be a nonnegative array.")

    if p0.size == p1.size:
        K = p0.size
    else:
        raise ValueError("'p0' and 'p1' need to be of the same size.")

    # get outlier ratios
    if hlp.is_nonnegative_scalar(eps):
        eps = eps, eps

    if not hlp.is_nonnegative_scalar(eps[0]) or eps[0] > 1.0:
        raise ValueError("outlier ratio 'eps0' must be between 0 and 1.")

    if not hlp.is_nonnegative_scalar(eps[1]) or eps[1] > 1.0:
        raise ValueError("outlier ratio 'eps1' must be between 0 and 1.")

    # initialize bands corresponding to outlier model
    p_min = np.zeros((2, K))
    p_max = np.zeros((2, K))
    p_min[0, :] = (1 - eps[0]) * p0
    p_min[1, :] = (1 - eps[1]) * p1
    p_max[0, :] = np.ones_like(p0) / dx
    p_max[1, :] = np.ones_like(p1) / dx

    # solve via density band algorithm
    q0, q1, llr, c, _ = density_band(p_min, p_max, dx, q0_init=p0, q1_init=p1)

    return q0, q1, llr, c
