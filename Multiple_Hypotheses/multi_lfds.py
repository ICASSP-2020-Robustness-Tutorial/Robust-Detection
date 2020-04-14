import sys
import os

import numpy as np
from scipy import optimize

sys.path.append(os.path.abspath("../Helper_Functions"))
import robust_detection_helpers as hlp


def density_band(
    f, df, P_min, P_max, dx, verbose=False, Q_init=np.nan, tol=1e-6, itmax=100,
):
    """
    Get least favourable densities for multiple hypotheses under density band
    uncertainty. The least favorable densities minimize the f-dissimilarity defined
    by the function f, with partial derivatives df. For details on how to define
    f and df see 'example.py' and the paper

    M. Fauß and A. M. Zoubir, "On the Minimization of Convex Functionals of
    Probability Distributions Under Band Constraints," in IEEE Transactions on
    Signal Processing, vol. 66, no. 6, pp. 1425-1437, March, 2018.

    INPUT
        f:              convex function that defines an f-dissimilarity
        df:             partial derivative of f
        P_min:          2xK matrix specifying the lower bounds
        P_max:          2xK matrix specifying the upper bounds
        dx:             grid size for numerical integraton

    OPTIONAL INPUT
        verbose:        display progress, defaults to false
        Q_init:         2xK matrix, initial guess for q0, defaults to uniform density
        tol:            vector-norm tolerance of fixed-point, defaults to 1e-6
        itmax:          maximum number of iterations, defaults to 100

    OUTPUT
        Q:              least favorable densities, NxK matrix
        c:              clipping constants c0, c1
        nit:            number of iterations
    """

    # sanity checks
    if P_min.shape == P_min.shape:
        N, K = np.shape(P_min)
    else:
        raise ValueError("'P_min' and 'P_max' must be of the same shape")

    if not hlp.is_valid_density_band(P_min, P_max, dx):
        raise ValueError("Invalid density bands")

    # initialize lfds
    Q = hlp.set_densities(Q_init, P_min, P_max, dx)

    # initialize clipping constants
    c, c_min, c_max = get_c(N, K, df, P_min, P_max)

    # initialize residuals
    residuals = get_residuals(N, K, df, c, Q, P_min, P_max, dx)

    # initialize counter
    nit = 0

    # disply prgoress, if verbose flag is set
    if verbose:
        print("\n")
        print("Iteration | Residual Objective | Residual Densities")
        print("----------|--------------------|-------------------")
        print(
            "%9d |     %.4e     |     %.4e"
            % (nit, np.sum(residuals), np.max(np.abs(np.sum(Q, 1)) * dx - 1))
        )

    # iteratively solve optimality conditions
    while np.sum(residuals) > tol and nit < itmax:

        # select denisty with largest residual
        n = np.argmax(residuals)

        # define function for root search
        def func(c):
            return (
                np.sum(
                    np.minimum(
                        P_max[n, :],
                        np.maximum(
                            df_inv(df, n, N, K, Q, P_min[n, :], P_max[n, :], c),
                            P_min[n, :],
                        ),
                    )
                )
                * dx
                - 1
            )

        # find root
        c[n] = optimize.brentq(func, c_min[n], c_max[n])

        # update lfd
        Q[n, :] = np.minimum(
            P_max[n, :],
            np.maximum(
                df_inv(df, n, N, K, Q, P_min[n, :], P_max[n, :], c[n]), P_min[n, :]
            ),
        )

        # update residuals
        residuals = get_residuals(N, K, df, c, Q, P_min, P_max, dx)

        # display progress
        nit += 1
        if verbose:
            print(
                "%9d |     %.4e     |     %.4e"
                % (nit, np.sum(residuals), np.max(np.abs(np.sum(Q, 1)) * dx - 1))
            )

    return Q, c, nit


def density_band_proximal(
    f,
    df,
    P_min,
    P_max,
    dx,
    verbose=False,
    Q_init=np.nan,
    tol=1e-6,
    itmax=100,
    itmax_prox=50,
):
    """
    Get least favourable densities for multiple hypotheses under density band
    uncertainty. The proximal version of the algorithm is more robust and works
    for functions f that are not strictly convex. However, it is significantly
    slower and should only be used with problems that the regular version fails
    to solve.

    INPUT
        f:              convex function that defines an f-dissimilarity
        df:             partial derivative of f
        P_min:          2xK matrix specifying the lower bounds
        P_max:          2xK matrix specifying the upper bounds
        dx:             grid size for numerical integraton

    OPTIONAL INPUT
        verbose:        display progress, defaults to false
        Q_init:         2xK matrix, initial guess for q0, defaults to uniform density
        tol:            vector-norm tolerance of fixed-point, defaults to 1e-6
        itmax:          maximum number of iterations, defaults to 100
        itmax_prox:     maximum number of proximal iterations, defaults to 50

    OUTPUT
        Q:              least favorable densities, NxK matrix
        c:              clipping constants c0, c1
        nit:            number of iterations
    """

    # sanity checks
    if P_min.shape == P_min.shape:
        N, K = np.shape(P_min)
    else:
        raise ValueError("'P_min' and 'P_max' must be of the same shape")

    # initialize lfds
    Q = hlp.set_densities(Q_init, P_min, P_max, dx)

    # initialize clipping constants
    c, c_min, c_max = get_c(N, K, df, P_min, P_max)

    # initialize residuals
    residuals = get_residuals(N, K, df, c, Q, P_min, P_max, dx)

    # initialize counter
    nit_prox = 0

    # disply prgoress, if verbose flag is set
    if verbose:
        print("\n")
        print("Proximal Iteration | Residual Objective | Residual Densities")
        print("-------------------|--------------------|-------------------")
        print(
            "%18d |     %.4e     |     %.4e"
            % (nit_prox, np.sum(residuals), np.max(np.abs(np.sum(Q, 1)) * dx - 1))
        )

    # iteratively solve non-proximal problems with augmented df
    while np.sum(residuals) > tol and nit_prox < itmax_prox:

        # define proximal objective
        def df_prox(n, k, X):
            return df(n, k, X) + X[n] - Q[n, k]

        # solve proximal problem
        Q, c, _ = density_band(f, df_prox, P_min, P_max, dx, Q_init=Q)

        # update residuals
        residuals = get_residuals(N, K, df, c, Q, P_min, P_max, dx)

        # display progress
        nit_prox += 1
        if verbose:
            print(
                "%18d |     %.4e     |     %.4e"
                % (nit_prox, np.sum(residuals), np.max(np.abs(np.sum(Q, 1)) * dx - 1))
            )

    return Q, c, nit_prox


def outliers(
    f, df, P, dx, eps, verbose=False, tol=1e-6, itmax=100, itmax_prox=50, proximal=False
):
    """
    Get least favourable densities for multiple hypotheses under outlier
    uncertainty. This function calls the density band version internally

        INPUT
        f:              convex function that defines an f-dissimilarity
        df:             partial derivative of f
        P:              2xK matrix specifying the nominal densities
        dx:             grid size for numerical integraton
        eps:            contamination ratio, either scalar of of size N

    OPTIONAL INPUT
        verbose:        display progress, defaults to false
        tol:            vector-norm tolerance of fixed-point, defaults to 1e-6
        itmax:          maximum number of iterations, defaults to 100
        itmax_prox:     maximum number of proximal iterations, defaults to 50
        proximal:       use proximal version of the algorithm

    OUTPUT
        Q:              least favorable densities, NxK matrix
        c:              clipping constants c0, c1
        nit:            number of iterations
    """

    # sanity checks
    if not np.all(P >= 0.0):
        raise ValueError("'P' must be a nonnegative matrix.")

    N, K = P.shape

    # get outlier ratios
    if np.isscalar(eps):
        eps = np.ones(N) * eps

    if not np.all([0 <= eps[n] <= 1 for n in range(N)]):
        raise ValueError("outlier ratio 'eps' must be between 0 and 1.")

    # initialize bands corresponding to outlier model
    P_min = np.zeros((N, K))
    for n in range(N):
        P_min[n, :] = (1 - eps[n]) * P[n, :]
    P_max = 2 * P + 0.1

    # define output outside the loop
    Q = P
    c = np.ones(K)
    nit = 0

    while True:
        # solve density band model
        if proximal:
            Q, c, nit = density_band_proximal(
                f, df, P_min, P_max, dx, verbose=verbose, Q_init=Q, tol=tol, itmax=itmax
            )
        else:
            Q, c, nit = density_band(
                f, df, P_min, P_max, dx, verbose=verbose, Q_init=Q, tol=tol, itmax=itmax
            )

        # check if upper bound binds, increase if necessary
        if np.any(Q == P_max):
            P_max *= 2
            print("\nRe-running with adjusted upper bounds")
        else:
            break

    return Q, c, nit


def outliers_proximal(
    f, df, P, dx, eps, verbose=False, tol=1e-6, itmax=100, itmax_prox=50
):
    """
    Get least favourable densities for multiple hypotheses under outlier
    uncertainty using the proximal algorithm. This is merely a wrapper around
    'outlier' with the proximal parameter set to True.
    """

    return outliers(
        f,
        df,
        P,
        dx,
        eps,
        verbose=verbose,
        tol=tol,
        itmax=itmax,
        itmax_prox=itmax_prox,
        proximal=True,
    )


# invert partial derivative numerically, such that df(n, q) = c
def df_inv(df, n, N, K, Q, pmin, pmax, c):
    qn = np.zeros(K)
    for k in np.arange(K):
        q = Q[:, k].reshape(N, 1)

        def func(qnk):
            if n == 0:
                q_new = np.vstack((qnk, q[1:]))
            elif n == N - 1:
                q_new = np.vstack((q[:-1], qnk))
            else:
                q_new = np.vstack((q[:n], qnk))
                q_new = np.vstack((q_new, q[n + 1:]))
            return df(n, k, q_new) - c

        if func(pmin[k]) > 0.0:
            qn[k] = pmin[k]
        elif func(pmax[k]) < 0.0:
            qn[k] = pmax[k]
        else:
            qn[k] = optimize.brentq(func, pmin[k], pmax[k])

    return qn


# calculate residuals
def get_residuals(N, K, df, c, Q, P_min, P_max, dx):
    residuals = np.zeros(N)
    for n in np.arange(N):
        delta_c = df(n, np.arange(K), Q) - c[n]
        u = np.minimum(delta_c, 0)
        v = np.maximum(delta_c, 0)
        residuals[n] = ((Q[n, :] - P_max[n, :]) @ u + (Q[n, :] - P_min[n, :]) @ v) * dx

    return residuals


# get inital c and bounds
def get_c(N, K, df, P_min, P_max):
    c_min = np.array([np.min(df(n, np.arange(K), P_min)) for n in range(N)])
    c_max = np.array([np.max(df(n, np.arange(K), P_max)) for n in range(N)])

    return (c_min + c_max) / 2, c_min - 0.1, c_max + 0.1
