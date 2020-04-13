import sys
import os

import numpy as np
from scipy import optimize
from scipy.optimize import fsolve

sys.path.append(os.path.abspath("../../Helper_Functions"))
import robust_detection_helpers as hlp


def density_band(
    f, df, P_min, P_max, dx, verbose=False, Q_init=np.nan, tol=1e-6, itmax=100,
):
    # sanity checks
    if P_min.shape == P_min.shape:
        N, K = np.shape(P_min)
    else:
        raise ValueError("'P_min' and 'P_max' must be of the same shape")

    if len(df) != N:
        raise ValueError("Length of 'df' must coincide with the number of hypotheses")

    for n in range(N):
        if not hlp.is_valid_density_band(P_min[n, :], P_max[n, :], dx):
            raise ValueError("Invalid density band under H%d." % n)

    # initialize lfds
    Q = np.zeros((N, K))

    # user defined initializations
    if not np.isscalar(Q_init):
        for n in range(Q_init.shape[0]):
            if np.any(np.isnan(Q_init[n, :])):
                Q[n, :] = P_min[n, :] + P_max[n, :]
                Q[n, :] = Q[n, :] / (np.sum(Q[n, :]) * dx)
            else:
                if np.all(Q_init[n, :] >= 0.0):
                    Q[n, :] = Q_init[n, :]
                else:
                    raise ValueError(
                        "User supplied initialization for q%d is invalid." % n
                    )
    else:
        Q = P_min + P_max
        Q = Q / (np.sum(Q, 1)[:, None] * dx)

    # initialize clipping constants
    c_min = np.array([np.min(df[n](np.arange(K), P_min)) for n in range(N)]) - 0.1
    c_max = np.array([np.max(df[n](np.arange(K), P_max)) for n in range(N)]) + 0.1
    c = (c_min + c_max) / 2

    # initialize residuals
    residuals = np.zeros(N)
    for n in np.arange(N):
        delta_c = df[n](np.arange(K), Q) - c[n]
        u = np.minimum(delta_c, 0)
        v = np.maximum(delta_c, 0)
        residuals[n] = ((Q[n, :] - P_max[n, :]) @ u + (Q[n, :] - P_min[n, :]) @ v) * dx

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
        # select hypothesis
        n = np.argmax(residuals)

        # define function for root search
        def func(c):
            return (
                np.sum(
                    np.minimum(
                        P_max[n, :],
                        np.maximum(
                            df_inv(df[n], n, N, K, Q, P_min[n, :], P_max[n, :], c),
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
                df_inv(df[n], n, N, K, Q, P_min[n, :], P_max[n, :], c[n]), P_min[n, :]
            ),
        )

        # update residuals
        for n in np.arange(N):
            delta_c = df[n](np.arange(K), Q) - c[n]
            u = np.minimum(delta_c, 0)
            v = np.maximum(delta_c, 0)
            residuals[n] = (
                (Q[n, :] - P_max[n, :]) @ u + (Q[n, :] - P_min[n, :]) @ v
            ) * dx

        # display progress
        nit += 1
        if verbose:
            print(
                "%9d |     %.4e     |     %.4e"
                % (nit, np.sum(residuals), np.max(np.abs(np.sum(Q, 1)) * dx - 1))
            )

    return Q, c, nit


# def density_band_proximal(
#     f, df, P_min, P_max, dx, verbose=False, Q_init=np.nan, tol=1e-6, itmax=100,
# ):
#     while True:
#         Q = P_min + P_max
#         Q = Q / (np.sum(Q, 1)[:, None] * dx)
        
#         for n in range(N):
#             def df_prox(k, X) = df[n](k, X) + X
    

def outliers(
    f, df, P, dx, eps, verbose=False, tol=1e-6, itmax=100,
):
    # sanity checks
    if not np.all(P >= 0.0):
        raise ValueError("'P' must be a nonnegative matrix.")

    N, K = P.shape

    # get outlier ratios
    if hlp.is_nonnegative_scalar(eps):
        eps = np.ones(N) * eps

    if not np.all([0 <= eps[n] <= 1 for n in range(N)]):
        raise ValueError("outlier ratio 'eps' must be between 0 and 1.")

    # initialize bands corresponding to outlier model
    P_min = np.zeros((N, K))
    P_max = np.zeros((N, K))
    for n in range(N):
        P_min[n, :] = (1 - eps[n]) * P[n, :]
        P_max = 2 * P + 0.1

    Q = P
    c = np.ones(K)
    nit = 0

    while True:
        Q, c, nit = density_band(
            f, df, P_min, P_max, dx, verbose=verbose, Q_init=Q, tol=tol, itmax=itmax
        )
        if np.any(Q == P_max):
            P_max *= 2
            print("\nRe-running with adjusted upper bounds")
        else:
            break

    return Q, c, nit


def lfds_band_df(f, df, Pmin, Pmax, dw, eps, verbose):

    N, K = np.shape(Pmin)
    U = np.zeros((N, K))
    V = np.zeros((N, K))

    Q = Pmin + Pmax
    Q = Q / (np.sum(Q, 1)[:, None] * dw)

    res = np.zeros(N) + eps
    cmin = np.zeros(N)
    cmax = np.zeros(N)

    for n in np.arange(N):
        cmin[n] = np.min(df[n](np.arange(K), Pmin)) - 0.1
        cmax[n] = np.max(df[n](np.arange(K), Pmax)) + 0.1

    c = (cmin + cmax) / 2

    for n in np.arange(N):
        U[n, :] = np.minimum(df[n](np.arange(K), Q) - c[n], 0)
        V[n, :] = np.maximum(df[n](np.arange(K), Q) - c[n], 0)
        res[n] = (
            np.dot(Q[n, :] - Pmax[n, :], U[n, :]) * dw
            + np.dot(Q[n, :] - Pmin[n, :], V[n, :]) * dw
        )

    it = 0
    itmax = 100

    if verbose:
        print("\n")
        print("Iteration | Residual Objective | Residual Densities")
        print("----------|--------------------|-------------------")
        print(
            "%9d |     %.4e     |     %.4e"
            % (it, np.sum(res), np.max(np.abs(np.sum(Q, 1)) * dw - 1))
        )

    while np.sum(res) > eps and it < itmax:

        it = it + 1

        n = np.argmax(res)

        def func(c):
            return (
                np.sum(
                    np.minimum(
                        Pmax[n, :],
                        np.maximum(
                            df_inv(df[n], n, N, K, Q, Pmin[n, :], Pmax[n, :], c),
                            Pmin[n, :],
                        ),
                    )
                )
                - 1 / dw
            )

        c[n] = fsolve(func, c[n])

        Q[n, :] = np.minimum(
            Pmax[n, :],
            np.maximum(
                df_inv(df[n], n, N, K, Q, Pmin[n, :], Pmax[n, :], c[n]), Pmin[n, :]
            ),
        )

        # update residuals
        for n in np.arange(N):
            U[n, :] = np.minimum(df[n](np.arange(K), Q) - c[n], 0)
            V[n, :] = np.maximum(df[n](np.arange(K), Q) - c[n], 0)
            res[n] = (
                np.dot(Q[n, :] - Pmax[n, :], U[n, :]) * dw
                + np.dot(Q[n, :] - Pmin[n, :], V[n, :]) * dw
            )

        # display progress
        if verbose:
            print(
                "%9d |     %.4e     |     %.4e"
                % (it, np.sum(res), np.max(np.abs(np.sum(Q, 1)) * dw - 1))
            )

    return Q


def df_inv(dfn, n, N, K, Q, pmin, pmax, c):

    qn = np.zeros(K)

    for k in np.arange(K):
        q = Q[:, k].reshape(3, 1)

        def func(qnk):
            if n == 0:
                q_new = np.vstack((qnk, q[1:]))
            elif n == N - 1:
                q_new = np.vstack((q[:-1], qnk))
            else:
                q_new = np.vstack((q[:n], qnk))
                q_new = np.vstack((q_new, q[n + 1 :]))
            return dfn(k, q_new) - c

        if func(pmin[k]) > 0.0:
            qn[k] = pmin[k]
        elif func(pmax[k]) < 0.0:
            qn[k] = pmax[k]
        else:
            qn[k] = optimize.brentq(func, pmin[k], pmax[k])  # fsolve(func, q[n])

    return qn
