import numpy as np

from scipy import optimize


def is_valid_density_band(P_min, P_max, dx):
    return (
        dx > 0.0
        and all(np.sum(P_min, 1) * dx <= 1)
        and all(np.sum(P_max, 1) * dx >= 1)
        and np.all(P_min >= 0)
        and np.all(P_min <= P_max)
    )


def get_pos_root(func, x_max_guess):
    x_max = x_max_guess
    while True:
        try:
            x0 = optimize.brentq(func, 0.0, x_max)
            return x0, x_max
        except ValueError:
            x_max *= 2


def set_densities(Q_init, P_min, P_max, dx):
    N, K = P_min.shape
    Q = np.zeros((N, K))
    if not np.isscalar(Q_init):
        for n in range(N):
            if np.any(np.isnan(Q_init[n, :])):
                a = (1 / dx - np.sum(P_min[n, :])) / (
                    np.sum(P_max[n, :]) - np.sum(P_min[n, :])
                )
                Q[n, :] = (1 - a) * P_min[n, :] + a * P_max[n, :]
            else:
                if np.all(Q_init[n, :] >= 0.0):
                    Q[n, :] = Q_init[n, :]
                else:
                    raise ValueError(
                        f"User supplied initialization for q{n} is invalid."
                    )
    else:
        a = (1 / dx - np.sum(P_min, 1)) / (np.sum(P_max, 1) - np.sum(P_min, 1))
        Q = (1 - a.reshape(N, 1)) * P_min + a.reshape(N, 1) * P_max

    return Q
