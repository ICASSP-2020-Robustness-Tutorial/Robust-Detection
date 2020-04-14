import numpy as np

from scipy import optimize


def is_valid_density_band(p_min, p_max, dx):
    return (
        dx > 0.0
        and np.sum(p_min) * dx <= 1
        and np.sum(p_max) * dx >= 1
        and np.all(p_min >= 0)
        and np.all(p_min <= p_max)
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
    if len(P_min.shape) == 1:
        if np.any(np.isnan(Q_init)):
            a = (1 / dx - np.sum(P_min)) / (np.sum(P_max) - np.sum(P_min))
            Q = (1 - a) * P_min + a * P_max
        else:
            if np.all(Q_init >= 0.0):
                Q = Q_init
            else:
                raise ValueError("User supplied initialization for q is invalid.")
    else:
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
