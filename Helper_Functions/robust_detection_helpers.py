import numpy as np

from scipy import optimize


def is_nonnegative_scalar(x):
    return np.isscalar(x) and np.isreal(x) and x >= 0


def is_positive_scalar(x):
    return np.isscalar(x) and np.isreal(x) and x > 0


def is_valid_density_band(p_min, p_max, dx):
    return (
        is_positive_scalar(dx)
        and np.all(np.isreal(p_min))
        and p_min.ndim == 1
        and np.all(np.isreal(p_max))
        and p_max.ndim == 1
        and p_min.size == p_max.size
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
