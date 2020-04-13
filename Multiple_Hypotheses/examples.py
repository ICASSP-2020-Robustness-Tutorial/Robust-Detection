import numpy as np
from scipy.stats import norm

import multi_lfds as mlfds

# Sample Space
dx = 0.01
x = np.arange(-10, 10, dx)
K = x.size

# Density Bands
N = 3
p1 = norm.pdf(x, -2, 3)
p2 = norm.pdf(x, 2, 2)
p3 = norm.pdf(x, 0.0, 1)

P = np.stack((p1, p2, p3))
Pmin = 0.8*P
Pmax = 1.2*P

# Choose objective function as weighted sum of KL devergences
a = np.array([0.5, 0.5])


def f(k, X):
    val = a[0]*np.log(X[2, k]/X[0, k])*X[2, k]
    val += a[1]*np.log(X[2, k]/X[1, k])*X[2, k]
    return val


# Derivative
def df_1(k, X):
    return -a[0]*X[2, k]/X[0, k]


def df_2(k, X):
    return -a[1]*X[2, k]/X[1, k]


def df_3(k, X):
    return 1 + a[0]*np.log(X[2, k]/X[0, k]) + a[1]*np.log(X[2, k]/X[1, k])


df = [df_1, df_2, df_3]

# get lfds
# Q, c, nit = multi_density_band(f, df, Pmin, Pmax, dx, verbose=True)

Q, c, nit = mlfds.outliers(f, df, P, dx, 0.1, verbose=True)
