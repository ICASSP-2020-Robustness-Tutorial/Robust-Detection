"""
Example for least favorable densities under band uncertainty
"""

import numpy as np
import matplotlib.pyplot as plt

import multi_lfds as mlfds

from scipy.stats import norm


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
P_min = 0.8*P
P_max = 1.2*P

# Choose objective function as weighted sum of KL divergences
a = np.array([0.5, 0.5])


# Objective function
def f(k, X):
    val = a[0]*np.log(X[2]/X[0])*X[2]
    val += a[1]*np.log(X[2]/X[1])*X[2]
    return val


# Derivative
def df(n, k, X):
    if n == 0:
        return -a[0]*X[2]/X[0]
    if n == 1:
        return -a[1]*X[2]/X[1]
    if n == 2:
        return 1 + a[0]*np.log(X[2]/X[0]) + a[1]*np.log(X[2]/X[1])
    else:
        raise ValueError(f"Invalid index n = {n}.")


"""
LFDs under density band uncertainty
"""

# denisty band model using regular algorithm
Q, I_val, c, nit = mlfds.density_band(f, df, P_min, P_max, dx, verbose=True)

# plot lfds
fig1, ax1 = plt.subplots()
for n in range(N):
    ax1.plot(x, Q[n, :], label=f"$q_{n}$")
ax1.set_title('Density band uncertainty - Regular')
legend = ax1.legend()

# denisty band model using proximal algorithm
Q, c, nit = mlfds.density_band_proximal(f, df, P_min, P_max, dx, verbose=True)

# plot lfds
fig2, ax2 = plt.subplots()
for n in range(N):
    ax2.plot(x, Q[n, :], label=f"$q_{n}$")
ax2.set_title('Density band uncertainty - Proximal')
legend = ax2.legend()


"""
LFDs under 10%  contamination (outliers)
"""

# # outlier model using regular algorithm
Q, c, nit = mlfds.outliers(f, df, P, dx, 0.1, verbose=True)

# plot lfds
fig3, ax3 = plt.subplots()
for n in range(N):
    ax3.plot(x, Q[n, :], label=f"$q_{n}$")
ax3.set_title('Outlier uncertainty - Regular')
legend = ax3.legend()


# # denisty band model using proximal algorithm
# # Warning, this is slow!
# Q, c, nit = mlfds.outliers_proximal(f, df, P, dx, 0.1, verbose=True)

# # plot lfds
# fig4, ax4 = plt.subplots()
# for n in range(N):
#     ax4.plot(x, Q[n, :], label=f"$q_{n}$")
# ax4.set_title('Outlier uncertainty - Proximal')
# legend = ax4.legend()
