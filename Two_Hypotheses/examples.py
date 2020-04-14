"""
Example for least favorable densities under band uncertainty
"""

import numpy as np
import matplotlib.pyplot as plt

import lfds

from scipy.stats import norm

# support
dx = 0.01
w = np.arange(-8, 8, dx)

# nominal densities
P = np.vstack((norm.pdf(w, -2, 2), norm.pdf(w, 1, 1)))

"""
LFDs under density band uncertainty
"""

# bands
P_min, P_max = 0.8 * P, 1.2 * P

# solve for LFDs
Q, llr, c, nit = lfds.density_band(P_min, P_max, dx, verbose=True)

# plot lfds
fig1, ax1 = plt.subplots()
ax1.plot(w, Q[0, :], label='$q_0$')
ax1.plot(w, Q[1, :], label='$q_1$')
ax1.set_title('Density band uncertainty')
legend = ax1.legend()

# plot log-likelihood ratio
fig2, ax2 = plt.subplots()
ax2.plot(w, llr, label='log-likelihood ratio')
ax2.set_title('Density band uncertainty')
legend = ax2.legend()

"""
LFDs under 10% and 15% contamination (outliers)
"""

eps = 0.1, 0.15

# solve for LFDs
Q, llr, c = lfds.outliers(P, dx, eps, verbose=True)

# plot lfds
fig3, ax3 = plt.subplots()
ax3.plot(w, Q[0, :], label='$q_0$')
ax3.plot(w, Q[1, :], label='$q_1$')
ax3.set_title('Epsilon contamination')
legend = ax3.legend()

# plot log-likelihood ratio
fig4, ax4 = plt.subplots()
ax4.plot(w, llr, label='log-likelihood ratio')
ax4.set_title('Epsilon contamination')
legend = ax4.legend()
