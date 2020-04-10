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
p0 = norm.pdf(w, -2, 2)
p1 = norm.pdf(w, 1, 1)

"""
LFDs under density band uncertainty
"""

# bands
p_min = np.zeros((2, w.size))
p_max = np.zeros((2, w.size))
p_min[0, :] = 0.8 * p0
p_min[1, :] = 0.8 * p1
p_max[0, :] = 1.2 * p0
p_max[1, :] = 1.2 * p1

# solve for LFDs
q0, q1, llr, c, nit = lfds.density_band(p_min, p_max, dx)

# plot lfds
fig1, ax1 = plt.subplots()
ax1.plot(w, q0, label='$q_0$')
ax1.plot(w, q1, label='$q_1$')
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
q0, q1, llr, c = lfds.outliers(p0, p1, dx, eps)

# plot lfds
fig3, ax3 = plt.subplots()
ax3.plot(w, q0, label='$q_0$')
ax3.plot(w, q1, label='$q_1$')
ax3.set_title('Epsilon contamination')
legend = ax3.legend()

# plot log-likelihood ratio
fig4, ax4 = plt.subplots()
ax4.plot(w, llr, label='log-likelihood ratio')
ax4.set_title('Epsilon contamination')
legend = ax4.legend()
