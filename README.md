# Robust Detection

This repository provides algorithms for calculating **least favorable densities (LFDs)** that can be used to implement minimax robust detectors. It is split into two main categories, namely LFDs for two hypotheses (binary) detection and LFDs for multi-hypotheses (m-ary) detection. The corresponding algorithms are detailed in 

[1] M. Fauß and A. M. Zoubir, "Old Bands, New Tracks—Revisiting the Band Model for Robust Hypothesis Testing," in IEEE Transactions on Signal Processing, vol. 64, no. 22, pp. 5875-5886, Nov, 2016.

[2] M. Fauß and A. M. Zoubir, "On the Minimization of Convex Functionals of Probability Distributions Under Band Constraints," in IEEE Transactions on Signal Processing, vol. 66, no. 6, pp. 1425-1437, March, 2018.

The uncertainty model used here is the **density band model**, which restricts feasible densities to be bounded from above and below. That is,

*p'(x) ≤ p(x) ≤ p''(x)*,

where *p* denotes the density function and *p'*, *p''* its lower and upper bound, respectively. The epsilon-contamination model (outlier model) is included as a special case. 

For the two-hypotheses case, the LFDs are independent of the cost funtion and the detection threshold. That is, they joinitly minimize all f-divergences. Hence, the LFDs in the two-hypotheses only depend on the uncertainty model. The repository provides Python and MATLAB implementations of an efficient algorithm for calculating these LFDs given upper and lower density bounds.

For the multi-hypotheses case this independence no longer holds. Hence, the LFDs are calculated for a particular cost function and threshold. The latter translate into an f-dissimilarity, which is then minimized under density band constraints. The repositoy provides Python, MATLAB and C implementations for an efficient algorithm for this task.


