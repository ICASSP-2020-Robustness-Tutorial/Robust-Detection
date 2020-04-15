# Robust Detection

This repository provides algorithms for calculating **least favorable densities (LFDs)** that can be used to implement minimax robust detectors. It is split into two main categories, namely LFDs for two hypotheses (binary) detection and LFDs for multi-hypotheses (m-ary) detection.

The uncertainty model used here is the **density band model**, which restricts feasible densities to be bounded from above and below. That is,

<img src="http://latex.codecogs.com/svg.latex?{p'(x) \leq p(x) \leq p''(x)}" border="0"/>

where *p* denotes the density function and *p'*, *p''* its lower and upper bound, respectively. The epsilon-contamination model (outlier model) is included as a special case. 

For the two-hypotheses case, the LFDs are independent of the cost funtion and the detection threshold. That is, they joinitly minimize all f-divergences. Hence, the LFDs in the two-hypotheses only depend on the uncertainty model. The repository provides Python and MATLAB implementations of an efficient algorithm for calculating these LFDs given upper and lower density bounds.

For the multi-hypotheses case this independence no longer holds. Hence, the LDSs are calculated for a particular cost function and threshold. The latter translate into an f-dissimilarity, which is then minimized under density band constraints. The repositoy provides Python, MATLAB and C implementations for an efficient algorithm for this task.


