function [q0, q1, llr, c] = lfds_outliers(p0, p1, dx, eps)
% Get least favourable densities for two hypotheses under epsilon contamination
% uncertainty as a specail case of band uncertainty.
% 
% INPUT
%      p0:             nominal density under H0, 1xK vector
%      p1:             nominal density under H1, 1xK vector
%      dx:             grid size for numerical integraton
%      eps:            outlier ration, can be a scalar or a 2-tuple (eps0, eps1)
% 
% OUTPUT
%      q0, q1:         least favorable densities
%      llr:            log-likelihood ratio of q1 and q0, log(q1/q0)
%      c:              clipping constants c0, c1
%      nit:            number of iterations

% sanity checks
if ~is_nonnegative_vector(p0)
    error("'p0' must be a nonnegative array.");
end

if ~is_nonnegative_vector(p1)
    error("'p1' must be a nonnegative array.");
end

if length(p0) == length(p1)
    K = length(p0);
else
    error("'p0' and 'p1' need to be of the same size.")
end

% get outlier ratios
if is_nonnegative_scalar(eps)
    eps = [eps, eps];
end

if ~is_nonnegative_scalar(eps(1)) || eps(1) > 1.0
    error("outlier ratio 'eps0' must be between 0 and 1.")
end

if ~is_nonnegative_scalar(eps(2)) || eps(2) > 1.0
    error("outlier ratio 'eps1' must be between 0 and 1.")
end

% initialize bands corresponding to outlier model
p_min(1, :) = (1-eps(1)) * p0;
p_min(2, :) = (1-eps(2)) * p1;
p_max(1, :) = ones(1, K)/dx;
p_max(2, :) = ones(1, K)/dx;

% solve via density band algorithm
[q0, q1, llr, c] = lfds_density_band(p_min, p_max, dx, 0.0, p0, p1);