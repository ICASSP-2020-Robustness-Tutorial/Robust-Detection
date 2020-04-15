function [Q, llr, c, nit] = lfds_outliers(P, dx, eps, varargin)
% Get least favourable densities for two hypotheses under epsilon
% contamination uncertainty (outliers). For details see:
%
% M. Fauß and A. M. Zoubir, "Old Bands, New Tracks—Revisiting the Band Model for Robust Hypothesis
% Testing," in IEEE Transactions on Signal Processing, vol. 64, no. 22, pp. 5875-5886, 15 Nov.15, 2016.
%
% INPUT
%   P:              nominal densities, 2xK matrix
%   dx:             grid size for numerical integraton 
%   eps:            outlier ratio, can be a scalar or a vector
%
% varargin
% | {1}:            display progress, defaults to false
% | {2}:            regularization parameter, defaults to 0.0
% | {3}:            initial guess for Q, defaults to weighted sum of lower
%                   and upper bound  
% | {4}:            tolerance of fixed-point in terms of sup-norm, defaults to 1e-6
% | {5}:            maximum number of iterations, defaults to 100
% | {6}:            order of vector norm used for convergence criterion, defaults to Inf   
%
% OUTPUT
%   Q:              least favorable densities
%   llr:            log-likelihood ratio of q1 and q0, log(q1/q0)
%   c:              clipping constants c0, c1
%   nit:            number of iterations

% add path to helper functions
addpath ../Helper_Functions

[N, K] = size(P);

% get outlier ratios
if isscalar(eps)
    eps = eps*ones(N, 1);
end

if any(eps < 0) || any(eps > 1)
    error("outlier ratio must be between 0 and 1.")
end

% initialize bands corresponding to outlier model
P_min = (1-repmat(eps(:), 1, K)) .* P;
P_max = ones(N, K)/dx;

% solve via density band algorithm
[Q, llr, c] = lfds_density_band(P_min, P_max, dx, varargin{:});