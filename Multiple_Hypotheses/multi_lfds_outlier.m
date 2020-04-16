function [Q, I, residuals, c, nit] = multi_lfds_outlier(f, df, f_param, P, dx, eps, varargin)
% Get least favourable densities for multiple hypotheses under outlier
% uncertainty. This function calls the density band version internally
% 
% INPUT
%     f:            convex function that defines an f-dissimilarity
%     df:           partial derivative of f
%     f_param:      dditional parameters of f and df
%     P:            2xK matrix specifying the nominal densities
%     dw:           grid size for numerical integraton
%     eps:          contamination ratio, either scalar of of size N
% 
% OPTIONAL INPUT
%     varargin
%     | {1}:        verbose {true, false}. Turn display of progress on or off.
%     |             Defaults to true
%     | {2}:        initial lfds, K x N dim. matrix. Provide initial densities for
%     |             the algorithm. 
%     | {3}:        tolerance, defaults to 1e-6.
%     | {4}:        maximum number of iterations, defaults to 100.
% 
% OUTPUT
%     Q:            K x N dim. matrix. least favorable distributions
%     I:            Objective function evaluated at Q
%     residuals:    N dim. vector of residuals
%     c:            N dim. vector of derivative thresholds
%     nit:          Number of iterations

% add path to helper functions
addpath ../Helper_Functions

% get dimensions
[N, K] = size(P);

% get outlier ratios
if isscalar(eps)
    eps = ones(1, N) * eps;
end

if any(eps < 0) || any(eps > 1)
    error('Outlier ratio must be between 0 and 1.');
end

% initialize bands corresponding to outlier model
P_min = repmat(1-eps(:), 1, K) .* P;
P_max = 2*P + 0.1;

while true
    
    % solve density band model
    [Q, I, residuals, c, nit] = multi_lfds_density_band(f, df, f_param, P_min, P_max, dx, varargin{:});

    % check if upper bound binds, increase if necessary
    if any(Q == P_max)
        P_max = 2*P_max;
        print("\nRe-running with adjusted upper bounds")
    else
        break
    end
end
