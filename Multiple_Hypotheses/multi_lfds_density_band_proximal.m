function [Q, I, res, c, nit_prox] = multi_lfds_density_band_proximal(f, df, f_param, P_min, P_max, dx, varargin)
% Get least favourable densities for multiple hypotheses under density band
% uncertainty using the proximal version of the algorithm detailed in
% 
% M. Fauß and A. M. Zoubir, "On the Minimization of Convex Functionals of
% Probability Distributions Under Band Constraints," in IEEE Transactions on
% Signal Processing, vol. 66, no. 6, pp. 1425-1437, March, 2018.
%
% INPUT
%     f:                function handle of objective function
%     df:               N dim. cell array of function handles of the partial
%                       derivatives
%     f_param:          additional parameters of f and df
%     Pmin, Pmax:       matrices with density band specifications, each row corresponds 
%                       to one density bound
%     dw:               scalar, grid size parameter
%
% OPTIONAL INPUT
%     varargin
%     | {1}:            verbose {true, false}. Turn display of progress on or off.
%     |                 Defaults to true
%     | {2}:            initial lfds, K x N dim. matrix. Provide initial densities for
%     |                 the algorithm. 
%     | {3}:            tolerance, defaults to 1e-6.
%     | {4}:            maximum number of iterations, defaults to 100.
%     | {5}:            maximum number of proximal iterations, defaults to 50.
%
% OUTPUT
%     Q:                K x N dim. matrix. least favorable distributions
%     I:                Objective function evaluated at Q
%     residuals:        N dim. vector of residuals
%     c:                N dim. vector of derivative thresholds
%     nit_prox:         Number of proximal iterations
%
% REMARK
%
% The objective function must be of the form
% 
% f(k, X, f_param) 
%
% and must accept vectors of positive integers as the first argument and N
% x K dim. matrices as the second argument, for example, f(1:K, Pmin). Use
% X(n,k) to access X within the function
%
% The derivative function must be of the form
% 
% df(n, k, X, f_param) 
%
% where the first index indicates the element of X w.r.t. which the partial 
% derivative is taken. Again, use X(n,k) to access X.
%
% See the examples for how to define f and df.

% add path to helper functions
addpath ../Helper_Functions

% defaults
verbose = false;
Q_init = NaN;
tol = 1e-6;
itmax_prox = 50;

% verbosity
if nargin >= 7 && ~isempty(varargin{1})
    verbose = varargin{1};
end

% user defined Q
if nargin >= 8 && ~isempty(varargin{2})
    Q_init = varargin{2};
end

% user defined tolerance
if nargin >= 9 && ~isempty(varargin{3})
    if varargin{3} > 0
        tol = varargin{3};
    else
        error('Tolerance must be a positive scalar.');
    end
end

% user defined number of iterations
if nargin >= 11 && ~isempty(varargin{5})
    if varargin{5} > 0
        itmax_prox = varargin{5};
    else
        error('Maximum number of iterations must be a positive scalar.');
    end
end

% get dimensions
[N, K] = size(P_min);

% initialize lfds
Q = set_densities(Q_init, P_min, P_max, dx);

% get bounds for c
c = get_c(N, K, df, f_param, P_min, P_max);

% update residuals
res = get_residuals(N, K, df, Q, f_param, P_max, P_min, dx, c);

% silence inner iteration
varargin{1} = false;

% display progress
nit_prox = 0;
if verbose
    fprintf("\n");
    fprintf("Proximal Iteration | Residual Objective | Residual Densities\n");
    fprintf("-------------------|--------------------|-------------------\n");
    fprintf("%18d |     %.4e     |     %.4e\n", nit_prox, sum(res), max(abs(sum(Q,2))*dx-1));
end

% iteratively solve non-proximal problems with augmented df
while sum(res) > tol && nit_prox < itmax_prox

        % define proximal objective
        df_prox = @(n, k, X, f_param) df(n, k, X, f_param) + X(n,k) - Q(n, k);

        % use LFDs form previous iteration as starting point
        varargin{2} = Q;
        
        % solve proximal problem
        [Q, I, ~, c] = multi_lfds_density_band(f, df_prox, f_param, P_min, P_max, dx, varargin{:});

        % update residuals
        res = get_residuals(N, K, df, Q, f_param, P_max, P_min, dx, c);

        % display progress
        nit_prox = nit_prox+1;
        if verbose
            fprintf("%18d |     %.4e     |     %.4e\n", nit_prox, sum(res), max(abs(sum(Q,2))*dx-1));
        end

end
