function [Q, I, res, c, nit] = multi_lfds_density_band(f, df, f_param, P_min, P_max, dx, varargin)
% Get least favourable densities for multiple hypotheses under density band uncertainty. 
% The least favorable densities minimize the f-dissimilarity defined
% by the function f, with partial derivatives df. For details see
% 
% M. Fauß and A. M. Zoubir, "On the Minimization of Convex Functionals of
% Probability Distributions Under Band Constraints," in IEEE Transactions on
% Signal Processing, vol. 66, no. 6, pp. 1425-1437, March, 2018.
%
% INPUT
%
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
%
% OUTPUT
%     Q:                K x N dim. matrix. least favorable distributions
%     I:                Objective function evaluated at Q
%     residuals:        N dim. vector of residuals
%     c:                N dim. vector of derivative thresholds
%     nit:              Number of iterations
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

% default values
verbose = false;
Q_init = NaN;
tol = 1e-6;
itmax = 100;

% get dimensions
[N, K] = size(P_min);

% sanity checks
if ~is_valid_density_band(P_min, P_max, dx)
    error("Invalid density bands.");
end

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
if nargin >= 10 && ~isempty(varargin{4})
    if varargin{4} > 0
        itmax = varargin{4};
    else
        error('Maximum number of iterations must be a positive scalar.');
    end
end

% initialize lfds
Q = set_densities(Q_init, P_min, P_max, dx);

% get bounds for c
[c, c_min, c_max] = get_c(N, K, df, f_param, P_min, P_max);

% update residuals
res = get_residuals(N, K, df, Q, f_param, P_max, P_min, dx, c);

% display progress
nit = 0;
if verbose
    fprintf("\n");
    fprintf("Iteration | Residual Objective | Residual Densities\n");
    fprintf("----------|--------------------|-------------------\n");
    fprintf("%9d |     %.4e     |     %.4e\n", nit, sum(res), max(abs(sum(Q,2))*dx-1));
end


while sum(res) > tol && nit < itmax
      
    % select coordinate with largest residual
    [~, n] = max(res);
    
    % determine c(n)
    func = @(c) sum( get_idf(df, n, K, Q, f_param, P_min(n,:), P_max(n,:), c) )*dx - 1;
    c(n) = fzero(func, [c_min(n), c_max(n)]);
    Q(n,:) = get_idf(df, n, K, Q, f_param, P_min(n,:), P_max(n,:), c(n));
    
    % update residuals
    res = get_residuals(N, K, df, Q, f_param, P_max, P_min, dx, c);
    
    % display progress
    nit = nit+1;
    if verbose
        fprintf("%9d |     %.4e     |     %.4e\n", nit, sum(res), max(abs(sum(Q,2))*dx-1));
    end
    
end

% evaluate objective function
I = sum(f(1:K, Q, f_param))*dx;


% get inverse function of partial derivative of f
function qn = get_idf(df, n, K, Q, f_param, pmin, pmax, c)
qn = zeros(1,K);
for k=1:K
    func = @(qn) df(n, k, [Q(1:n-1,:); repmat(qn,1,K); Q(n+1:end,:)], f_param) - c;
    if func(pmin(k)) >= 0
        qn(k) = pmin(k);
    elseif func(pmax(k)) <= 0
        qn(k) = pmax(k);
    else
        qn(k) = fzero(func, [pmin(k), pmax(k)]);
    end
end
