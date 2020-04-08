function [q0, q1, llr, c, nit] = lfds_density_band(p_min, p_max, alpha, dx, varargin)
% Get least favourable densities for two hypotheses under density band uncertainty
% For details see:
%
% M. Fauß and A. M. Zoubir, "Old Bands, New Tracks—Revisiting the Band Model for Robust Hypothesis
% Testing," in IEEE Transactions on Signal Processing, vol. 64, no. 22, pp. 5875-5886, 15 Nov.15, 2016.
%
% INPUT
% p_min:            2xK vector specifying the lower bounds
% p_max:            2xK vector specifying the upper bounds
% dx:               grid size for numerical integraton 
%
% varargin
% | {1}:            initial guess for q0, defaults to uniform density
% | {2}:            initial guess for q1, defaults to uniform density
% | {3}:            vector norm used for convergence criterion, defaults to p = Inf     
% | {4}:            tolerance of fixed-point in terms of sup-norm, defaults to 1e-6
% | {5}:            maximum number of iterations, defaults to 100
%
% OUTPUT
% q0, q1:           least favorable densities
% llr:              log-likelihood ratio of q1 and q0, log(q1/q0)
% c:                clipping constants c0, c1
% nit:              number of iterations

% access to helper functions
addpath ../../Helper_Functions

% sanity checks
if size(p_min, 1) ~= 1 || size(p_max, 1) ~= 1
    if size(p_min, 2) == size(p_max, 2)
        K = size(p_min, 2);
        
        % use more explicit names
        p0_min = p_min(1,:); p0_max = p_max(1,:);
        p1_min = p_min(2,:); p1_max = p_max(2,:);
    else
        error("'p_min' and 'p_max' must be of size 2xK");
    end
end

if ~is_nonnegative_scalar(alpha)
    error("The parameter 'alpha' must be a nonegative scalar");
end

if ~is_valid_density_band(p0_min, p0_max, dx)
    error('Invalid density band under H0.');
end
if ~is_valid_density_band(p1_min, p1_max, dx)
    error('Invalid density band under H1.');
end

% default values
tol = 1e-6;
itmax = 100;
p = Inf;
c0 = 1; c1 = 1;
q0_new = ones(1, K)*dx;
q1_new = ones(1, K)*dx;

% user defined initialization for q0
if nargin >= 5 && ~isempty(varargin{1})
    arg = varargin{1};
    if is_nonnegative_vector(arg) && length(arg) == K
        q0_new = arg;
    else
        error('User supplied initialization for q0 is invalid.');
    end
end
    
% user defined initialization for q1
if nargin >= 6 && ~isempty(varargin{2})
    arg = varargin{2};
    if is_nonnegative_vector(arg) && length(arg) == K
        q1_new = arg;
    else
        error('User supplied initialization for q1 is invalid.');
    end
end

% user defined vector norm
if nargin >= 7 && ~isempty(varargin{3})
    if is_positive_scalar(varargin{3})
        p = varargin{3};
    else
        error('Vector norm parameter must be a positive scalar.');
    end
end

% user defined tolerance
if nargin >= 8 && ~isempty(varargin{4})
    if is_positive_scalar(varargin{4})
        tol = varargin{4};
    else
        error('Tolerance must be a positive scalar.');
    end
end

% user defined number of iterations
if nargin >= 9 && ~isempty(varargin{5})
    if is_positive_scalar(varargin{5})
        itmax = varargin{5};
    else
        error('Maximum number of iterations must be a positive scalar.');
    end
end

% initialize counters
dist = Inf;
nit = 0;

% solve fixed-point equation iteratively
while dist > tol && nit < itmax
    
    % assigne updated lfds
    q0 = q0_new;
    q1 = q1_new;
    
    % update q0
    func0 = @(c0) sum(min(p0_max, max(c0*(alpha*q0 + q1), p0_min))) - 1/dx;
    c0 = fzero(func0,c0);
    q0_new = min(p0_max, max(c0*(alpha*q0 + q1), p0_min));
        
    % update q1 using q0_new (!)
    func1 = @(c1) sum(min(p1_max, max(c1*(q0_new + alpha*q1), p1_min))) - 1/dx;
    c1 = fzero(func1, c1);
    q1_new = min(p1_max, max(c1*(q0_new + alpha*q1), p1_min));
    
    % calculate sup-norm
    dist = max( vecnorm(q0_new-q0, p), vecnorm(q1_new-q1, p) );
      
    % count iterations
    nit = nit+1;
       
end

% check results
if nit == itmax
    warning([int2str(itmax) ' iterations exeeded, possible numerical problem!']);
elseif vecnorm(q1-q0, p) < tol
    disp('   Overlapping densities!');
end

% scaling factors
c(1:2) = [c0 c1];

% log-likelihood ratio
llr = log(q1./q0);
