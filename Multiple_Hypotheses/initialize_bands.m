 function [verbose, skip_sanity_check, Q, eps, itmax, itmax_prox] = initialize_lfds_band(P_min, P_max, dw, nargin, argin)

% verbose
arg = 1;
if nargin >= arg && ~isempty(argin{arg})
    if islogical(argin{arg})
        verbose = argin{arg};
    else
        error('Verbose must be a logical.');
    end
else
    verbose = true;
end

% skip sanity checks
arg = 2;
if nargin >= arg && ~isempty(argin{arg})
    if islogical(argin{arg})
        skip_sanity_check = argin{arg};
    else
        error('Skip sanity checks must be a logical.');
    end
else
    skip_sanity_check = false;
end

% initial lfds Q
arg = 3;
if nargin >= arg && ~isempty(argin{arg})
    if all(argin{arg}(:) >= 0)
        Q = argin{arg};
    else
        error('Invalid specifiactions for intial lfds.');
    end
else
    [~, K] = size(P_min);
    a = (1/dw-sum(P_min,2))./(sum(P_max,2)-sum(P_min,2));
    Q = repmat((1-a), 1, K).*P_min + repmat(a, 1, K).*P_max;
end

% tolerance
arg = 4;
if nargin >= arg && ~isempty(argin{arg})
    if isscalar(argin{arg}) && argin{arg} > 0
        eps = argin{arg}(1);
    else
        error('Precision must be a positive sccalar.');
    end
else
    eps = 1e-6;
end

% maximum number of iterations
arg = 5;
if nargin >= arg && ~isempty(argin{arg})
    if isscalar(argin{arg}) && argin{arg} > 0
        itmax = argin{arg};
    else
        error('Maximum number of iterations must be positive.');
    end
else
    itmax = 100;
end

% maximum number of proximal iterations
arg = 6;
if nargin >= arg && ~isempty(argin{arg})
    if isscalar(argin{arg}) && argin{arg} > 0
        itmax = argin{arg};
    else
        error('Maximum number of iterations must be positive.');
    end
else
    itmax_prox = 50;
end

