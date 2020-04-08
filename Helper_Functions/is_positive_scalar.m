function value = is_positive_scalar(x)

value = isscalar(x) && isreal(x) && x > 0;