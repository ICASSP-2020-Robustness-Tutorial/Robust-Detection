function value = is_nonnegative_scalar(x)

value = isscalar(x) && isreal(x) && x >= 0;