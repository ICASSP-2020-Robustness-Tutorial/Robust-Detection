function value = is_nonnegative_vector(x)

value = isreal(x) && isvector(x) && all(x >= 0);