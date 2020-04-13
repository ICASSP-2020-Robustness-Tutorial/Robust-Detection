#include <gsl/gsl_math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


// Objective function
double
weighted_kl(const gsl_vector *x,
            const size_t      k,
            const void       *params)
{
    const gsl_vector *alpha = (gsl_vector*) params;
    size_t N = x->size;

    double result = 0;
    for (size_t n = 0; n < N - 1; n++)
        result += gsl_vector_get(alpha, n) * log(gsl_vector_get(x, N - 1) / gsl_vector_get(x, n));
    result *= gsl_vector_get(x, N - 1);

    return result;
}


// Objective function derivatives
double
weighted_kl_derivative(const size_t      n,
                       const gsl_vector *x,
                       const size_t      k,
                       const void       *params)
{
    const gsl_vector *alpha = (gsl_vector*) params;
    size_t N = x->size;

    switch(n) {
	case 0:
            return -gsl_vector_get(alpha, 0) * gsl_vector_get(x, N - 1) / gsl_vector_get(x, 0);
        case 1:
            return -gsl_vector_get(alpha, 1) * gsl_vector_get(x, N - 1) / gsl_vector_get(x, 1);
        case 2: {
            double result = 0;
            for (size_t n = 0; n < N-1; n++) {
                result += gsl_vector_get(alpha, n) * log(gsl_vector_get(x, N - 1) / gsl_vector_get(x, n));
            }
            result += 1;

            return result;
        }
    }
    return 0.0;
}
