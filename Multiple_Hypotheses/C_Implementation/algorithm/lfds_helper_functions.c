#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


double
lfds_vector_sum(const gsl_vector *v)
{
    double v_sum = 0;
    for (size_t i = 0; i < v->size; i++) {
        v_sum += gsl_vector_get(v, i);
    }

    return v_sum;
}


int
lfds_vector_finite(const gsl_vector *v)
{
    for (size_t n = 0; n < v->size; n++) {
        if (!gsl_finite(gsl_vector_get(v, n))) {
            return 0;
        }
    }

    return 1;
}


int
lfds_matrix_geq(const gsl_matrix *m1,
                const gsl_matrix *m2)
{
    if (m1->size1 != m2->size1 || m1->size2 != m2->size2) {
        GSL_ERROR(gsl_strerror(GSL_EBADLEN), GSL_EBADLEN);
    }

    for (size_t row = 0; row < m1->size1; row++) {
        for (size_t col = 0; col < m1->size2; col++) {
            if (gsl_matrix_get(m1, row, col) < gsl_matrix_get(m2, row, col)) {
                return 0;
            }
        }
    }

    return 1;
}


int
lfds_matrix_finite(const gsl_matrix *m)
{
    for (size_t row = 0; row < m->size1; row++) {
        for (size_t col = 0; col < m->size2; col++) {
            if (!gsl_finite(gsl_matrix_get(m, row, col))) {
                return 0;
            }
        }
    }

    return 1;
}
