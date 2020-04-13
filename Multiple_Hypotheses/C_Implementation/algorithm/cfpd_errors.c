#include "cfpd_errors.h"


const char *
cfpd_strerror(const int cfpd_errno)
{
    switch (cfpd_errno) {
    case CFPD_CONTINUE:
        return "continue";
    case CFPD_SOLVED:
        return "solved";
    case CFPD_MAXITER:
        return "exceeded max number of iterations";
    case CFPD_NOPROG:
        return "iteration is not making progress towards solution";
    case CFPD_INVALID_F:
        return "no objective function specified";
    case CFPD_INVALID_BANDS:
        return "invalid band specifications";
    case CFPD_INVALID_P:
        return "invalid densities";
    case CFPD_INVALID_MU:
        return "invalid mu";
    case CFPD_INVALID_TOLERANCES:
        return "invalid tolerances";
    case CFPD_INVALID_C:
        return "invalid c";
    case CFPD_NONCONVEX:
        return "nonconvex objective function";
    case CFPD_FAILURE:
        return "something went wrong, this should not have happend";
    default:
        return "unknown error code" ;
    }
}
