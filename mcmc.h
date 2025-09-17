#ifndef MCMC_H
#define MCMC_H

#include "matrix.h"

void update_parameters_fixed_dim(Matrix* X, Matrix* Y, const Matrix* B);

#endif /* MCMC_H */