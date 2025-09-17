#ifndef MODEL_H
#define MODEL_H

#include "matrix.h"

typedef struct {
    Matrix* X;           /* G x C */
    Matrix* Y;           /* C x S */
    const Matrix* B;     /* G x S */
    /* Hyperparameters */
    double sigma2;       /* observation variance */
    double alpha_x;      /* Gamma prior shape for X */
    double beta_x;       /* Gamma prior rate for X */
    double alpha_y;      /* Gamma prior shape for Y */
    double beta_y;       /* Gamma prior rate for Y */
} ModelParameters;

int init_model_from_B(ModelParameters* params, const Matrix* B, int C_init);

#endif /* MODEL_H */