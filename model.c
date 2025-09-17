#include "model.h"
#include "stats.h"
#include "rng.h"
#include <stdlib.h>
#include <math.h>

/* Initialize model parameters given B's dimensions: choose initial C (e.g., 2) and sample X,Y from priors */
int init_model_from_B(ModelParameters* params, const Matrix* B, int C_init) {
    if (!params || !B || C_init <= 0) return 0;
    int G = B->rows;
    int S = B->cols;
    params->X = create_matrix(G, C_init);
    params->Y = create_matrix(C_init, S);
    params->B = B;
    if (!params->X || !params->Y) {
        if (params->X) free_matrix(params->X);
        if (params->Y) free_matrix(params->Y);
        return 0;
    }
    if (params->sigma2 <= 0.0) params->sigma2 = 1.0;
    if (params->alpha_x <= 0.0) params->alpha_x = 1.0;
    if (params->beta_x  <= 0.0) params->beta_x  = 1.0;
    if (params->alpha_y <= 0.0) params->alpha_y = 1.0;
    if (params->beta_y  <= 0.0) params->beta_y  = 1.0;

    /* Sample X and Y from priors */
    for (int g = 0; g < G; ++g) {
        for (int c = 0; c < C_init; ++c) {
            double x = sample_gamma(params->alpha_x, params->beta_x);
            params->X->data[g][c] = (x <= 0.0 || isnan(x)) ? 1e-3 : x;
        }
    }
    for (int c = 0; c < C_init; ++c) {
        for (int s = 0; s < S; ++s) {
            double y = sample_gamma(params->alpha_y, params->beta_y);
            params->Y->data[c][s] = (y <= 0.0 || isnan(y)) ? 1e-3 : y;
        }
    }
    return 1;
}
