#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "matrix.h"
#include "rng.h"
#include "stats.h"
#include "mcmc.h"
#include "rjmcmc.h"
#include "io.h"
#include "model.h"

static void accumulate_mean(Matrix* acc, const Matrix* src, double w) {
    for (int i = 0; i < acc->rows; ++i) {
        for (int j = 0; j < acc->cols; ++j) {
            acc->data[i][j] += w * src->data[i][j];
        }
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s B.csv [iters=1000] [burn=200] [seed=1234] [C_init=2]\n", argv[0]);
        return 1;
    }
    const char* bpath = argv[1];
    int iters = (argc > 2) ? atoi(argv[2]) : 1000;
    int burn  = (argc > 3) ? atoi(argv[3]) : 200;
    unsigned long seed = (argc > 4) ? strtoul(argv[4], NULL, 10) : 1234UL;
    int C_init = (argc > 5) ? atoi(argv[5]) : 2;
    if (iters <= 0) iters = 1000;
    if (burn < 0) burn = 0;
    if (burn > iters) burn = iters / 2;
    if (C_init <= 0) C_init = 2;

    init_rng(seed);

    Matrix* B = load_csv_matrix(bpath);
    if (!B) {
        fprintf(stderr, "Failed to load CSV: %s\n", bpath);
        return 2;
    }

    ModelParameters params = {0};
    params.sigma2 = 1.0;
    params.alpha_x = 1.0; params.beta_x = 1.0;
    params.alpha_y = 1.0; params.beta_y = 1.0;
    if (!init_model_from_B(&params, B, C_init)) {
        fprintf(stderr, "Failed to init model\n");
        free_matrix(B);
        return 3;
    }

    /* Running means for X and Y after burn-in; note that C can change in RJMCMC, so we keep last dims only */
    Matrix* X_mean = create_matrix(params.X->rows, params.X->cols);
    Matrix* Y_mean = create_matrix(params.Y->rows, params.Y->cols);
    int mean_count = 0;

    for (int t = 0; t < iters; ++t) {
        /* RJ moves: choose birth or death with 0.5 */
        if (random_double() < 0.5) perform_birth_step(&params); else perform_death_step(&params);

        /* Fixed-dim updates */
        update_parameters_fixed_dim(params.X, params.Y, params.B);

        if (t >= burn) {
            /* If dims changed, reset means to current dims */
            if (X_mean->rows != params.X->rows || X_mean->cols != params.X->cols) {
                free_matrix(X_mean);
                X_mean = create_matrix(params.X->rows, params.X->cols);
            }
            if (Y_mean->rows != params.Y->rows || Y_mean->cols != params.Y->cols) {
                free_matrix(Y_mean);
                Y_mean = create_matrix(params.Y->rows, params.Y->cols);
            }
            accumulate_mean(X_mean, params.X, 1.0);
            accumulate_mean(Y_mean, params.Y, 1.0);
            mean_count++;
        }
    }

    if (mean_count > 0) {
        double inv = 1.0 / (double)mean_count;
        for (int i = 0; i < X_mean->rows; ++i)
            for (int j = 0; j < X_mean->cols; ++j)
                X_mean->data[i][j] *= inv;
        for (int i = 0; i < Y_mean->rows; ++i)
            for (int j = 0; j < Y_mean->cols; ++j)
                Y_mean->data[i][j] *= inv;
    }

    save_csv_matrix(X_mean, "X_mean.csv");
    save_csv_matrix(Y_mean, "Y_mean.csv");

    /* cleanup */
    free_matrix(X_mean);
    free_matrix(Y_mean);
    free_matrix(params.X);
    free_matrix(params.Y);
    free_matrix(B);
    return 0;
}
