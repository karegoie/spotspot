#include <math.h>
#include <stddef.h>
#include <stdio.h>

#include "matrix.h"
#include "rng.h"
#include "stats.h"

/* Hyperparameters and proposal scales */
static const double SIGMA2 = 1.0;            /* observation variance for Gaussian likelihood */
static const double PROP_SD_X = 0.05;        /* proposal std for X entries */
static const double PROP_SD_Y = 0.05;        /* proposal std for Y entries */
static const double PRIOR_ALPHA_X = 1.0;     /* Gamma(shape=alpha, rate=beta) for X */
static const double PRIOR_BETA_X  = 1.0;
static const double PRIOR_ALPHA_Y = 1.0;     /* Gamma prior for Y */
static const double PRIOR_BETA_Y  = 1.0;

/* Local standard normal using Box-Muller with caching */
static double randn_local(void) {
    static int has_spare = 0;
    static double spare;
    if (has_spare) { has_spare = 0; return spare; }
    double u1, u2, r2;
    do {
        u1 = 2.0 * random_double() - 1.0;
        u2 = 2.0 * random_double() - 1.0;
        r2 = u1*u1 + u2*u2;
    } while (r2 >= 1.0 || r2 == 0.0);
    double f = sqrt(-2.0 * log(r2) / r2);
    spare = u2 * f; has_spare = 1;
    return u1 * f;
}

/* Compute initial product P = X * Y */
static Matrix* compute_product(const Matrix* X, const Matrix* Y) {
    return matrix_multiply(X, Y);
}

/* Log-likelihood contribution for a single element using Gaussian with variance SIGMA2 */
static inline double ll_element(double b, double mean) {
    double r = b - mean;
    return -0.5 * (r * r) / SIGMA2;
}

void update_parameters_fixed_dim(Matrix* X, Matrix* Y, const Matrix* B) {
    if (!X || !Y || !B) return;
    if (X->cols != Y->rows || B->rows != X->rows || B->cols != Y->cols) return;

    const int G = X->rows;     /* genes */
    const int C = X->cols;     /* cells */
    const int S = Y->cols;     /* spots */

    Matrix* P = compute_product(X, Y); /* current mean matrix */
    if (!P) return;

    /* Sweep X entries */
    for (int g = 0; g < G; ++g) {
        for (int c = 0; c < C; ++c) {
            double x_old = X->data[g][c];
            double x_prop = x_old + PROP_SD_X * randn_local();
            if (x_prop <= 0.0) {
                continue; /* gamma prior support: reject */
            }

            /* Prior terms */
            double lp_old = log_gamma_pdf(x_old, PRIOR_ALPHA_X, PRIOR_BETA_X);
            double lp_new = log_gamma_pdf(x_prop, PRIOR_ALPHA_X, PRIOR_BETA_X);

            /* Likelihood terms only row g is affected */
            double ll_old = 0.0, ll_new = 0.0;
            double dx = x_prop - x_old;
            for (int s = 0; s < S; ++s) {
                double ycs = Y->data[c][s];
                double m_old = P->data[g][s];
                double m_new = m_old + dx * ycs;
                double bgs = B->data[g][s];
                ll_old += ll_element(bgs, m_old);
                ll_new += ll_element(bgs, m_new);
            }

            double log_acc = (ll_new + lp_new) - (ll_old + lp_old);
            double u = random_double();
            if (log(u) < log_acc) {
                /* accept: update X and P */
                X->data[g][c] = x_prop;
                for (int s = 0; s < S; ++s) {
                    P->data[g][s] += dx * Y->data[c][s];
                }
            }
        }
    }

    /* Sweep Y entries */
    for (int c = 0; c < C; ++c) {
        for (int s = 0; s < S; ++s) {
            double y_old = Y->data[c][s];
            double y_prop = y_old + PROP_SD_Y * randn_local();
            if (y_prop <= 0.0) {
                continue; /* gamma prior support: reject */
            }

            double lp_old = log_gamma_pdf(y_old, PRIOR_ALPHA_Y, PRIOR_BETA_Y);
            double lp_new = log_gamma_pdf(y_prop, PRIOR_ALPHA_Y, PRIOR_BETA_Y);

            /* Likelihood: column s is affected across g */
            double ll_old = 0.0, ll_new = 0.0;
            double dy = y_prop - y_old;
            for (int g = 0; g < G; ++g) {
                double xgc = X->data[g][c];
                double m_old = P->data[g][s];
                double m_new = m_old + dy * xgc;
                double bgs = B->data[g][s];
                ll_old += ll_element(bgs, m_old);
                ll_new += ll_element(bgs, m_new);
            }

            double log_acc = (ll_new + lp_new) - (ll_old + lp_old);
            double u = random_double();
            if (log(u) < log_acc) {
                /* accept: update Y and P */
                Y->data[c][s] = y_prop;
                for (int g = 0; g < G; ++g) {
                    P->data[g][s] += dy * X->data[g][c];
                }
            }
        }
    }

    free_matrix(P);
}
