#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#include "model.h"
#include "rng.h"
#include "stats.h"

/* Helper: compute log-likelihood under Gaussian with variance sigma2 for given X,Y,B */
static double compute_log_likelihood(const Matrix* X, const Matrix* Y, const Matrix* B, double sigma2) {
    if (!X || !Y || !B) return -INFINITY;
    if (X->cols != Y->rows || B->rows != X->rows || B->cols != Y->cols) return -INFINITY;
    Matrix* P = matrix_multiply(X, Y);
    if (!P) return -INFINITY;
    double ll = 0.0;
    const double inv2s = 1.0 / (2.0 * sigma2);
    for (int i = 0; i < B->rows; ++i) {
        for (int j = 0; j < B->cols; ++j) {
            double r = B->data[i][j] - P->data[i][j];
            ll += -inv2s * (r * r);
        }
    }
    free_matrix(P);
    return ll;
}

/* Helper: log-prior for X and Y with Gamma priors */
static double compute_log_prior(const Matrix* X, const Matrix* Y, double ax, double bx, double ay, double by) {
    double lp = 0.0;
    for (int g = 0; g < X->rows; ++g) {
        for (int c = 0; c < X->cols; ++c) {
            double x = X->data[g][c];
            lp += log_gamma_pdf(x, ax, bx);
        }
    }
    for (int c = 0; c < Y->rows; ++c) {
        for (int s = 0; s < Y->cols; ++s) {
            double y = Y->data[c][s];
            lp += log_gamma_pdf(y, ay, by);
        }
    }
    return lp;
}

/* Beta(a,b) via Gamma draws */
static double sample_beta(double a, double b) {
    double g1 = sample_gamma(a, 1.0);
    double g2 = sample_gamma(b, 1.0);
    double s = g1 + g2;
    if (s <= 0.0) return 0.5;
    return g1 / s;
}

/* Resize helpers for X (GxC) and Y (CxS) adding/removing a column/row at index idx */
static int insert_column(Matrix* X, int idx) {
    int G = X->rows, C = X->cols;
    /* realloc row-wise arrays: we must reallocate each row to hold one more entry */
    for (int g = 0; g < G; ++g) {
        double* newrow = (double*)realloc(X->data[g], (size_t)(C + 1) * sizeof(double));
        if (!newrow) return 0;
        X->data[g] = newrow;
        for (int j = C; j > idx; --j) {
            X->data[g][j] = X->data[g][j - 1];
        }
        X->data[g][idx] = 0.0;
    }
    X->cols = C + 1;
    return 1;
}

static int remove_column(Matrix* X, int idx) {
    int G = X->rows, C = X->cols;
    for (int g = 0; g < G; ++g) {
        for (int j = idx; j < C - 1; ++j) {
            X->data[g][j] = X->data[g][j + 1];
        }
        double* newrow = (double*)realloc(X->data[g], (size_t)(C - 1) * sizeof(double));
        if (!newrow && C - 1 > 0) return 0;
        X->data[g] = newrow;
    }
    X->cols = C - 1;
    return 1;
}

static int insert_row(Matrix* Y, int idx) {
    int C = Y->rows, S = Y->cols;
    /* allocate new row pointers */
    double** newdata = (double**)malloc((size_t)(C + 1) * sizeof(double*));
    if (!newdata) return 0;
    for (int i = 0; i < idx; ++i) newdata[i] = Y->data[i];
    newdata[idx] = (double*)calloc((size_t)S, sizeof(double));
    if (!newdata[idx]) { free(newdata); return 0; }
    for (int i = idx; i < C; ++i) newdata[i + 1] = Y->data[i];
    free(Y->data);
    Y->data = newdata;
    Y->rows = C + 1;
    return 1;
}

static int remove_row(Matrix* Y, int idx) {
    int C = Y->rows, S = Y->cols;
    (void)S; /* unused */
    /* free the target row and compact */
    free(Y->data[idx]);
    for (int i = idx; i < C - 1; ++i) {
        Y->data[i] = Y->data[i + 1];
    }
    double** newdata = NULL;
    if (C - 1 > 0) {
        newdata = (double**)realloc(Y->data, (size_t)(C - 1) * sizeof(double*));
        if (!newdata) return 0;
    } else {
        free(Y->data);
    }
    Y->data = newdata;
    Y->rows = C - 1;
    return 1;
}

void perform_birth_step(ModelParameters* params) {
    if (!params || !params->X || !params->Y || !params->B) return;
    Matrix* X = params->X;
    Matrix* Y = params->Y;
    const Matrix* B = params->B;
    int G = X->rows, C = X->cols, S = Y->cols;

    /* Propose: new column for X sampled from Gamma prior, and split one existing Y row using Beta(u) */
    int c_star = (int)(random_double() * C);
    if (c_star < 0) c_star = 0;
    if (c_star >= C) c_star = C - 1;

    /* Sample new X column from Gamma prior */
    if (!insert_column(X, C)) return; /* append at end */
    for (int g = 0; g < G; ++g) {
        double xg = sample_gamma(params->alpha_x, params->beta_x);
        X->data[g][C] = (isnan(xg) || xg <= 0.0) ? 1e-6 : xg;
    }

    /* Split Y row c_star into two rows: scale by u and (1-u) */
    double u = sample_beta(1.0, 1.0); /* symmetric Beta(1,1) */
    if (!insert_row(Y, C)) { /* append at end */
        /* rollback X */
        remove_column(X, C);
        return;
    }
    /* Y rows were 0..C-1; we added new row at index C (last), now split c_star */
    for (int s = 0; s < S; ++s) {
        double val = Y->data[c_star][s];
        Y->data[c_star][s] = u * val;
        Y->data[C][s] = (1.0 - u) * val;
    }

    /* Compute log-acceptance ratio */
    double ll_new = compute_log_likelihood(X, Y, B, params->sigma2);
    double lp_new = compute_log_prior(X, Y, params->alpha_x, params->beta_x, params->alpha_y, params->beta_y);

    /* Build old model temporarily by merging back for evaluation */
    /* Copy Y and X for old */
    Matrix* X_old = create_matrix(G, C); /* old C */
    Matrix* Y_old = create_matrix(C, S);
    if (!X_old || !Y_old) { if (X_old) free_matrix(X_old); if (Y_old) free_matrix(Y_old); return; }
    for (int g = 0; g < G; ++g) for (int c = 0; c < C; ++c) X_old->data[g][c] = X->data[g][c];
    for (int c = 0; c < C; ++c) for (int s = 0; s < S; ++s) Y_old->data[c][s] = Y->data[c][s];
    /* Merge the two rows back in Y_old */
    for (int s = 0; s < S; ++s) {
        Y_old->data[c_star][s] = Y->data[c_star][s] + Y->data[C][s];
    }
    double ll_old = compute_log_likelihood(X_old, Y_old, B, params->sigma2);
    double lp_old = compute_log_prior(X_old, Y_old, params->alpha_x, params->beta_x, params->alpha_y, params->beta_y);

    free_matrix(X_old);
    free_matrix(Y_old);

    /* Proposal and Jacobian terms: Beta(1,1) split has Jacobian |det| = product over s of val (cancels when using proportional Y); we approximate as 0 (symmetric) */
    double log_q_ratio = 0.0; /* symmetric Beta(1,1) */
    double log_jacobian = 0.0; /* simple split without rescaling keeps sum */

    double log_acc = (ll_new + lp_new) - (ll_old + lp_old) + log_q_ratio + log_jacobian;
    double uacc = random_double();
    if (log(uacc) < log_acc) {
        /* accept: already in new state with C+1 */
        (void)0;
    } else {
        /* reject: rollback */
        /* merge the two Y rows back */
        for (int s = 0; s < S; ++s) {
            Y->data[c_star][s] += Y->data[C][s];
        }
        remove_row(Y, C);
        remove_column(X, C);
    }
}

void perform_death_step(ModelParameters* params) {
    if (!params || !params->X || !params->Y || !params->B) return;
    Matrix* X = params->X;
    Matrix* Y = params->Y;
    const Matrix* B = params->B;
    int G = X->rows, C = X->cols, S = Y->cols;
    if (C <= 1) return; /* cannot remove below 1 */

    /* Randomly pick a component to remove */
    int c_del = (int)(random_double() * C);
    if (c_del < 0) c_del = 0;
    if (c_del >= C) c_del = C - 1;
    int c_partner = (c_del == 0) ? 1 : 0; /* arbitrary partner to receive mass */

    /* Build proposed new state (C-1): merge Y rows and drop X column */
    Matrix* X_new = create_matrix(G, C - 1);
    Matrix* Y_new = create_matrix(C - 1, S);
    if (!X_new || !Y_new) { if (X_new) free_matrix(X_new); if (Y_new) free_matrix(Y_new); return; }
    /* Copy X columns except c_del */
    for (int g = 0; g < G; ++g) {
        int t = 0;
        for (int c = 0; c < C; ++c) {
            if (c == c_del) continue;
            X_new->data[g][t++] = X->data[g][c];
        }
    }
    /* Copy Y rows except c_del, merging into c_partner target index */
    int trow = 0;
    for (int c = 0; c < C; ++c) {
        if (c == c_del) continue;
        for (int s = 0; s < S; ++s) {
            double v = Y->data[c][s];
            if (c == c_partner) v += Y->data[c_del][s];
            Y_new->data[trow][s] = v;
        }
        trow++;
    }

    /* Compute log acc ratio comparing new (C-1) to old C */
    double ll_new = compute_log_likelihood(X_new, Y_new, B, params->sigma2);
    double lp_new = compute_log_prior(X_new, Y_new, params->alpha_x, params->beta_x, params->alpha_y, params->beta_y);
    double ll_old = compute_log_likelihood(X, Y, B, params->sigma2);
    double lp_old = compute_log_prior(X, Y, params->alpha_x, params->beta_x, params->alpha_y, params->beta_y);

    /* Symmetric proposal approx */
    double log_q_ratio = 0.0;
    double log_jacobian = 0.0;

    double log_acc = (ll_new + lp_new) - (ll_old + lp_old) + log_q_ratio + log_jacobian;
    double uacc = random_double();
    if (log(uacc) < log_acc) {
        /* accept: move to new state */
        /* Replace X and Y in-place via resize */
        /* First update Y: copy back */
        for (int i = 0; i < Y->rows; ++i) {
            if (Y->data[i]) free(Y->data[i]);
        }
        free(Y->data);
        Y->rows = Y_new->rows; Y->cols = Y_new->cols; Y->data = Y_new->data;
        /* prevent double free */ Y_new->data = NULL; free_matrix(Y_new);

        /* Now update X similarly */
        for (int g = 0; g < X->rows; ++g) free(X->data[g]);
        free(X->data);
        X->rows = X_new->rows; X->cols = X_new->cols; X->data = X_new->data;
        X_new->data = NULL; free_matrix(X_new);
    } else {
        /* reject */
        free_matrix(X_new);
        free_matrix(Y_new);
    }
}
