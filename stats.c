#include "stats.h"
#include "rng.h"

#include <math.h>
#include <stddef.h>

/* Internal helpers */
static inline double uniform01(void) {
    return random_double(); /* in [0,1) */
}

/* Box-Muller transform for standard normal N(0,1) with value caching */
static double randn(void) {
    static int has_spare = 0;
    static double spare = 0.0;
    if (has_spare) {
        has_spare = 0;
        return spare;
    }

    double u, v, s;
    do {
        u = 2.0 * uniform01() - 1.0; /* (-1,1) */
        v = 2.0 * uniform01() - 1.0; /* (-1,1) */
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);

    double mul = sqrt(-2.0 * log(s) / s);
    spare = v * mul;
    has_spare = 1;
    return u * mul;
}

/* Marsaglia-Tsang method for Gamma(alpha, rate=beta) */
double sample_gamma(double alpha, double beta) {
    if (!(alpha > 0.0) || !(beta > 0.0)) {
        return NAN;
    }

    if (alpha == 1.0) {
        /* Exponential with rate beta */
        double u;
        do { u = uniform01(); } while (u <= 0.0);
        return -log(u) / beta;
    }

    if (alpha < 1.0) {
        /* Boosting method: sample from alpha+1 and scale by U^{1/alpha} */
        double u;
        double x;
        do { u = uniform01(); } while (u <= 0.0);
        x = sample_gamma(alpha + 1.0, 1.0); /* rate=1 */
        return (x * pow(u, 1.0 / alpha)) / beta;
    }

    /* alpha >= 1 */
    const double d = alpha - 1.0 / 3.0;
    const double c = 1.0 / sqrt(9.0 * d);

    while (1) {
        double z = randn();
        double u = uniform01();

        double x = 1.0 + c * z;
        if (x <= 0.0) continue;
        double v = x * x * x; /* (1 + c z)^3 */

        /* Squeeze test */
        double z2 = z * z;
        double lhs = log(u);
        double rhs = 0.5 * z2 + d * (1.0 - v + log(v));
        if (lhs < -rhs) {
            return (d * v) / beta;
        }
    }
}

double log_gamma_pdf(double x, double alpha, double beta) {
    if (!(alpha > 0.0) || !(beta > 0.0) || !(x > 0.0)) {
        return -INFINITY;
    }
    return alpha * log(beta) - lgamma(alpha) + (alpha - 1.0) * log(x) - beta * x;
}

void sample_dirichlet(const double* alpha, int n, double* output) {
    if (!alpha || !output || n <= 0) return;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        double a = alpha[i];
        double g = sample_gamma(a, 1.0); /* rate=1 */
        output[i] = g;
        sum += g;
    }
    if (sum <= 0.0) {
        /* fallback: uniform if all zeros or invalid */
        double invn = 1.0 / (double)n;
        for (int i = 0; i < n; ++i) output[i] = invn;
        return;
    }
    double inv = 1.0 / sum;
    for (int i = 0; i < n; ++i) {
        output[i] *= inv;
    }
}
