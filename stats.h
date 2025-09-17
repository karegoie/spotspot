#ifndef STATS_H
#define STATS_H

#ifdef __cplusplus
extern "C" {
#endif

/* Gamma(shape=alpha, rate=beta) */
double sample_gamma(double alpha, double beta);

double log_gamma_pdf(double x, double alpha, double beta);

/* Dirichlet(alpha[0..n-1]) -> output[0..n-1], sums to 1 */
void sample_dirichlet(const double* alpha, int n, double* output);

#ifdef __cplusplus
}
#endif

#endif /* STATS_H */
