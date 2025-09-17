#ifndef RNG_H
#define RNG_H

#ifdef __cplusplus
extern "C" {
#endif

/* Initialize the RNG with a seed */
void init_rng(unsigned long seed);

/* Return a double in [0.0, 1.0) */
double random_double(void);

#ifdef __cplusplus
}
#endif

#endif /* RNG_H */
