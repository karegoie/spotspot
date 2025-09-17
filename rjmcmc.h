#ifndef RJMCMC_H
#define RJMCMC_H

#include "model.h"

void perform_birth_step(ModelParameters* params);
void perform_death_step(ModelParameters* params);

#endif /* RJMCMC_H */