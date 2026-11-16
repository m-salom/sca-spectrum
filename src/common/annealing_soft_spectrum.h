#ifndef ANNEALING_SOFT_SPECTRUM_H
#define ANNEALING_SOFT_SPECTRUM_H

#include "parameters.h"
#include <stdint.h>

typedef struct {
  long long iterations;
  int restarts;
  int exploration_restarts;
  double start_temp;
  double end_temp;
  int init_top_k;
  int swap_trials;
  float swap_proba;
  uint64_t seed;
} anneal_params_t;

int annealing(double S_soft[PARAM_N / 2][max_mult_y], uint64_t *h, uint64_t *s,
              const anneal_params_t *params);

#endif