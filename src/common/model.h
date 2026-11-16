#pragma once

#include "parameters.h"
#include <stdint.h>

extern double means[PARAM_N / 2];
extern int8_t S_y[PARAM_N / 2];

void shifting_values(double *syndromes, long int nb_dec, double a, double b);
void compute_model(double avg_dist_model[max_mult_y],
                   double std_dist_model[max_mult_y],
                   long double count_dist_model[max_mult_y], double noise,
                   long int nb_dec, double a, double b, double sigma0);
