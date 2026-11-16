#pragma once

#include <stdbool.h>
#include <stdint.h>

#include "parameters.h"

int distance(int i, int j, int r);
void spectrum(int8_t *S, uint16_t *index, int r, int d);
void get_soft_spectrum(double S_soft[PARAM_N / 2][max_mult_y],
                       double *syndromes, long int nb_dec, double *means,
                       double sigma, double a, double b, double sigma0);