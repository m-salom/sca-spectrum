#pragma once

#include "poly.h"

poly_t normalize(poly_t h);
int reconstruct_soft(double S_soft[PARAM_N / 2][max_mult_y], long int nb_dec,
                     double sigma, poly_t h_normalized);

static const double thresholds_110k[9] = {
    0, 0, 0, -15., -20., -30., -50., -70., -85.}; // -90 general
static const double thresholds_190k[9] = {0,    0,    0,    -10., -15.,
                                          -20., -30., -35., -40.}; //-50
