#pragma once

#include "parameters.h"

typedef struct poly {
  int length, weight, alloc;
  uint16_t *index;
} *poly_t;

int compare_floats(const void *a, const void *b);
int compare_decrease(const void *a, const void *b, void *context);
poly_t poly_alloc(int length, int weight);
void poly_rand(poly_t h);
void poly_free(poly_t h);
double log_score_poly(poly_t h, double S[PARAM_N / 2][max_mult_y]);
poly_t poly_shift_union_soft(double *res, uint16_t *V, poly_t g, poly_t f,
                             int l, double S[PARAM_N / 2][max_mult_y]);
