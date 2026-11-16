#define _GNU_SOURCE

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "poly.h"
#include "spectrum.h"

int compare_floats(const void *a, const void *b) {
  double fa = *(const double *)a;
  double fb = *(const double *)b;

  if (fa < fb)
    return 1;
  if (fa > fb)
    return -1;
  return 0;
}

int compare_decrease(const void *a, const void *b, void *context) {
  const double *scores = context;

  int idx1 = *(const int *)a;
  int idx2 = *(const int *)b;

  return (scores[idx1] < scores[idx2]) - (scores[idx1] > scores[idx2]);
}

poly_t poly_alloc(int length, int weight) {
  poly_t h = malloc(sizeof(struct poly));
  h->length = length;
  h->weight = weight;
  h->index = (uint16_t *)calloc(weight, sizeof(uint16_t));
  return h;
}

static int rand_int(int n) { return random() % n; }

static void rand_sparse(uint16_t *index, int length, int weight) {
  for (int i = weight - 1; i >= 0; --i) {
    int j = i + rand_int(length - i);
    for (int l = weight - 1; l > i; --l) {
      if (index[l] == j) {
        j = i;
        break;
      }
    }
    index[i] = j;
  }
}

void poly_rand(poly_t h) { rand_sparse(h->index, h->length, h->weight); }

void poly_free(poly_t h) {
  free(h->index);
  free(h);
}

/*
   Compute the score of a polynomial h using a sum of log
*/
double log_score_poly(poly_t h, double S[PARAM_N / 2][max_mult_y]) {
  int8_t S_h[PARAM_N / 2];
  spectrum(S_h, h->index, h->length, h->weight);

  double score = 0.0;

  for (int dist = 0; dist < PARAM_N / 2; dist++) {
    int mult_start = S_h[dist];
    double p_sum = 0.0;

    if (mult_start > max_mult_y)
      continue;

    for (int mult = mult_start; mult < max_mult_y; mult++) {
      p_sum += S[dist][mult];
    }

    score += log(p_sum + 1e-10);
  }

  return score;
}

/*
  Compute the extension of h (monomials sorted by their score).
*/
poly_t poly_shift_union_soft(double *ext_scores, uint16_t *ext_monomials,
                             poly_t g, poly_t f, int l,
                             double S[PARAM_N / 2][max_mult_y]) {
  int i, j, k, r = g->length;
  poly_t h = poly_alloc(r, g->weight + f->weight + 1);

  memcpy(h->index, g->index, g->weight * sizeof(g->index[0]));
  /* append coeffs of f, avoiding duplicates */
  for (i = 0, k = g->weight; i < f->weight; ++i) {
    h->index[k] = (f->index[i] + l) % r;
    for (j = 0; j < g->weight; ++j)
      if (h->index[j] == h->index[k])
        break;
    if (j == g->weight)
      ++k;
  }

  h->weight = k + 1;

  for (i = 0; i < r; i++) {
    for (j = 0; j < k; j++) {
      if (h->index[j] == i) {
        ext_scores[i] =
            10.; // positions already in h so surely in the extension
        break;
      }
    }

    if (j == k) {
      h->index[k] = i;
      ext_scores[i] = log_score_poly(h, S);
    }
  }

  int index[r];
  for (int i = 0; i < r; i++)
    index[i] = i;

  qsort_r(index, r, sizeof(int), compare_decrease, ext_scores);

  for (int i = 0; i < r; i++) {
    ext_monomials[i] = index[i];
  }

  qsort(ext_scores, r, sizeof(double), compare_floats);

  h->weight--;

  return h;
}