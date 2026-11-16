#define _GNU_SOURCE
#include <omp.h>
#include <signal.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#include "guess_key.h"
#include "model.h"
#include "poly.h"
#include "spectrum.h"

static double get_time() {
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + t.tv_usec * 1e-6;
}

static int compare(const void *a, const void *b, void *context) {
  const double *scores = context;

  int i = *(const int *)a;
  int j = *(const int *)b;

  return (scores[i] > scores[j]) - (scores[i] < scores[j]);
}

static int int16comp(const void *a, const void *b) {
  return *((uint16_t *)a) - *((uint16_t *)b);
}

typedef struct {
  poly_t poly;
  double ext_w;
  double ext_wm1;
  uint16_t ext_monomials[PARAM_OMEGA];
} offset_score_t;

static int compare_offset_ext_desc(const void *a, const void *b) {
  double sa = ((offset_score_t *)a)->ext_wm1;
  double sb = ((offset_score_t *)b)->ext_wm1;
  return (sb > sa) - (sb < sa);
}

static int valid_offset_soft(uint16_t *V, poly_t g, uint16_t delta,
                             double S[PARAM_N / 2][max_mult_y], int depth,
                             double *T) {
  int i, j, l, k = 0, r = g->length;
  int d = g->weight;
  double res[r];
  poly_t h = poly_alloc(r, d + 2);
  memcpy(h->index, g->index, d * sizeof(uint16_t));

  for (l = 0; l < r; ++l) {
    j = (l < r - l) ? l : (r - l);

    if ((l == 0) || (S[j - 1][0] != 1.)) {
      j = (l + delta) % r;

      // Check for duplicate positions
      bool duplicate = (l == j);

      for (i = 0; i < d && !duplicate; ++i) {
        if (g->index[i] == l || g->index[i] == j)
          duplicate = true;
      }

      if (duplicate) {
        res[l] = -10e10;
        continue;
      }

      h->index[d] = l;
      h->index[d + 1] = j;
      res[l] = log_score_poly(h, S);

    } else {
      res[l] = -10e10;
    }
  }

  int index[PARAM_N];

  for (i = 0; i < PARAM_N; i++)
    index[i] = i;
  qsort_r(index, PARAM_N, sizeof(int), compare_decrease, res);

  for (i = 0; i < r; i++)
    V[i] = index[i];
  qsort(res, r, sizeof(double), compare_floats);

  poly_free(h);
  /*
  Global threshold for the score of polynomials (below this score, it is not
  likely we are exploring a promising branch.)
  */
  if (res[0] < -90.)
    return 0;

  if (depth == 1)
    return 30;
  if (depth == 2)
    return 10;

  if (T != NULL && depth >= 3 && depth < 9) {
    if (res[0] <
        T[depth]) // thresholds for different sets of parameters in guess_key.h
      return 0;
  }

  for (i = 0; i < 10; i++) {
    if (res[i] - res[i + 1] > 1.5) {
      k = i + 1;
      break;
    }
  }

  if (k == 0) {
    for (i = 0; i < r - 1; i++) {
      if (res[i] < -2.) {
        k = i;
        break;
      }
    }
  }

  return k;
}

static poly_t guess_key_rec_soft(uint16_t *delta, poly_t *f, int depth,
                                 int max_depth, poly_t g,
                                 double S[PARAM_N / 2][max_mult_y], int weight,
                                 double gap_threshold, double *T) {

  int nb_offset, r = g->length;
  poly_t key = NULL;
  uint16_t V[r];

  if (depth >= max_depth)
    return NULL;

  if (f[depth] == NULL) { // initialize, only once per depth
    f[depth] = poly_alloc(r, 2);
    f[depth]->index[0] = 0;
    f[depth]->index[1] = delta[depth];
  }

  nb_offset = valid_offset_soft(V, g, delta[depth], S, depth, T);

  if (nb_offset <= 0)
    return NULL;

  offset_score_t *offset_scores =
      calloc((size_t)nb_offset, sizeof(*offset_scores));

  if (offset_scores == NULL)
    return NULL;

#pragma omp parallel for
  for (int i = 0; i < nb_offset; ++i) {
    double Ext[r];
    uint16_t Ext_monomials[r];

    poly_t p = poly_shift_union_soft(Ext, Ext_monomials, g, f[depth], V[i], S);

    offset_scores[i].poly = p;
    offset_scores[i].ext_w = (p != NULL) ? Ext[weight] : -1e11;
    offset_scores[i].ext_wm1 = (p != NULL) ? Ext[weight - 1] : -1e11;

    if (p != NULL) {
      memcpy(offset_scores[i].ext_monomials, Ext_monomials,
             (size_t)weight * sizeof(uint16_t));
    }
  }

  qsort(offset_scores, nb_offset, sizeof(*offset_scores),
        compare_offset_ext_desc);

  for (int i = 0; i < nb_offset; ++i) {
    offset_score_t *o = &offset_scores[i];
    poly_t p = o->poly;

    if (p == NULL)
      continue;

    if (o->ext_wm1 - o->ext_w > gap_threshold) {
      key = poly_alloc(r, weight);

      if (key != NULL) {
        memcpy(key->index, o->ext_monomials, (size_t)weight * sizeof(uint16_t));
      }
      break;
    }

    if (p->weight == weight) {
      key = poly_alloc(r, weight);

      if (key != NULL) {
        memcpy(key->index, p->index, (size_t)weight * sizeof(uint16_t));
      }
      break;
    }

    key = guess_key_rec_soft(delta, f, depth + 1, max_depth, p, S, weight,
                             gap_threshold, T);
    if (key != NULL)
      break;
  }

  for (int i = 0; i < nb_offset; ++i) {
    poly_free(offset_scores[i].poly);
  }

  free(offset_scores);

  return key;
}

static poly_t guess_key_soft(uint16_t *delta, double S[PARAM_N / 2][max_mult_y],
                             int r, int max_depth, int weight,
                             double gap_threshold, double *T) {
  int i;
  poly_t f[max_depth];
  poly_t h;

  for (i = 0; i < max_depth; ++i) {
    f[i] = NULL;
  }

  f[0] = poly_alloc(r, 2);
  f[0]->index[0] = 0;
  f[0]->index[1] = delta[0];

  h = guess_key_rec_soft(delta, f, 1, max_depth, f[0], S, weight, gap_threshold,
                         T);
  for (i = 0; i < max_depth; ++i) {
    if (f[i] != NULL)
      poly_free(f[i]);
  }
  return h;
}

static void sample_delta_soft(uint16_t *delta, double *S, int nb_pos) {

  int index[(PARAM_N) / 2];

  for (int i = 0; i < (PARAM_N) / 2; i++) {
    index[i] = i;
  }

  qsort_r(index, PARAM_N / 2, sizeof(int), compare, S);

  for (int i = 0; i < nb_pos; i++) {
    delta[i] = index[i] + 1;
  }
}

poly_t normalize(poly_t h) {
  int i, j, l, offset = 0, r = h->length, d = h->weight;

  poly_t h_normalized = poly_alloc(r, d);
  int8_t S[r / 2];

  spectrum(S, h->index, r, d);
  /* l = smallest distance with multiplicity 1 */
  /* conjecture: at least one distance has multiplicity 1 */
  for (l = 0; S[l] != 1; ++l)
    ;

  for (i = 0; i < d; ++i) {
    for (j = 0; j < i; ++j) {
      if (distance(h->index[i], h->index[j], r) == l) {
        break;
      }
    }

    if (i != j) {
      if (h->index[i] == (h->index[j] + l + 1) % r)
        offset = h->index[j];
      else
        offset = h->index[i];
      break;
    }
  }

  for (i = 0; i < d; ++i)
    h_normalized->index[i] = (r - offset + h->index[i]) % r;

  qsort(h_normalized->index, d, sizeof(uint16_t), int16comp);

  return h_normalized;
}

static poly_t mirror_poly(poly_t p) {
  int w = p->weight;
  int r = p->length;

  poly_t m = poly_alloc(r, w);
  if (!m)
    return NULL;

  for (int i = 0; i < w; i++) {
    m->index[i] = (r - p->index[i]) % r;
  }

  return m;
}

static bool is_key(poly_t guess_key, poly_t key_normalized) {
  if (guess_key == NULL)
    return false;
  poly_t guess_key_normalized = normalize(guess_key);
  bool success = true;

  for (int i = 0; i < PARAM_OMEGA; i++) {
    if (guess_key_normalized->index[i] != key_normalized->index[i]) {
      success = false;
      break;
    }
  }

  if (!success) {
    success = true;
    poly_t key_mirror = mirror_poly(guess_key);
    poly_t key_mirror_normalized = normalize(key_mirror);
    for (int i = 0; i < PARAM_OMEGA; i++) {
      if (key_mirror_normalized->index[i] != key_normalized->index[i]) {
        success = false;
        break;
      }
    }
    poly_free(key_mirror);
    poly_free(key_mirror_normalized);
  }

  poly_free(guess_key_normalized);

  return success;
}

int reconstruct_soft(double S_soft[PARAM_N / 2][max_mult_y], long int nb_dec,
                     double sigma, poly_t h_normalized) {
  double gap_threshold = 6.0;
  uint16_t delta[PARAM_OMEGA];

  sample_delta_soft(delta, means, PARAM_OMEGA);

  double T_l[9];
  double *T = NULL;

  if (sigma == 0 && nb_dec == 110000) {
    memcpy(T_l, thresholds_110k, sizeof(T_l));
    T = T_l;
  } else if (sigma == 0 && nb_dec == 190000) {
    memcpy(T_l, thresholds_190k, sizeof(T_l));
    T = T_l;
  }

  for (int i = 0; i < PARAM_N / 2; i++) {
    for (int j = 0; j < max_mult_y; j++) {
      if (S_soft[i][j] < 0 || S_soft[i][j] > 1) {
        printf("problem!\n");
      }
    }
  }
  printf("\n ********* Depth-first search reconstruction *********\n\n");
  double t0 = get_time();
  poly_t key = guess_key_soft(delta, S_soft, PARAM_N, PARAM_OMEGA, PARAM_OMEGA,
                              gap_threshold, T);
  double t1 = get_time();
  printf("\nTime:\t%g s\n", t1 - t0);

  if (key != NULL) {
    printf("Recovered key: ");
    for (int i = 0; i < PARAM_OMEGA; i++) {
      printf("%d ", key->index[i]);
    }

    if (is_key(key, h_normalized)) {
      printf("\n\n");
      printf("success soft\n");
      poly_free(key);
      return EXIT_SUCCESS;
    }

    else {
      poly_free(key);
      return EXIT_FAILURE;
    }

  } else {
    return EXIT_FAILURE;
  }
}