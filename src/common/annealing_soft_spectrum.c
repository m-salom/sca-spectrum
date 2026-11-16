#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "annealing_soft_spectrum.h"
#include "gf2x.h"
#include "hqc.h"
#include "spectrum.h"

#define N_DIST (PARAM_N / 2)
#define EPS_PROB 1e-300
#define ELITE_SIZE 20

typedef struct {
  uint64_t state;
} rng_t;

static uint64_t splitmix64_next(rng_t *rng) {
  uint64_t z = (rng->state += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

static double rng_uniform01(rng_t *rng) {
  uint64_t x = splitmix64_next(rng);
  return (x >> 11) * (1.0 / 9007199254740992.0);
}

static size_t rng_index(rng_t *rng, size_t n) {
  return (size_t)(rng_uniform01(rng) * (double)n);
}

typedef struct {
  double *data;
} soft_spectrum_t;

static inline double soft_get(const soft_spectrum_t *S, int d, int m) {
  if (d < 0 || d >= N_DIST) {
    return EPS_PROB;
  }

  if (m < 0 || m >= max_mult_y) {
    return EPS_PROB;
  }

  double p = S->data[(size_t)d * max_mult_y + (size_t)m];

  if (!isfinite(p) || p <= 0.0) {
    return EPS_PROB;
  }

  return p;
}

typedef struct {
  uint16_t support[PARAM_OMEGA];
  uint8_t *in_support;
  uint16_t *mu;
  double score;
} state_t;

static bool state_alloc(state_t *st) {
  st->in_support = (uint8_t *)calloc(PARAM_N, sizeof(uint8_t));
  st->mu = (uint16_t *)calloc(N_DIST, sizeof(uint16_t));

  if (st->in_support == NULL || st->mu == NULL) {
    free(st->in_support);
    free(st->mu);
    st->in_support = NULL;
    st->mu = NULL;
    return false;
  }

  st->score = -INFINITY;
  return true;
}

static void state_free(state_t *st) {
  free(st->in_support);
  free(st->mu);
  st->in_support = NULL;
  st->mu = NULL;
}

/** Compute score of a candidate from its distance multiplicities
 *  as the sum of log posterior probabilities over distances. */
static double score_from_mu(const uint16_t *mu, const soft_spectrum_t *S) {
  double score = 0.0;

  for (int d = 0; d < N_DIST; ++d) {
    int m = (int)mu[d];
    double p = soft_get(S, d, m);
    score += log(p);
  }

  return score;
}

static int compare_u16(const void *a, const void *b) {
  uint16_t x = *(const uint16_t *)a;
  uint16_t y = *(const uint16_t *)b;
  return (x > y) - (x < y);
}

typedef struct {
  uint16_t support_idx[PARAM_OMEGA];
  double score;
  bool valid;
} EliteState;

static bool same_support_idx(const uint16_t a[PARAM_OMEGA],
                             const uint16_t b[PARAM_OMEGA]) {
  uint16_t a_sorted[PARAM_OMEGA];
  uint16_t b_sorted[PARAM_OMEGA];

  memcpy(a_sorted, a, PARAM_OMEGA * sizeof(uint16_t));
  memcpy(b_sorted, b, PARAM_OMEGA * sizeof(uint16_t));

  qsort(a_sorted, PARAM_OMEGA, sizeof(uint16_t), compare_u16);
  qsort(b_sorted, PARAM_OMEGA, sizeof(uint16_t), compare_u16);

  return memcmp(a_sorted, b_sorted, PARAM_OMEGA * sizeof(uint16_t)) == 0;
}

/** Store high-scoring candidate supports. */
static void add_to_elite(EliteState elite[ELITE_SIZE], const state_t *state) {
  double reference_score = state->score;

  for (size_t i = 0; i < ELITE_SIZE; i++) {
    if (elite[i].valid &&
        same_support_idx(elite[i].support_idx, state->support)) {
      if (reference_score > elite[i].score) {
        memcpy(elite[i].support_idx, state->support,
               PARAM_OMEGA * sizeof(uint16_t));
        elite[i].score = state->score;
      }
      return;
    }
  }

  size_t target = SIZE_MAX;

  for (size_t i = 0; i < ELITE_SIZE; i++) {
    if (!elite[i].valid) {
      target = i;
      break;
    }
  }

  if (target == SIZE_MAX) {
    target = 0;

    for (size_t i = 1; i < ELITE_SIZE; i++) {
      if (elite[i].score < elite[target].score) {
        target = i;
      }
    }

    if (reference_score <= elite[target].score) {
      return;
    }
  }

  memcpy(elite[target].support_idx, state->support,
         PARAM_OMEGA * sizeof(uint16_t));
  elite[target].score = state->score;
  elite[target].valid = true;
}

static size_t choose_valid_elite(const EliteState elite[ELITE_SIZE],
                                 rng_t *rng) {
  uint16_t valid_indices[ELITE_SIZE];
  size_t n_valid = 0;

  for (size_t i = 0; i < ELITE_SIZE; i++) {
    if (elite[i].valid) {
      valid_indices[n_valid++] = i;
    }
  }

  if (n_valid == 0) {
    return SIZE_MAX;
  }

  return valid_indices[rng_index(rng, n_valid)];
}

static size_t choose_best_elite(const EliteState elite[ELITE_SIZE]) {
  uint16_t valid_indices[ELITE_SIZE];
  size_t n_valid = 0;

  for (size_t i = 0; i < ELITE_SIZE; i++) {
    if (elite[i].valid) {
      valid_indices[n_valid++] = i;
    }
  }

  if (n_valid == 0) {
    return SIZE_MAX;
  }

  double max_score = -INFINITY;
  uint16_t best_idx = 0;

  for (size_t i = 0; i < n_valid; i++) {
    if (elite[valid_indices[i]].score > max_score) {
      max_score = elite[valid_indices[i]].score;
      best_idx = valid_indices[i];
    }
  }

  return best_idx;
}

static void compute_spectrum(state_t *state) {
  memset(state->mu, 0, N_DIST * sizeof(uint16_t));

  for (size_t i = 0; i < PARAM_OMEGA; i++) {
    uint16_t p1 = state->support[i];

    for (size_t j = i + 1; j < PARAM_OMEGA; j++) {
      uint16_t p2 = state->support[j];
      int d = distance(p1, p2, PARAM_N);
      state->mu[d]++;
    }
  }
}

static void load_elite_state(state_t *state, const EliteState *elite_state,
                             const soft_spectrum_t *S) {
  memset(state->in_support, 0, PARAM_N * sizeof(bool));

  for (size_t i = 0; i < PARAM_OMEGA; i++) {
    state->support[i] = elite_state->support_idx[i];
    state->in_support[state->support[i]] = true;
  }

  compute_spectrum(state);
  state->score = score_from_mu(state->mu, S);
}

/** Calculate the score increment for inserting one support position.
 *  Only the affected distance multiplicities are evaluated. */
static double delta_add_position(uint16_t candidate, const uint16_t *support,
                                 int support_size, const uint16_t *mu,
                                 const soft_spectrum_t *S) {
  int touched[PARAM_OMEGA];
  int touched_count = 0;
  uint16_t increments[PARAM_OMEGA];

  for (int i = 0; i < support_size; ++i) {
    int d = distance(candidate, support[i], PARAM_N);
    if (d < 0 || d >= N_DIST) {
      continue;
    }

    int pos = -1;
    for (int k = 0; k < touched_count; ++k) {
      if (touched[k] == d) {
        pos = k;
        break;
      }
    }

    if (pos == -1) {
      touched[touched_count] = d;
      increments[touched_count] = 1;
      touched_count++;
    }

    else {
      increments[pos]++;
    }
  }

  double delta = 0.0;

  for (int k = 0; k < touched_count; ++k) {
    int d = touched[k];
    int old_m = mu[d];
    int new_m = old_m + increments[k];

    delta += log(soft_get(S, d, new_m));
    delta -= log(soft_get(S, d, old_m));
  }

  return delta;
}

typedef struct {
  uint16_t pos;
  double delta;
} candidate_t;

static int compare_candidate_desc(const void *a, const void *b) {
  const candidate_t *ca = (const candidate_t *)a;
  const candidate_t *cb = (const candidate_t *)b;

  if (ca->delta < cb->delta) {
    return 1;
  }
  if (ca->delta > cb->delta) {
    return -1;
  }
  return 0;
}

/** Build an initial support greedily, sampling among top-ranked additions.
 *  Position zero fixes the arbitrary cyclic-rotation symmetry. */
static bool initialize(state_t *st, const soft_spectrum_t *S, rng_t *rng,
                       int top_k) {
  memset(st->in_support, 0, PARAM_N * sizeof(uint8_t));
  memset(st->mu, 0, N_DIST * sizeof(uint16_t));

  st->support[0] = 0;
  st->in_support[0] = 1;
  int support_size = 1;

  candidate_t *candidates =
      (candidate_t *)malloc((PARAM_N - 1) * sizeof(candidate_t));
  if (candidates == NULL) {
    return false;
  }

  while (support_size < PARAM_OMEGA) {
    int n_candidates = 0;
    /* Evaluate all absent positions and rank them according
     * to their incremental contribution to the score. */
    for (int p = 1; p < PARAM_N; ++p) {
      if (st->in_support[p]) {
        continue;
      }

      candidates[n_candidates].pos = (uint16_t)p;
      candidates[n_candidates].delta =
          delta_add_position((uint16_t)p, st->support, support_size, st->mu, S);
      n_candidates++;
    }

    qsort(candidates, (size_t)n_candidates, sizeof(candidate_t),
          compare_candidate_desc);

    if (top_k > n_candidates) {
      top_k = n_candidates;
    }

    /*
     * Randomly select one of the top_k positions using
     * exponentially decreasing weights based on their rank.
     */
    const double rank_temperature = 5.0;
    double total_weight = 0.0;

    for (int k = 0; k < top_k; ++k) {
      total_weight += exp(-(double)k / rank_temperature);
    }

    double u = rng_uniform01(rng) * total_weight;
    int chosen_rank = top_k - 1;
    double acc = 0.0;

    for (int k = 0; k < top_k; ++k) {
      acc += exp(-(double)k / rank_temperature);
      if (u <= acc) {
        chosen_rank = k;
        break;
      }
    }

    uint16_t chosen = candidates[chosen_rank].pos;

    for (int i = 0; i < support_size; ++i) {
      int d = distance(chosen, st->support[i], PARAM_N);
      if (d >= 0 && d < N_DIST) {
        st->mu[d]++;
      }
    }

    st->support[support_size] = chosen;
    st->in_support[chosen] = 1;
    support_size++;
  }

  free(candidates);

  st->score = score_from_mu(st->mu, S);
  return true;
}

/** Evaluate replacing one position without modifying the state. */
static double swap_delta(const state_t *st, int remove_idx, uint16_t add_pos,
                         const soft_spectrum_t *S) {
  uint16_t remove_pos = st->support[remove_idx];
  int touched[2 * PARAM_OMEGA];
  int delta_mu[2 * PARAM_OMEGA];
  int touched_count = 0;

  for (int i = 0; i < PARAM_OMEGA; ++i) {
    if (i == remove_idx) {
      continue;
    }

    uint16_t other = st->support[i];
    int d_remove = distance(remove_pos, other, PARAM_N);
    int d_add = distance(add_pos, other, PARAM_N);

    if (d_remove >= 0 && d_remove < N_DIST) {
      int pos = -1;
      for (int k = 0; k < touched_count; ++k) {
        if (touched[k] == d_remove) {
          pos = k;
          break;
        }
      }

      if (pos == -1) {
        touched[touched_count] = d_remove;
        delta_mu[touched_count] = -1;
        touched_count++;
      }

      else {
        delta_mu[pos]--;
      }
    }

    if (d_add >= 0 && d_add < N_DIST) {
      int pos = -1;
      for (int k = 0; k < touched_count; ++k) {
        if (touched[k] == d_add) {
          pos = k;
          break;
        }
      }

      if (pos == -1) {
        touched[touched_count] = d_add;
        delta_mu[touched_count] = +1;
        touched_count++;
      }

      else {
        delta_mu[pos]++;
      }
    }
  }

  double delta_score = 0.0;

  for (int k = 0; k < touched_count; ++k) {
    int d = touched[k];
    int old_m = st->mu[d];
    int new_m = old_m + delta_mu[k];

    if (new_m < 0) {
      return -INFINITY;
    }

    delta_score += log(soft_get(S, d, new_m));
    delta_score -= log(soft_get(S, d, old_m));
  }

  return delta_score;
}

/** Apply a swap and update the distance spectrum accordingly. */
static void swap(state_t *st, int remove_idx, uint16_t add_pos) {
  uint16_t remove_pos = st->support[remove_idx];

  for (int i = 0; i < PARAM_OMEGA; ++i) {
    if (i == remove_idx) {
      continue;
    }

    uint16_t other = st->support[i];

    int d_remove = distance(remove_pos, other, PARAM_N);
    if (d_remove >= 0 && d_remove < N_DIST) {
      st->mu[d_remove]--;
    }

    int d_add = distance(add_pos, other, PARAM_N);
    if (d_add >= 0 && d_add < N_DIST) {
      st->mu[d_add]++;
    }
  }

  st->in_support[remove_pos] = 0;
  st->in_support[add_pos] = 1;
  st->support[remove_idx] = add_pos;
}

static uint16_t random_absent_position(const state_t *st, rng_t *rng) {
  for (;;) {
    uint16_t p = (uint16_t)(rng_index(rng, PARAM_N));
    if (!st->in_support[p]) {
      return p;
    }
  }
}

static void copy_state_solution(uint16_t dst_support[PARAM_OMEGA],
                                const state_t *src) {
  memcpy(dst_support, src->support, PARAM_OMEGA * sizeof(uint16_t));
}

static void print_support(const uint16_t support[PARAM_OMEGA]) {
  uint16_t tmp[PARAM_OMEGA];
  memcpy(tmp, support, sizeof(tmp));

  qsort(tmp, PARAM_OMEGA, sizeof(uint16_t), compare_u16);

  for (int i = 0; i < PARAM_OMEGA; ++i) {
    printf("%u", (unsigned)tmp[i]);
    if (i + 1 < PARAM_OMEGA) {
      printf(" ");
    }
  }
  printf("\n");
}

/** Evaluate random absent positions and select the best move. */
static uint16_t guided_add_position(const state_t *st, int remove_idx,
                                    rng_t *rng, const soft_spectrum_t *S,
                                    int n_trials, double *best_delta) {
  uint16_t best_pos = 0;
  *best_delta = -INFINITY;

  for (int k = 0; k < n_trials; k++) {
    uint16_t p = random_absent_position(st, rng);

    double delta = swap_delta(st, remove_idx, p, S);

    if (delta > *best_delta) {
      *best_delta = delta;
      best_pos = p;
    }
  }

  return best_pos;
}

static void multiply_h_by_sparse_y(const uint64_t h[VEC_N_256_SIZE_64],
                                   const uint16_t y[PARAM_OMEGA],
                                   uint64_t hy[VEC_N_256_SIZE_64]) {
  __m256i y_avx[VEC_N_256_NUM_WORDS] = {0};
  __m256i hy_avx[VEC_N_256_NUM_WORDS] = {0};

  uint64_t *y_words = (uint64_t *)y_avx;

  for (int i = 0; i < PARAM_OMEGA; i++) {
    uint16_t pos = y[i];
    y_words[pos >> 6] |= 1ULL << (pos & 63);
  }

  __m256i h_avx[VEC_N_256_NUM_WORDS] = {0};
  memcpy(h_avx, h, VEC_N_SIZE_64 * sizeof(uint64_t));

  vect_mul(hy_avx, h_avx, y_avx);

  memcpy(hy, hy_avx, sizeof(hy_avx));
}

static void rotate_left_one_u64(uint64_t v[VEC_N_256_SIZE_64]) {
  const int last_word = (PARAM_N - 1) >> 6;
  const int last_bit = (PARAM_N - 1) & 63;

  uint64_t wrap = (v[last_word] >> last_bit) & 1ULL;

  uint64_t carry = 0;

  for (int i = 0; i <= last_word; i++) {
    uint64_t current = v[i];
    uint64_t next_carry = current >> 63;
    v[i] = (current << 1) | carry;
    carry = next_carry;
  }

  if (wrap) {
    v[0] |= 1ULL;
  }

  if ((PARAM_N & 63) != 0) {
    uint64_t mask = (1ULL << (PARAM_N & 63)) - 1ULL;
    v[last_word] &= mask;
  }
}

static int public_weight_from_hy_u64(const uint64_t hy[VEC_N_256_SIZE_64],
                                     const uint64_t s[VEC_N_256_SIZE_64]) {
  int weight = 0;

  for (int i = 0; i < VEC_N_256_SIZE_64; i++) {
    uint64_t x = s[i] ^ hy[i];

    if (i == VEC_N_256_SIZE_64 - 1 && (PARAM_N & 63) != 0) {
      uint64_t mask = (1ULL << (PARAM_N & 63)) - 1ULL;

      x &= mask;
    }

    weight += __builtin_popcountll(x);
  }

  return weight;
}

/** Validate a candidate using the public-key relation s = x + h*y.
 *  Test cyclic shifts and the mirrored support; the residual must have
 *  the expected Hamming weight. */
static bool is_key(const uint16_t candidate[PARAM_OMEGA], uint64_t *h,
                   uint64_t *s, int *found_shift, bool *found_mirror,
                   uint16_t recovered_y[PARAM_OMEGA]) {
  uint64_t hy[VEC_N_256_SIZE_64];

  multiply_h_by_sparse_y(h, candidate, hy);

  for (int shift = 0; shift < PARAM_N; shift++) {
    int weight = public_weight_from_hy_u64(hy, s);

    if (weight == PARAM_OMEGA) {
      *found_shift = shift;
      *found_mirror = false;

      for (int i = 0; i < PARAM_OMEGA; i++) {
        recovered_y[i] = (candidate[i] + shift) % PARAM_N;
      }

      return true;
    }

    rotate_left_one_u64(hy);
  }

  uint16_t mirror_support[PARAM_OMEGA];

  for (int i = 0; i < PARAM_OMEGA; i++) {
    mirror_support[i] = (PARAM_N - candidate[i]) % PARAM_N;
  }

  multiply_h_by_sparse_y(h, mirror_support, hy);

  for (int shift = 0; shift < PARAM_N; shift++) {
    int weight = public_weight_from_hy_u64(hy, s);

    if (weight == PARAM_OMEGA) {
      *found_shift = shift;
      *found_mirror = true;

      for (int i = 0; i < PARAM_OMEGA; i++) {
        recovered_y[i] = (mirror_support[i] + shift) % PARAM_N;
      }

      return true;
    }

    rotate_left_one_u64(hy);
  }

  return false;
}

/** Run simulated annealing exploration followed by restarts from elite
 * candidates. Returns true only when public-key validation succeeds. */
static bool simulated_annealing(const soft_spectrum_t *S,
                                const anneal_params_t *params,
                                uint16_t best_support[PARAM_OMEGA],
                                double *best_score_out, uint64_t *h,
                                uint64_t *s) {
  state_t st;
  if (!state_alloc(&st)) {
    fprintf(stderr, "Cannot allocate state.\n");
    return false;
  }

  rng_t rng;
  rng.state = params->seed;
  EliteState elite[ELITE_SIZE];
  memset(elite, 0, sizeof(elite));
  double global_best_score = -INFINITY;

  for (int restart = 0; restart < params->restarts; ++restart) {
    bool elite_restart = restart >= params->exploration_restarts;
    /* Initialize candidate support from exploration restart
       using the most promising positions according to the estimated
       distance spectrum */
    if (!elite_restart) {
      if (!initialize(&st, S, &rng, params->init_top_k)) {
        state_free(&st);
        return false;
      }
    }
    /* Initialize candidate support by loading an elite state */
    else {
      if (restart == params->exploration_restarts) {
        size_t elite_idx = choose_best_elite(elite);
        if (elite_idx == SIZE_MAX) {
          fprintf(stderr, "No elite candidate available for restart.\n");
          state_free(&st);
          return false;
        }
        load_elite_state(&st, &elite[elite_idx], S);
        printf("Restart %u ELITE from elite %lu\n", restart, elite_idx);
      } else {
        size_t elite_idx = choose_valid_elite(elite, &rng);
        if (elite_idx != SIZE_MAX) {
          load_elite_state(&st, &elite[elite_idx], S);
          printf("Restart %u ELITE from elite %lu\n", restart, elite_idx);
        }
      }
    }

    add_to_elite(elite, &st);
    double local_best_score = st.score;
    uint16_t local_best_support[PARAM_OMEGA];
    copy_state_solution(local_best_support, &st);

    for (long long it = 0; it < params->iterations; ++it) {
      /* Update temperature using a geometric cooling */
      double progress = (params->iterations <= 1)
                            ? 1.0
                            : (double)it / (double)(params->iterations - 1);
      double temp = params->start_temp *
                    pow(params->end_temp / params->start_temp, progress);
      double delta;
      int remove_idx;
      uint16_t add_pos;

      if (!elite_restart) {
        /* Evaluate score difference from a random swap for exploration restarts
         */
        remove_idx = (int)rng_index(&rng, PARAM_OMEGA);
        add_pos = random_absent_position(&st, &rng);
        delta = swap_delta(&st, remove_idx, add_pos, S);
      }

      else {
        remove_idx = (int)rng_index(&rng, PARAM_OMEGA);
        /* With probability swap_proba, evaluate score difference from a random
           swap for elite restarts. */
        if (elite_restart && rng_uniform01(&rng) < params->swap_proba) {
          add_pos = random_absent_position(&st, &rng);
          delta = swap_delta(&st, remove_idx, add_pos, S);
        }
        /* Else, evaluate score differences between one random position inside
           the support and swap_trials different positions outside and choose
           the best move for elite restarts. */
        else {
          add_pos = guided_add_position(&st, remove_idx, &rng, S,
                                        params->swap_trials, &delta);
        }
      }
      bool accept = false;

      /* Swaps that increase the score are always accepted. */
      if (delta >= 0.0) {
        accept = true;
      }

      /* Swaps that decrease the score are accepted with
         probability depending on the current temperature. */
      else if (isfinite(delta) && temp > 0.0) {
        double prob = exp(delta / temp);
        if (rng_uniform01(&rng) < prob) {
          accept = true;
        }
      }

      if (accept) {
        swap(&st, remove_idx, add_pos);
        st.score += delta;

        if (st.score > local_best_score) {
          /* Add supports with highest scores to elite states.*/
          add_to_elite(elite, &st);
          local_best_score = st.score;
          copy_state_solution(local_best_support, &st);
          int shift;
          bool mirror;
          uint16_t validated_y[PARAM_OMEGA];
          /* Check solution using public-key verification. */
          bool success =
              is_key(local_best_support, h, s, &shift, &mirror, validated_y);

          if (success) {
            printf("Best support = ");
            print_support(local_best_support);
            memcpy(best_support, validated_y, PARAM_OMEGA * sizeof(uint16_t));

            if (mirror) {
              printf("\nSUCCESS: public validation, mirror + shift %d\n",
                     shift);
            }

            else {
              printf("\nSUCCESS: public validation, shift %d\n", shift);
            }

            state_free(&st);
            *best_score_out = local_best_score;

            return true;
          }
        }
      }
    }

    printf("restart %d/%d: init+SA best_score = %.12f\n", restart + 1,
           params->restarts, local_best_score);

    if (local_best_score > global_best_score) {
      global_best_score = local_best_score;
      memcpy(best_support, local_best_support, PARAM_OMEGA * sizeof(uint16_t));
    }
  }

  state_free(&st);
  *best_score_out = global_best_score;
  return false;
}

/** Run simulated annealing based key reconstruction from a soft distance
 * spectrum.
 *
 * @param S_soft Soft spectrum containing probabilities indexed by distance and
 * multiplicity.
 * @param h Public key h for verification.
 * @param s Public key s for verification.
 * @param params Empirical search parameters.
 * @return EXIT_SUCCESS if a candidate passes public-key validation.
 */
int annealing(double S_soft[PARAM_N / 2][max_mult_y], uint64_t *h, uint64_t *s,
              const anneal_params_t *params) {

  printf("\n ********* Simulated Annealing Reconstruction *********\n\n");

  printf("iterations        = %lld\n", params->iterations);
  printf("restarts          = %d\n", params->restarts);
  printf("start_temp        = %.8g\n", params->start_temp);
  printf("end_temp          = %.8g\n", params->end_temp);
  printf("init_top_k        = %d\n", params->init_top_k);
  printf("swap trials       = %d\n", params->swap_trials);
  printf("swap probability  = %f\n", params->swap_proba);
  printf("seed              = %llu\n", (unsigned long long)params->seed);

  soft_spectrum_t S;
  S.data = &S_soft[0][0];

  uint16_t best_support[PARAM_OMEGA];
  double best_score = -INFINITY;

  bool ok = simulated_annealing(&S, params, best_support, &best_score, h, s);

  if (!ok) {
    return EXIT_FAILURE;
  }

  printf("\nBest_score = %.12f\n", best_score);
  printf("Recovered_secret y = ");
  print_support(best_support);

  return EXIT_SUCCESS;
}
