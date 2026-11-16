#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "err_attack.h"
#include "gf2x.h"
#include "hqc.h"
#include "parameters.h"
#include "parsing.h"
#include "spectrum.h"
#include "vector.h"

static bool save_logs = false;
static unsigned verbose_synd = 7;
static int max_mult;
static long nb_samples;
static double total_weight;
static long count_synd[PARAM_N / 2];
static double total_synd[PARAM_N / 2];
static double total_synd2[PARAM_N / 2];
static int8_t S_r[PARAM_N / 2];
int8_t S_y[PARAM_N / 2];
double means[PARAM_N / 2];

void init_synd_attack(uint16_t *y, int length) {
  int j;
  nb_samples = 0;
  total_weight = 0;
  max_mult = 0;

  spectrum(S_y, y, length, PARAM_OMEGA);

  memset(count_synd, 0, (PARAM_N / 2) * sizeof(long));
  memset(total_synd, 0, (PARAM_N / 2) * sizeof(double));
  memset(total_synd2, 0, (PARAM_N / 2) * sizeof(double));

  for (j = 0; j < length / 2; ++j)
    if (S_y[j] > max_mult)
      max_mult = S_y[j];

  if (max_mult >= max_mult_y) {
    fprintf(stderr, "Unsupported spectrum multiplicity: %d (max %d)\n",
            max_mult, max_mult_y - 1);
    exit(EXIT_FAILURE);
  }
}

static double normal_random() {
  double u1 = (rand() + 1.0) / (RAND_MAX + 1.0);
  double u2 = (rand() + 1.0) / (RAND_MAX + 1.0);
  return sqrt(-2.0 * log(u1)) * cos(2 * M_PI * u2);
}

// N(0, sigma^2)
static double normal_random_sigma(double sigma) {
  return sigma * normal_random();
}

void set_synd_attack_logging(bool enabled) { save_logs = enabled; }

#define WEIGHTED_SYND_W

void log_synd_attack(double *syndromes, uint16_t *e, double syndrome_weight,
                     long int nb_dec, double sigma) {
  int j;
  nb_samples++;

  if (sigma > 0) {
    double noise = normal_random_sigma(sigma);
    syndrome_weight += noise;
  }

  total_weight += syndrome_weight;

  spectrum(S_r, e, PARAM_N, PARAM_OMEGA_R);
  syndromes[nb_dec] = syndrome_weight;

  for (j = 0; j < PARAM_N / 2; ++j) {
    if (S_r[j] > 0) {
#ifdef WEIGHTED_SYND_W // weighted average syndrome weight
      count_synd[j] += S_r[j];
      total_synd[j] += S_r[j] * syndrome_weight;
      total_synd2[j] += S_r[j] * syndrome_weight * syndrome_weight;
#else // average syndrome weight
      count_synd[j]++;
      total_synd[j] += syndrome_weight;
      total_synd2[j] += syndrome_weight * syndrome_weight;
#endif
    }
  }
}

void display_synd_attack(uint8_t *entropy, double sigma) {
  int j;
  int blocklength = PARAM_N;
  char logfilename[100];
  long count_dist[max_mult_y];
  double avg_dist[max_mult_y];
  double std_dist[max_mult_y];
  unsigned long N = 0;
  double T = 0;
  FILE *file = NULL;

  if (save_logs) {
    sprintf(logfilename, "hqc_synd_attack_%ld_%d%d", nb_samples, entropy[0],
            entropy[1]);
    file = fopen(logfilename, "a");
    if (file == NULL) {
      perror("fopen");
      return;
    }
  }

  memset(count_dist, 0, max_mult_y * sizeof(long));
  memset(avg_dist, 0, max_mult_y * sizeof(double));
  memset(std_dist, 0, max_mult_y * sizeof(double));
  if (file != NULL) {
    fprintf(file, "# %d %d %d %ld\n", blocklength, PARAM_OMEGA,
            PARAM_OMEGA_E * 2, nb_samples);
    fprintf(file, "# sigma : %f\n", sigma);
    fprintf(file, "# entropy: ");
    for (int j = 0; j < SEED_BYTES; j++) {
      fprintf(file, "%d, ", entropy[j]);
    }
    fprintf(file, "\n");

    if (verbose_synd & 1) {
#ifdef WEIGHTED_SYND_W // weighted average syndrome weight
      fprintf(file, "# weighted average syndrome weight (weighted with "
                    "multiplicity in error)\n");
#else
      fprintf(file,
              "# average syndrome weight (when multiplicity in error > 0)\n");
#endif
    }
  }

  for (j = 0; j < blocklength / 2; ++j) {
    if (S_y[j] < max_mult_y) {
      double x = total_synd[j] / count_synd[j];
      avg_dist[S_y[j]] += x;
      std_dist[S_y[j]] += x * x;
      count_dist[S_y[j]]++;
      N += 1;
      T += x;
    }
  }

  T /= N;

  if (verbose_synd & 2) {
    for (j = 0; j <= max_mult; ++j) {
      if (count_dist[j] > 0) {
        avg_dist[j] /= count_dist[j];
        std_dist[j] /= count_dist[j];
        std_dist[j] -= avg_dist[j] * avg_dist[j];
        if (count_dist[j] == 1)
          std_dist[j] = 0;
        if (file != NULL) {
          fprintf(file, "# %d\t%-16.16g\t%5ld\t%-16.16g\n", j, avg_dist[j],
                  count_dist[j], sqrt(std_dist[j]));
        }
      }
    }
    if (file != NULL) {
      fprintf(file, "# \t\t%-16.16g\t%5ld\n", T, N);
      fprintf(file, "# \t\t%-16.16g\n", total_weight / nb_samples);
    }
  }

  for (j = 0; j < blocklength / 2; ++j) {
    double x = total_synd[j] / count_synd[j];
    double y = total_synd2[j] / count_synd[j];
    y -= x * x;
    means[j] = x;

    if ((verbose_synd & 4) && file != NULL)
      fprintf(file, "%-16.16g\t%-16.16g\t%d\t%d\t%d\n", x, sqrt(y), 0, j,
              S_y[j]);
  }

  if (file != NULL) {
    fclose(file);
  }
}

// Regenerate the secret polynomials from HQC secret and public keys.
void recover_x_y(__m256i *y, __m256i *x, unsigned char *sk, unsigned char *pk) {
  // Retrieve y from secret key
  uint8_t dk_pke[SEED_BYTES] = {0};
  memcpy(dk_pke, sk + PUBLIC_KEY_BYTES, SEED_BYTES);
  hqc_dk_pke_from_string(y, dk_pke);

  // Retrieve x
  __m256i h[VEC_N_256_NUM_WORDS] = {0};
  uint64_t s[VEC_N_256_SIZE_64] = {0};
  __m256i tmp_s[VEC_N_256_NUM_WORDS] = {0};

  hqc_ek_pke_from_string((uint64_t *)h, s, pk);
  vect_mul(tmp_s, h, y);
  vect_add((uint64_t *)x, (uint64_t *)tmp_s, s, VEC_N_256_SIZE_64);
}

// Compute r2*x + r1*y + e in the binary cyclic polynomial ring.
void compute_syndrome(uint64_t *s, __m256i *x, __m256i *y, __m256i *e,
                      __m256i *r1, __m256i *r2) {
  __m256i tmp1[VEC_N_256_NUM_WORDS] = {0};
  __m256i tmp2[VEC_N_256_NUM_WORDS] = {0};
  uint64_t tmp3[VEC_N_256_SIZE_64] = {0};

  vect_mul(tmp1, r2, x);
  vect_mul(tmp2, r1, y);
  vect_add(tmp3, (uint64_t *)tmp1, (uint64_t *)tmp2, VEC_N_256_SIZE_64);
  vect_add(s, tmp3, (uint64_t *)e, VEC_N_256_SIZE_64);
}

void extract_bit_positions(const uint64_t *input, long int input_len,
                           uint16_t *output, long int max_output_len) {
  long int out_index = 0;
  for (long int word = 0; word < input_len; word++) {
    uint64_t value = input[word];
    for (int bit = 0; bit < 64; bit++) {
      if ((value >> bit) & 1) {
        if (out_index >= max_output_len) {
          fprintf(stderr,
                  "Erreur : poids de Hamming trop élevé (%ld bits à 1, max = "
                  "%ld)\n",
                  out_index, max_output_len);
          exit(EXIT_FAILURE);
        }
        output[out_index++] = word * 64 + bit;
      }
    }
  }
}