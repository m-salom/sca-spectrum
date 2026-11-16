#define _DEFAULT_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#include "annealing_soft_spectrum.h"
#include "api.h"
#include "err_attack.h"
#include "guess_key.h"
#include "hqc.h"
#include "model.h"
#include "parameters.h"
#include "parsing.h"
#include "poly.h"
#include "spectrum.h"
#include "symmetric.h"
#include "vector.h"

enum Mode { SA, DFS, BOTH };

int error_attack(double *syndromes, uint16_t *index, long nb_dec, double sigma,
                 uint64_t *h, uint64_t *s) {
  uint64_t syndrome[VEC_N_256_SIZE_64] = {0};
  uint8_t pk[PUBLIC_KEY_BYTES];
  uint8_t sk[SECRET_KEY_BYTES];
  uint8_t entropy[SEED_BYTES];
  __m256i y[VEC_N_256_NUM_WORDS] = {0};
  __m256i x[VEC_N_256_NUM_WORDS] = {0};
  __m256i r2[VEC_N_256_NUM_WORDS];
  __m256i r1[VEC_N_256_NUM_WORDS];
  __m256i e[VEC_N_256_NUM_WORDS];
  uint16_t y_sparse[PARAM_OMEGA];
  uint16_t r1_sparse[PARAM_OMEGA_R];
  uint8_t m[PARAM_SECURITY_BYTES] = {0};

  // Generate public et private key
  random_bytes(entropy, sizeof(entropy));
  prng_init(entropy, NULL, SEED_BYTES, 0);
  crypto_kem_keypair(pk, sk);
  hqc_ek_pke_from_string(h, s, pk);

  // Create seed
  uint8_t K_theta[SHARED_SECRET_BYTES + SEED_BYTES] = {0};
  uint8_t theta[SEED_BYTES] = {0};
  uint8_t hash_ek_kem[SEED_BYTES] = {0};
  ciphertext_kem_t c_kem_t = {0};
  shake256_xof_ctx theta_xof_ctx = {0};

  prng_get_bytes(c_kem_t.salt, SALT_BYTES);
  hash_h(hash_ek_kem, pk);
  hash_g(K_theta, hash_ek_kem, m, c_kem_t.salt);
  memcpy(theta, K_theta + SEED_BYTES, SEED_BYTES);
  xof_init(&theta_xof_ctx, theta, SEED_BYTES);

  // Recover secret key y
  recover_x_y(y, x, sk, pk);
  extract_bit_positions((uint64_t *)y, VEC_N_SIZE_64, y_sparse, PARAM_OMEGA);
  init_synd_attack(y_sparse, PARAM_N);

  for (int i = 0; i < PARAM_OMEGA; i++)
    index[i] = y_sparse[i];

  for (long int i = 0; i < nb_dec; ++i) {
    memset(r1, 0, sizeof(r1));
    memset(r2, 0, sizeof(r2));
    memset(e, 0, sizeof(e));

    vect_sample_fixed_weight2(&theta_xof_ctx, r1, PARAM_OMEGA_R);
    extract_bit_positions((uint64_t *)r1, VEC_N_SIZE_64, r1_sparse,
                          PARAM_OMEGA_R);

    vect_sample_fixed_weight2(&theta_xof_ctx, r2, PARAM_OMEGA_R);
    vect_sample_fixed_weight2(&theta_xof_ctx, e, PARAM_OMEGA_E);

    compute_syndrome(syndrome, x, y, e, r1, r2);

    double syndrome_weight = 0;
    for (int j = 0; j < VEC_N_SIZE_64; ++j) {
      syndrome_weight += __builtin_popcountll(syndrome[j]);
    }

    log_synd_attack(syndromes, r1_sparse, syndrome_weight, i, sigma);
  }

  display_synd_attack(entropy, sigma);
  return 0;
}

double get_time() {
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + t.tv_usec * 1e-6;
}

static void print_usage(const char *program) {
  printf("Usage: %s [OPTIONS]\n\n", program);

  printf("HQC secret-key reconstruction from soft spectrum.\n");

  printf("Options:\n");
  printf("  -h,  --help                  Show this help message\n");
  printf("  -m,  --mode                  Reconstruction mode: simulated "
         "annealing (SA) or depth first search (DFS) or BOTH (default: SA)\n");
  printf("  -N,  --samples               Number of observations (default: "
         "50000)\n");
  printf("  -s,  --sigma                 Measurement noise (default: 0)\n");
  printf("  -i,  --iterations            Iterations per restart (default: "
         "200000)\n");
  printf("  -r,  --restarts              Number of restarts (default: 30)\n");
  printf("  -e,  --exploration-restarts  Exploration restarts (default: 20)\n");
  printf("  -Ts, --start-temp            Initial temperature (default: 2.0)\n");
  printf(
      "  -Tf, --end-temp              Final temperature (default: 0.0001)\n");
  printf("  -k,  --init-top-k            Top-ranked initialization choices "
         "(default: 30)\n");
  printf("  -T,  --swap-trials           Guided move trials (default: 500)\n");
  printf("  -P,  --swap-proba            Random swap probability (default: "
         "0.8)\n");

  printf("Example:\n");
  printf("  %s -N 30000 -i 200000 -r 30\n", program);
}

int main(int argc, char **argv) {
  long int nb_dec = 50000; // default number of decoding samples
  double sigma = 0.0;
  enum Mode mode = SA;
  // Fixed model parameters
  double a = 1.0;
  double b = 0.0;
  double sigma0 = 94.4;
  srand(time(NULL));

  anneal_params_t params = {.iterations = 200000,
                            .restarts = 30,
                            .exploration_restarts = 20,
                            .start_temp = 2.0,
                            .end_temp = 0.0001,
                            .init_top_k = 30,
                            .swap_trials = 500,
                            .swap_proba = 0.8f,
                            .seed = (uint64_t)time(NULL)};

  // Optional command-line arguments
  for (int i = 1; i < argc; i++) {
    const char *arg = argv[i];

    // Display help and exit.
    if (strcmp(arg, "-h") == 0 || strcmp(arg, "--help") == 0) {
      print_usage(argv[0]);
      return EXIT_SUCCESS;
    }

    if (i + 1 >= argc) {
      fprintf(stderr, "Missing value for %s\n", arg);
      return EXIT_FAILURE;
    }

    const char *value = argv[++i];
    if (strcmp(arg, "-m") == 0 || strcmp(arg, "--mode") == 0) {
      if (strcmp(value, "SA") == 0) {
        mode = SA;
      } else if (strcmp(value, "DFS") == 0) {
        mode = DFS;
      } else if (strcmp(value, "BOTH") == 0) {
        mode = BOTH;
      } else {
        fprintf(stderr, "Invalid mode: %s (expected SA, DFS or BOTH)\n", value);
        return EXIT_FAILURE;
      }

    } else if (strcmp(arg, "-N") == 0 || strcmp(arg, "--samples") == 0) {
      nb_dec = atol(value);

    } else if (strcmp(arg, "-s") == 0 || strcmp(arg, "--sigma") == 0) {
      sigma = atof(value);

    } else if (strcmp(arg, "-i") == 0 || strcmp(arg, "--iterations") == 0) {
      params.iterations = atoll(value);

    } else if (strcmp(arg, "-r") == 0 || strcmp(arg, "--restarts") == 0) {
      params.restarts = atoi(value);

    } else if (strcmp(arg, "-e") == 0 ||
               strcmp(arg, "--exploration-restarts") == 0) {
      params.exploration_restarts = atoi(value);

    } else if (strcmp(arg, "-Ts") == 0 || strcmp(arg, "--start-temp") == 0) {
      params.start_temp = atof(value);

    } else if (strcmp(arg, "-Tf") == 0 || strcmp(arg, "--end-temp") == 0) {
      params.end_temp = atof(value);

    } else if (strcmp(arg, "-k") == 0 || strcmp(arg, "--init-top-k") == 0) {
      params.init_top_k = atoi(value);

    } else if (strcmp(arg, "-T") == 0 || strcmp(arg, "--swap-trials") == 0) {
      params.swap_trials = atoi(value);

    } else if (strcmp(arg, "-P") == 0 || strcmp(arg, "--swap-proba") == 0) {
      params.swap_proba = atof(value);

    } else {
      fprintf(stderr, "Unknown option: %s\n", arg);
      print_usage(argv[0]);
      return EXIT_FAILURE;
    }
  }

  if (nb_dec <= 0 || sigma < 0.0 || params.iterations <= 0 ||
      params.restarts <= 0 || params.exploration_restarts <= 0 ||
      params.exploration_restarts > params.restarts ||
      params.start_temp <= 0.0 || params.end_temp <= 0.0 ||
      params.init_top_k <= 0 || params.swap_trials <= 0 ||
      params.swap_proba < 0.0 || params.swap_proba > 1.0) {

    fprintf(stderr, "Invalid parameters.\n");
    return EXIT_FAILURE;
  }

  printf("\n ********* Collecting error weights *********\n\n");
  printf("Number of error weights = %ld\n", nb_dec);
  printf("Sigma                   = %f\n", sigma);

  double *syndromes = malloc(nb_dec * sizeof(double));
  if (syndromes == NULL) {
    perror("malloc failed");
    exit(1);
  }

  uint64_t h[VEC_N_256_SIZE_64] = {0};
  uint64_t s[VEC_N_256_SIZE_64] = {0};
  uint16_t index[PARAM_OMEGA];

  int status = error_attack(syndromes, index, nb_dec, sigma, h, s);

  if (status != EXIT_SUCCESS) {
    fprintf(stderr, "Failed to process observations.\n");
    free(syndromes);
    return EXIT_FAILURE;
  }

  poly_t key = poly_alloc(PARAM_N, PARAM_OMEGA);

  for (int i = 0; i < PARAM_OMEGA; i++) {
    key->index[i] = index[i];
  }

  printf("Secret key y : ");
  for (int i = 0; i < PARAM_OMEGA; i++) {
    printf("%d ", key->index[i]);
  }
  printf("\n");

  // Compute soft spectrum using model
  double S_soft[PARAM_N / 2][max_mult_y];
  get_soft_spectrum(S_soft, syndromes, nb_dec, means, sigma, a, b, sigma0);

  if (mode == SA || mode == BOTH) {
    // Simulated annealing
    double t0 = get_time();
    int sa_status = annealing(S_soft, h, s, &params);
    double t1 = get_time();

    printf("\nTime:\t%g s\n", t1 - t0);

    if (sa_status != EXIT_SUCCESS) {
      fprintf(stderr, "Simulated annealing failed.\n");
    }
  }

  if (mode == DFS || mode == BOTH) {
    poly_t key_normalized = normalize(key);
    // Depth-first search
    int dfs_status = reconstruct_soft(S_soft, nb_dec, sigma, key_normalized);
    poly_free(key_normalized);
    if (dfs_status != EXIT_SUCCESS) {
      fprintf(stderr, "Depth first reconstruction failed.\n");
    }
  }

  poly_free(key);
  free(syndromes);

  return 0;
}
