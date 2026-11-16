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

#define CCA_TOTAL_TRACES 410000L

enum Mode { SA, DFS, BOTH };

double get_time() {
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + t.tv_usec * 1e-6;
}

int cmp(const void *a, const void *b) {
  long x = *(long *)a;
  long y = *(long *)b;
  return (x > y) - (x < y);
}

int error_attack_cca(double *syndromes, uint16_t *index, long nb_dec,
                     uint64_t *h, uint64_t *s) {
  FILE *fo = fopen("data/cca/estimated_error_weights_cca.bin", "rb");

  if (!fo) {
    perror("fopen");
    return 1;
  }

  long total = CCA_TOTAL_TRACES;

  long *indices = malloc(total * sizeof(long));

  for (long i = 0; i < total; i++)
    indices[i] = i;

  for (long i = 0; i < nb_dec; i++) {
    long j = i + rand() % (total - i);
    long tmp = indices[i];
    indices[i] = indices[j];
    indices[j] = tmp;
  }

  qsort(indices, nb_dec, sizeof(long), cmp);

  float *weights = malloc(nb_dec * sizeof(float));
  long current = 0;
  long next_target = indices[current];
  float value;
  long idx = 0;

  while (fread(&value, sizeof(float), 1, fo) == 1) {
    if (idx == next_target) {
      weights[current] = value;
      current++;
      if (current == nb_dec)
        break;
      next_target = indices[current];
    }
    idx++;
  }
  fclose(fo);

  if (current != nb_dec) {
    fprintf(stderr, "Error: only %ld/%ld observations read\n", current, nb_dec);
    return EXIT_FAILURE;
  }

  uint8_t global_seed[SEED_BYTES] = {146, 150, 2,   31,  226, 68,  64,  230,
                                     105, 100, 81,  64,  90,  139, 122, 68,
                                     108, 67,  162, 191, 234, 48,  122, 242,
                                     18,  61,  26,  250, 71,  97,  250, 123};
  prng_init(global_seed, NULL, SEED_BYTES, 0);

  uint8_t pk[PUBLIC_KEY_BYTES];
  uint8_t sk[SECRET_KEY_BYTES];
  __m256i y[VEC_N_256_NUM_WORDS] = {0};
  __m256i x[VEC_N_256_NUM_WORDS] = {0};
  uint16_t y_sparse[PARAM_OMEGA];

  crypto_kem_keypair(pk, sk);
  hqc_ek_pke_from_string(h, s, pk);

  recover_x_y(y, x, sk, pk);

  extract_bit_positions((uint64_t *)y, VEC_N_SIZE_64, y_sparse, PARAM_OMEGA);

  init_synd_attack(y_sparse, PARAM_N);
  for (int i = 0; i < PARAM_OMEGA; i++)
    index[i] = y_sparse[i];

  FILE *f = fopen("data/cca/r1_vectors_cca.bin", "rb");
  if (f == NULL) {
    perror("Cannot open CCA r1 dataset");
    fprintf(stderr, "Download the dataset from Zenodo and place it "
                    "in data/cca/. See README.md.\n");
    return EXIT_FAILURE;
  }

  uint64_t r1[VEC_N_SIZE_64];
  uint16_t r1_sparse[PARAM_OMEGA_R];

  current = 0;
  idx = 0;
  next_target = indices[current];

  rewind(f);

  while (fread(r1, sizeof(uint64_t), VEC_N_SIZE_64, f) == VEC_N_SIZE_64) {
    if (idx == next_target) {
      extract_bit_positions(r1, VEC_N_SIZE_64, r1_sparse, PARAM_OMEGA_R);
      log_synd_attack(syndromes, r1_sparse, weights[current], current, 0);

      current++;
      if (current == nb_dec)
        break;

      next_target = indices[current];
    }

    idx++;
  }

  fclose(f);

  if (current != nb_dec) {
    fprintf(stderr, "Error: only %ld/%ld observations read\n", current, nb_dec);
    return EXIT_FAILURE;
  }

  uint8_t entropy[48];
  random_bytes(entropy, sizeof(entropy));
  display_synd_attack(entropy, 0);

  free(weights);
  free(indices);

  return 0;
}

static void print_usage(const char *program) {
  printf("Usage: %s [OPTIONS]\n\n", program);

  printf("HQC secret-key reconstruction from soft spectrum.\n");

  printf("Options:\n");
  printf("  -h,  --help                  Show this help message\n");
  printf("  -m,  --mode                  Reconstruction mode: simulated "
         "annealing (SA) or depth first search (DFS) or BOTH (default: SA)\n");
  printf("  -N,  --samples               Number of observations (default: "
         "50000 (max 410000))\n");
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

  if (nb_dec <= 0 || nb_dec > CCA_TOTAL_TRACES || params.iterations <= 0 ||
      params.restarts <= 0 || params.exploration_restarts <= 0 ||
      params.exploration_restarts > params.restarts ||
      params.start_temp <= 0.0 || params.end_temp <= 0.0 ||
      params.init_top_k <= 0 || params.swap_trials <= 0 ||
      params.swap_proba < 0.0 || params.swap_proba > 1.0) {

    fprintf(stderr, "Invalid parameters.\n");
    return EXIT_FAILURE;
  }

  printf("\n ********* Computing error weights *********\n\n");
  printf("Number of error weights = %ld\n", nb_dec);

  double *syndromes = malloc(nb_dec * sizeof(double));
  if (syndromes == NULL) {
    perror("malloc failed");
    exit(1);
  }

  uint64_t h[VEC_N_256_SIZE_64] = {0};
  uint64_t s[VEC_N_256_SIZE_64] = {0};
  uint16_t index[PARAM_OMEGA];

  int status = error_attack_cca(syndromes, index, nb_dec, h, s);

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
  get_soft_spectrum(S_soft, syndromes, nb_dec, means, 0, a, b, sigma0);

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
    int dfs_status = reconstruct_soft(S_soft, nb_dec, 0, key_normalized);
    poly_free(key_normalized);
    if (dfs_status != EXIT_SUCCESS) {
      fprintf(stderr, "Depth first reconstruction failed.\n");
    }
  }

  poly_free(key);
  free(syndromes);

  return 0;
}
