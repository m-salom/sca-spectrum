#include "model.h"
#include <math.h>

/* log(n choose k) */
static long double log_nCk_ld(int n, int k) {
  if (k < 0 || k > n)
    return -INFINITY;
  if (k > n - k)
    k = n - k;
  return lgammal((long double)n + 1.0L) - lgammal((long double)k + 1.0L) -
         lgammal((long double)(n - k) + 1.0L);
}

static long double S_hqc(int n, int w, int t) {
  int limit = (w < t) ? w : t;
  long double max_log = -INFINITY;

  for (int i = 1; i <= limit; i += 2) {
    int ti = t - i;
    if (ti < 0 || ti > n - w)
      continue;
    long double lnum = log_nCk_ld(w, i) + log_nCk_ld(n - w, ti);
    if (lnum > max_log)
      max_log = lnum;
  }

  if (max_log == -INFINITY)
    return 0.0L;

  long double sum = 0.0L;
  for (int i = 1; i <= limit; i += 2) {
    int ti = t - i;
    if (ti < 0 || ti > n - w)
      continue;
    long double lnum = log_nCk_ld(w, i) + log_nCk_ld(n - w, ti);
    sum += expl(lnum - max_log);
  }

  long double logsumexp = max_log + logl(sum);
  long double logden = log_nCk_ld(n, t);
  long double result = expl(logsumexp - logden);

  return result;
}

static double WHQC(long int n, int w, int t, int nu) {
  return ((n - 2 * (w - nu) - nu) * S_hqc(n - 2, w, t - 2) +
          2 * (w - nu) * (1 - S_hqc(n - 2, w - 1, t - 2)) +
          nu * S_hqc(n - 2, w - 2, t - 2));
}

/*
  Note: in the paper, the model is adjusted by centering the observations
  subtracting the observed average syndrome weight and the model subtracting the
  average syndrome weight. Here, we only center the observations subtracting
  both the average syndrome weight and the observed average syndrome weight but
  it amounts to the same thing.
*/
void shifting_values(double *syndromes, long int nb_dec, double a, double b) {
  long double som = 0.0;
  long double Shqc = S_hqc(PARAM_N, PARAM_OMEGA, PARAM_OMEGA_R);
  double w_prime = PARAM_OMEGA_E + (PARAM_N - 2 * PARAM_OMEGA_E) * Shqc;
  double coeff = (1.0 / nb_dec);
  double term = (1 - 2 * w_prime / PARAM_N) * PARAM_N * Shqc;

  for (long int i = 0; i < nb_dec; i++)
    som += coeff * (a * (w_prime + term) + b) - coeff * syndromes[i];

  for (int i = 0; i < PARAM_N / 2; i++) {
    means[i] += som;
  }
}

/*
  Compute parameters of the model: mean, standard deviation and multiplicity
  distribution.
*/
void compute_model(double avg_dist_model[max_mult_y],
                   double std_dist_model[max_mult_y],
                   long double count_dist_model[max_mult_y], double noise,
                   long int nb_dec, double a, double b, double sigma0) {
  double w_prime =
      PARAM_OMEGA_E + (PARAM_N - 2 * PARAM_OMEGA_E) *
                          S_hqc(PARAM_N, PARAM_OMEGA, PARAM_OMEGA_R);
  double coeff = (1 - 2 * w_prime / PARAM_N);
  double std_val = (sigma0 * sigma0 + noise * noise) / nb_dec;
  long double logC_N1 = log_nCk_ld(PARAM_N - 1, PARAM_OMEGA - 1);

  for (int m = 0; m < max_mult_y; m++) {
    std_dist_model[m] = std_val;
    avg_dist_model[m] =
        a * (w_prime + coeff * WHQC(PARAM_N, PARAM_OMEGA, PARAM_OMEGA_R, m)) +
        b;
    count_dist_model[m] = expl(
        log_nCk_ld(PARAM_OMEGA, m) +
        log_nCk_ld(PARAM_N - PARAM_OMEGA - 1, PARAM_OMEGA - m - 1) - logC_N1);
  }
}