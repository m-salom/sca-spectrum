#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "model.h"
#include "spectrum.h"

/* min(|i - j|, r - |i - j|)
   assumes 0 <= i,j < r (no need to reduce abs(i - j) modulo r) */
int distance(int i, int j, int r) {
  return ((r - abs(r - 2 * abs(i - j))) / 2) - 1;
}

void spectrum(int8_t *S, uint16_t *index, int r, int d) {
  memset(S, 0, (r / 2) * sizeof(int8_t));

  for (int i = 0; i < d; ++i) {
    for (int j = 0; j < i; ++j) {
      S[distance(index[i], index[j], r)]++;
    }
  }
}

static double normal_pdf(double x, double mu, double sigma) {
  return (1.0 / (sigma * sqrt(2 * M_PI))) *
         exp(-0.5 * pow((x - mu) / sigma, 2));
}

// Compute the probability of a distance to belong to each class (multiplicity)
static void bayesian_inference(double x, double *means, double *std_dev,
                               long double *count_dist, int num_classes,
                               double *posteriors) {
  double denominator = 0.0;
  double f_k[num_classes];

  for (int i = 0; i < num_classes; i++) {
    f_k[i] = normal_pdf(x, means[i], sqrt(std_dev[i]));
    denominator += f_k[i] * count_dist[i];
  }

  for (int i = 0; i < num_classes; i++)
    posteriors[i] = (f_k[i] * count_dist[i]) / denominator;
}

/*
    Compute spectrum using soft information:
    Each line i corresponds to a distance
    Each column j corresponds to the probability of the multiplicity j
*/
static void spectrum_soft_info(double S[PARAM_N / 2][max_mult_y], double *means,
                               double *avg_dist, double *std_devs,
                               long double *count_dist, int r) {
  double posteriors[max_mult_y];

  for (int i = 0; i < r / 2; i++) {
    bayesian_inference(means[i], avg_dist, std_devs, count_dist, max_mult_y,
                       posteriors);
    for (int j = 0; j < max_mult_y; j++) {
      S[i][j] = posteriors[j];
    }
  }
}

void get_soft_spectrum(double S_soft[PARAM_N / 2][max_mult_y],
                       double *syndromes, long int nb_dec, double *means,
                       double sigma, double a, double b, double sigma0) {
  double avg_dist_model[max_mult_y];
  double std_dist_model[max_mult_y];
  long double count_dist_model[max_mult_y];

  shifting_values(syndromes, nb_dec, a, b);
  compute_model(avg_dist_model, std_dist_model, count_dist_model, sigma, nb_dec,
                a, b, sigma0);
  spectrum_soft_info(S_soft, means, avg_dist_model, std_dist_model,
                     count_dist_model, PARAM_N);
}