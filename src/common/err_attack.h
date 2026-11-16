#pragma once

#include <immintrin.h>
#include <stdbool.h>

void init_synd_attack(uint16_t *x, int length);
void set_synd_attack_logging(bool enabled);
void log_synd_attack(double *syndromes, uint16_t *e, double syndrome_weight,
                     long int nb_dec, double sigma);
void display_synd_attack(uint8_t *entropy, double sigma);
void recover_x_y(__m256i *y, __m256i *x, unsigned char *sk, unsigned char *pk);
void compute_syndrome(uint64_t *s, __m256i *x, __m256i *y, __m256i *e,
                      __m256i *r1, __m256i *r2);
void extract_bit_positions(const uint64_t *input, long int input_len,
                           uint16_t *output, long int max_output_len);
