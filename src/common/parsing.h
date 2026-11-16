/**
 * @file parsing.h
 * @brief Header file for parsing.c
 */

#ifndef HQC_PARSING_H
#define HQC_PARSING_H

#include "data_structures.h"
#include "parameters.h"
#include <immintrin.h>
#include <stdint.h>

void hqc_dk_pke_from_string(__m256i *y, const uint8_t *dk_pke);
void hqc_ek_pke_from_string(uint64_t *h, uint64_t *s, const uint8_t *ek_pke);

#endif // HQC_PARSING_H
