/**
 * @file    data_structures.h
 * @brief   HQC PKE and KEM ciphertext data structures.
 */

#ifndef HQC_DATA_STRUCTURES_H
#define HQC_DATA_STRUCTURES_H

#include "parameters.h"
#include <stdint.h>
#include <x86intrin.h>

/**
 * @brief Public-key encryption ciphertext.
 *
 * - u: length VEC_N_256_SIZE_64 words
 * - v: length VEC_N1N2_256_SIZE_64 words
 */
typedef struct {
  uint64_t u[VEC_N_256_SIZE_64];    /**< first vector half */
  uint64_t v[VEC_N1N2_256_SIZE_64]; /**< second vector half */
} ciphertext_pke_t;

/**
 * @brief Key-encapsulation mechanism ciphertext.
 *
 * Wraps a PKE ciphertext along with the salt used in the KEM:
 * - c_pke: the ciphertext for PKE
 * - salt:  additional randomness (SALT_BYTES bytes)
 */
typedef struct {
  ciphertext_pke_t c_pke;   /**< embedded PKE ciphertext */
  uint8_t salt[SALT_BYTES]; /**< per-encapsulation salt */
} ciphertext_kem_t;

#endif // HQC_DATA_STRUCTURES_H
