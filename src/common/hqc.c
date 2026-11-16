/**
 * @file hqc.c
 * @brief Implementation of hqc.h
 */

#include "hqc.h"
#include "crypto_memset.h"
#include "gf2x.h"
#include "parameters.h"
#include "symmetric.h"
#include "vector.h"
#include <stdint.h>
#include <string.h>

/**
 * @brief Generates a key pair for the HQC public-key encryption (PKE) scheme.
 *
 * This function creates a public encryption key (`ek_pke`) and a private
 * decryption key (`dk_pke`) for use in the HQC PKE scheme. The key generation
 * process is seeded with the given `seed` input.
 *
 * @param[out] ek_pke  Pointer to the buffer that will receive the encryption
 * key.
 * @param[out] dk_pke  Pointer to the buffer that will receive the decryption
 * key.
 * @param[in]  seed    Pointer to the seed used to deterministically generate
 * the key pair.
 *
 */
void hqc_pke_keygen(uint8_t *ek_pke, uint8_t *dk_pke, uint8_t *seed) {
  uint8_t keypair_seed[2 * SEED_BYTES] = {0};
  uint8_t *seed_dk = keypair_seed;
  uint8_t *seed_ek = keypair_seed + SEED_BYTES;
  shake256_xof_ctx dk_xof_ctx = {0};
  shake256_xof_ctx ek_xof_ctx = {0};

  static __m256i h[VEC_N_256_NUM_WORDS];
  static __m256i y[VEC_N_256_NUM_WORDS];
  static __m256i x[VEC_N_256_NUM_WORDS];
  static uint64_t s[VEC_N_256_SIZE_64] = {0};
  static __m256i tmp[VEC_N_256_NUM_WORDS] = {0};

#ifdef __STDC_LIB_EXT1__
  memset_s(x, 0, (VEC_N_256_NUM_WORDS) * sizeof(__m256i));
  memset_s(y, 0, (VEC_N_256_NUM_WORDS) * sizeof(__m256i));
  memset_s(h, 0, (VEC_N_256_NUM_WORDS) * sizeof(__m256i));
#else
  memset(x, 0, (VEC_N_256_NUM_WORDS) * sizeof(__m256i));
  memset(y, 0, (VEC_N_256_NUM_WORDS) * sizeof(__m256i));
  memset(h, 0, (VEC_N_256_NUM_WORDS) * sizeof(__m256i));
#endif
  // Derive keypair seeds
  hash_i(keypair_seed, seed);

  // Compute decryption key
  xof_init(&dk_xof_ctx, seed_dk, SEED_BYTES);
  vect_sample_fixed_weight1(&dk_xof_ctx, y, PARAM_OMEGA);
  vect_sample_fixed_weight1(&dk_xof_ctx, x, PARAM_OMEGA);

  // Compute encryption key
  xof_init(&ek_xof_ctx, seed_ek, SEED_BYTES);
  vect_set_random(&ek_xof_ctx, (uint64_t *)h);
  vect_mul(tmp, y, h);
  vect_add(s, (uint64_t *)x, (uint64_t *)tmp, VEC_N_256_SIZE_64);

  // Parse encryption key to string
  memcpy(ek_pke, seed_ek, SEED_BYTES);
  memcpy(ek_pke + SEED_BYTES, s, VEC_N_SIZE_BYTES);

  // Parse decryption key to string
  memcpy(dk_pke, seed_dk, SEED_BYTES);

  // Zeroize sensitive data
  memset_zero(keypair_seed, sizeof keypair_seed);
  memset_zero(x, sizeof x);
  memset_zero(y, sizeof y);
  memset_zero(tmp, sizeof tmp);
  memset_zero(&dk_xof_ctx, sizeof dk_xof_ctx);
}