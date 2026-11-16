/**
 * @file parsing.c
 * @brief Functions to parse secret key, public key and ciphertext of the HQC
 * scheme
 */

#include "parsing.h"
#include "crypto_memset.h"
#include "symmetric.h"
#include "vector.h"
#include <stdint.h>
#include <string.h>

/**
 * @brief Deserializes a decryption key into its internal vectorized form.
 *
 * @param[out] y     Pointer to the output buffer where the internal vectorized
 * key will be stored.
 * @param[in]  dk_pke   Pointer to the serialized decryption key.
 */
void hqc_dk_pke_from_string(__m256i *y, const uint8_t *dk_pke) {
  shake256_xof_ctx dk_xof_ctx = {0};
  xof_init(&dk_xof_ctx, dk_pke, SEED_BYTES);
  vect_sample_fixed_weight1(&dk_xof_ctx, y, PARAM_OMEGA);

  // Zeroize sensitive data
  memset_zero(&dk_xof_ctx, sizeof dk_xof_ctx);
}

/**
 * @brief Deserializes an encryption key into its internal representation.
 *
 * @param[out] h        Pointer to the output buffer for `h` the first internal
 * component of the key.
 * @param[out] s        Pointer to the output buffer for `s` the second internal
 * component of the key.
 * @param[in]  ek_pke   Pointer to the serialized encryption key.
 */
void hqc_ek_pke_from_string(uint64_t *h, uint64_t *s, const uint8_t *ek_pke) {
  shake256_xof_ctx ek_xof_ctx = {0};

  xof_init(&ek_xof_ctx, ek_pke, SEED_BYTES);
  vect_set_random(&ek_xof_ctx, h);

  memcpy(s, ek_pke + SEED_BYTES, VEC_N_SIZE_BYTES);
}