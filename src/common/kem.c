/**
 * @file kem.c
 * @brief Implementation of api.h
 */

#include "api.h"
#include "crypto_memset.h"
#include "hqc.h"
#include "parameters.h"
#include "symmetric.h"
#include "vector.h"
#include <stdint.h>
#include <string.h>

/**
 * @brief Generates a keypair for the KEM (Key Encapsulation Mechanism) scheme.
 *
 * This function generates a public/private keypair used for key encapsulation
 * and decapsulation. The encapsulation key (`ek`) is used to encapsulate a
 * shared secret, while the decapsulation key
 * (`dk`) is used to recover it.
 *
 * @param[out] ek_kem Pointer to the output buffer where the encapsulation key
 * will be stored.
 * @param[out] dk_kem Pointer to the output buffer where the decapsulation key
 * will be stored.
 *
 * @return 0 on success.
 *
 * @pre The PRNG **must be seeded** with ::prng_init() before calling this
 * function.
 * @warning This function calls ::prng_get_bytes() to sample `seed_kem`. If the
 * PRNG has not been properly seeded beforehand, the generated keys will be
 * insecure/predictable.
 * @note An example of correct seeding is provided in `main_hqc.c` (see
 * `init_randomness()`), which seeds the PRNG using `syscall(SYS_getrandom,
 * ...)` (32 bytes) by default..
 * @see prng_init, prng_get_bytes, main_hqc.c
 */
int crypto_kem_keypair(uint8_t *ek_kem, uint8_t *dk_kem) {
  uint8_t seed_kem[SEED_BYTES] = {0};
  uint8_t sigma[PARAM_SECURITY_BYTES] = {0};
  uint8_t seed_pke[SEED_BYTES] = {0};
  shake256_xof_ctx ctx_kem;

  uint8_t ek_pke[PUBLIC_KEY_BYTES] = {0};
  uint8_t dk_pke[SEED_BYTES] = {0};

  // Sample seed_kem
  prng_get_bytes(seed_kem, SEED_BYTES);

  // Compute seed_pke and randomness sigma
  xof_init(&ctx_kem, seed_kem, SEED_BYTES);
  xof_get_bytes(&ctx_kem, seed_pke, SEED_BYTES);
  xof_get_bytes(&ctx_kem, sigma, PARAM_SECURITY_BYTES);

  // Compute HQC-PKE keypair
  hqc_pke_keygen(ek_pke, dk_pke, seed_pke);

  // Compute HQC-KEM keypair
  memcpy(ek_kem, ek_pke, PUBLIC_KEY_BYTES);
  memcpy(dk_kem, ek_kem, PUBLIC_KEY_BYTES);
  memcpy(dk_kem + PUBLIC_KEY_BYTES, dk_pke, SEED_BYTES);
  memcpy(dk_kem + PUBLIC_KEY_BYTES + SEED_BYTES, sigma, PARAM_SECURITY_BYTES);
  memcpy(dk_kem + PUBLIC_KEY_BYTES + SEED_BYTES + PARAM_SECURITY_BYTES,
         seed_kem, SEED_BYTES);

  // Zeroize sensitive data
  memset_zero(seed_kem, sizeof seed_kem);
  memset_zero(sigma, sizeof sigma);
  memset_zero(seed_pke, sizeof seed_pke);
  memset_zero(dk_pke, sizeof dk_pke);

  return 0;
}