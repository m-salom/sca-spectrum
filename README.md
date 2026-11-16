# HQC secret-key reconstruction from soft spectrum

Research code for reconstructing the secret sparse polynomial of **HQC-1** from a soft spectrum. The repository supports simulated error-weight observations and two sets of previously estimated side-channel error weights. It compares two reconstruction methods:

- **SA** — simulated annealing (`src/common/annealing_soft_spectrum.c`)
- **DFS** — depth-first search (`src/common/guess_key.c`)

The side-channel experiments do not acquire traces or perform side-channel profiling: they start from precomputed estimated error weights and the corresponding `r1` vectors.

> **Paper:** Magali Salom, Angélique Lopez, Valentin Vasseur & Gabriel Zaid: Recovering HQC Secret Key from Distance Spectrum using Side-Channel Attack.

## Overview

| Dataset/source | Executable | Number of observations (`-N`) |
| --- | --- | --- |
| HQC simulation | `bin/simulation` | Number of simulated error weights |
| CCA side-channel data | `bin/cca` | Number of selected trace-derived weight estimates (max. 410,000) |
| KCA side-channel data | `bin/kca` | Number of selected trace-derived weight estimates (max. 1,394,000) |

Each executable constructs a soft spectrum and runs **SA**, **DFS**, or **BOTH** on that spectrum. The code uses HQC-1 parameters (`PARAM_N = 17669`, `PARAM_OMEGA = 66`, `PARAM_OMEGA_R = 75`, `PARAM_OMEGA_E = 75`). The same source tree contains shared HQC routines, model/spectrum construction and both solvers.

## Requirements

- Linux with GCC, `make` and the C math library
- x86-64 CPU supporting **AVX2** and **PCLMULQDQ** (the Makefile enables `-mavx2 -mpclmul`)
- GCC OpenMP support (`-fopenmp`)
- Experimental binary data for the CCA and KCA executables, placed at the paths below

## Build

From the repository root:

```bash
make                  # Build all three executables
make simulation       # Build only simulation
make cca              # Build only CCA
make kca              # Build only KCA
make clean            # Remove binaries and intermediate object files
```

Expected output executables: `bin/simulation`, `bin/cca`, `bin/kca`.

## Experimental data

The CCA and KCA executables expect these files, relative to the **repository root**:

```text
data/
├── cca/
│   ├── estimated_error_weights_cca.bin
│   └── r1_vectors_cca.bin
└── kca/
    ├── estimated_error_weights_kca.bin
    └── r1_vectors_kca.bin
```

The estimated weights are read as consecutive native C `float` values. The `r1` files are read as consecutive packed vectors, each containing `VEC_N_SIZE_64` native C `uint64_t` words. For each source, entry *i* in the weight file must correspond to entry *i* in the `r1` file. The programs select `N` distinct indices randomly (without replacement), sort them, then read the matching entries sequentially.

### Experimental Data

The estimated error weights used in the CCA and KCA experiments are included in this GitHub repository.

Due to their size, the corresponding $r_1$ vectors are hosted separately on Zenodo:

**Dataset DOI:** https://doi.org/10.5281/zenodo.24085606

To reproduce the experiments, download the files from Zenodo and place them in the following directories:

```text
data/
├── cca/
│   ├── estimated_error_weights_cca.bin
│   └── r1_vectors_cca.bin
└── kca/
    ├── estimated_error_weights_kca.bin
    └── r1_vectors_kca.bin
```

The estimated error-weight files are already provided in the repository. Only the $r_1$ vector files need to be downloaded.

The dataset contains intermediate experimental data, not raw side-channel traces. The provided files are sufficient to reproduce the soft-spectrum construction and subsequent key-recovery experiments.

**Important:** The experimental programs regenerate the key associated with the data using a fixed seed in their `main_cca.c` / `main_kca.c`; changing the HQC implementation or the seed can break this correspondence. The code assumes compatible binary layout and endianness with the machine used to produce the files.

## Usage

Run commands from the repository root so relative paths resolve correctly.

### Simulation

```bash
./bin/simulation -m SA -N 30000
./bin/simulation -m DFS -N 190000
./bin/simulation -m BOTH -N 200000
./bin/simulation -m SA -N 30000 -s 2.0
```

`-s` / `--sigma` controls the additive measurement noise used by the simulation (default `0`).

### CCA

```bash
./bin/cca -m SA -N 50000
./bin/cca -m DFS -N 100000
./bin/cca -m BOTH -N 100000
```

The CCA dataset contains 410,000 candidate observations; `-N` must not exceed that value.

### KCA

```bash
./bin/kca -m SA -N 1000000
./bin/kca -m DFS -N 1000000
./bin/kca -m BOTH -N 1000000
```

The KCA dataset contains 1,394,000 candidate observations; `-N` must not exceed that value.

**Depth-First Search (DFS).** Although the DFS implementation is available for KCA experiments, no successful key recovery was observed in 200 trials under the evaluated experimental conditions. Simulated annealing is therefore the primary reconstruction method used for KCA.

### Command-line options

| Option | Meaning |
| --- | --- |
| `-h`, `--help` | Print help |
| `-m`, `--mode` | `SA`, `DFS` or `BOTH` (case-sensitive) |
| `-N`, `--samples` | Number of observations |
| `-s`, `--sigma` | Additional Gaussian noise **(simulation only)** |
| `-i`, `--iterations` | SA iterations per restart |
| `-r`, `--restarts` | SA restarts |
| `-e`, `--exploration-restarts` | SA exploration restarts |
| `-Ts`, `--start-temp` | SA initial temperature |
| `-Tf`, `--end-temp` | SA final temperature |
| `-k`, `--init-top-k` | SA initialization candidate limit |
| `-T`, `--swap-trials` | SA guided-move trials |
| `-P`, `--swap-proba` | SA random-swap probability |

## Interpretation of results

To reproduce success-rate curves from the paper, the intended experimental procedure is to repeat each `(source, N, algorithm)` combination many times, record the reconstruction result, and report the proportion of successful trials. 

## Repository layout

```text
.
├── Makefile
├── lib/fips202/                  # FIPS 202 implementation
├── src/
│   ├── simu/main_hqc.c          # Simulated error-weight observations
│   ├── cca/main_cca.c           # CCA weight estimates -> soft spectrum
│   ├── kca/main_kca.c           # KCA weight estimates -> soft spectrum
│   └── common/
│       ├── annealing_soft_spectrum.c  # SA key reconstruction
│       ├── guess_key.c                # DFS key reconstruction
│       ├── model.c                    # Statistical error-weight model
│       ├── spectrum.c                 # Spectrum and posterior computation
│       ├── err_attack.c               # Observation aggregation
│       └── ...                        # HQC and utility routines
└── data/                       # External experimental data
```

## Code Attribution and Licensing

This project builds upon existing implementations of HQC and a depth-first search key-reconstruction algorithm.

- **HQC implementation:** The HQC cryptographic operations are based on the [HQC implementation](https://gitlab.com/pqc-hqc/hqc). The original implementation is released into the public domain and incorporates FIPS202 code derived from the PQClean project.

- **Depth-First Search (DFS):** The key-reconstruction implementation is adapted from the code accompanying https://link.springer.com/chapter/10.1007/978-3-032-22695-2_4, available at https://github.com/m-salom/key_rec_spectrum. The original code is distributed under the MIT License.

The simulated annealing implementation, soft-spectrum extensions, and integration of side-channel-derived error-weight estimates were developed for this project.

The original copyright and license notices of the third-party components are retained. The original HQC implementation is declared to be in the public domain, while the reused DFS implementation is subject to the MIT License.

The original works should be cited when using the corresponding implementations.