
CC = gcc

BUILD_DIR = bin/build

# Executables
TARGET_SIM = bin/simulation
TARGET_CCA = bin/cca
TARGET_KCA = bin/kca

# Common source files
COMMON_SRC = \
    src/common/vector.c \
    src/common/gf2x.c \
    src/common/parsing.c \
    src/common/hqc.c \
    src/common/kem.c \
    src/common/symmetric.c \
    src/common/crypto_memset.c \
    src/common/poly.c \
    src/common/spectrum.c \
    src/common/err_attack.c \
    src/common/model.c \
    src/common/annealing_soft_spectrum.c \
    src/common/guess_key.c \
    lib/fips202/fips202.c

# Specific source files

SIM_SRC = src/simu/main_hqc.c
CCA_SRC = src/cca/main_cca.c
KCA_SRC = src/kca/main_kca.c

# Object files
COMMON_OBJ = $(patsubst %.c,$(BUILD_DIR)/%.o,$(COMMON_SRC))
SIM_OBJ = $(patsubst %.c,$(BUILD_DIR)/%.o,$(SIM_SRC))
CCA_OBJ = $(patsubst %.c,$(BUILD_DIR)/%.o,$(CCA_SRC))
KCA_OBJ = $(patsubst %.c,$(BUILD_DIR)/%.o,$(KCA_SRC))

OBJ = $(COMMON_OBJ) $(SIM_OBJ) $(CCA_OBJ) $(KCA_OBJ)
DEP = $(OBJ:.o=.d)

# Compilation flags
CPPFLAGS = -Isrc/simu -Isrc/common -Isrc/cca -Isrc/kca -Ilib/fips202

CFLAGS = -O3 -funroll-all-loops -flto \
         -pedantic -Wall -Wextra -mavx2 \
         -mpclmul -fopenmp -MMD -MP

LDFLAGS =
LDLIBS = -lm

.PHONY: all simulation cca kca clean

# Build everything
all: simulation cca kca

# Individual targets
simulation: $(TARGET_SIM)
cca: $(TARGET_CCA)
kca: $(TARGET_KCA)

# Link simulation executable
$(TARGET_SIM): $(COMMON_OBJ) $(SIM_OBJ)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

# Link CCA executable
$(TARGET_CCA): $(COMMON_OBJ) $(CCA_OBJ)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

# Link KCA executable
$(TARGET_KCA): $(COMMON_OBJ) $(KCA_OBJ)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

# Compile source files
$(BUILD_DIR)/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# Automatically track header dependencies
-include $(DEP)

# Clean
clean:
	rm -rf $(BUILD_DIR) $(TARGET_SIM) $(TARGET_CCA) $(TARGET_KCA)
