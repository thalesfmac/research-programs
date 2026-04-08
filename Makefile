# Compilador
# FC := gfortran
FC := h5fc

# Executável final
TARGET := cavityaa_rgf.out

# Diretórios
SRC_DIR := src
APP_DIR := app
OBJ_DIR := build/obj
MOD_DIR := build/mod
BIN_DIR := build/bin

# Flags
FFLAGS := -O2 -Wall -Wextra -J$(MOD_DIR) -I$(MOD_DIR)
LDFLAGS :=
# LDLIBS := -llapack -lblas -lhdf5_fortran -lhdf5
LDLIBS := -llapack -lblas

# Fontes
SRC_FILES := \
	$(SRC_DIR)/precision.f90 \
	$(SRC_DIR)/constants.f90 \
	$(SRC_DIR)/lapack_blas.f90 \
	$(SRC_DIR)/array_io.f90 \
	$(SRC_DIR)/hdf5_io.f90 \
	$(SRC_DIR)/rng_utils.f90 \
	$(SRC_DIR)/matrix_operations.f90 \
	$(SRC_DIR)/peierls_operator.f90 \
	$(SRC_DIR)/lead_green_function.f90 \
	$(SRC_DIR)/transmittance.f90 \
	$(SRC_DIR)/disordered_systems.f90 \
	$(APP_DIR)/cavityaa_rgf.f90

# Objetos
OBJ_FILES := $(patsubst %.f90,$(OBJ_DIR)/%.o,$(notdir $(SRC_FILES)))

# Regra principal
all: dirs $(BIN_DIR)/$(TARGET)

# Linkedição
$(BIN_DIR)/$(TARGET): $(OBJ_FILES)
	$(FC) $(FFLAGS) $(OBJ_FILES) -o $@ $(LDFLAGS) $(LDLIBS)

# Compilação dos arquivos de src
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | dirs
	$(FC) $(FFLAGS) -c $< -o $@

# Compilação do main em app
$(OBJ_DIR)/%.o: $(APP_DIR)/%.f90 | dirs
	$(FC) $(FFLAGS) -c $< -o $@

# Criar diretórios
dirs:
	mkdir -p $(OBJ_DIR) $(MOD_DIR) $(BIN_DIR)

# Limpeza
clean:
	rm -rf build

.PHONY: all clean dirs