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
TEST_DIR := tests

# Flags
FFLAGS := -g -fcheck=all -fbacktrace -Og -Wall -Wextra -J$(MOD_DIR) -I$(MOD_DIR)
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

TEST_SRC := $(wildcard $(TEST_DIR)/*.f90)
TEST_EXE := $(patsubst $(TEST_DIR)/%.f90,$(BIN_DIR)/%.test,$(TEST_SRC))

# Regra principal
all: dirs $(BIN_DIR)/$(TARGET)

test: $(TEST_EXE)
	@for ex in $(TEST_EXE); do \
		echo "Running $$ex..."; \
		./$$ex; \
	done

# Linkedição
$(BIN_DIR)/$(TARGET): $(OBJ_FILES)
	$(FC) $(FFLAGS) $(OBJ_FILES) -o $@ $(LDFLAGS) $(LDLIBS)

$(BIN_DIR)/%.test: $(OBJ_DIR)/$(TEST_DIR)/%.o $(filter-out $(OBJ_DIR)/cavityaa_rgf.o, $(OBJ_FILES))
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

# Compilação dos arquivos de src
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | dirs
	$(FC) $(FFLAGS) -c $< -o $@

# Compilação do main em app
$(OBJ_DIR)/%.o: $(APP_DIR)/%.f90 | dirs
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/$(TEST_DIR)/%.o: $(TEST_DIR)/%.f90 | dirs
	$(FC) $(FFLAGS) -c $< -o $@

# Criar diretórios
dirs:
	mkdir -p $(OBJ_DIR) $(MOD_DIR) $(BIN_DIR) $(OBJ_DIR)/$(TEST_DIR)

# Limpeza
clean:
	rm -rf build

.PHONY: all clean dirs test