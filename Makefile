# HDF5 wrapper
FC := h5fc

# Fortran compiler
HDF5_FC ?= gfortran
export HDF5_FC

# Compilation profile: debug or release
BUILD ?= release

# Executável final
TARGET := cavityaa_rgf.out

# Diretórios
SRC_DIR := src
APP_DIR := app
OBJ_DIR := build/obj
MOD_DIR := build/mod
BIN_DIR := build/bin

# Common flags
FFLAGS_COMMON := -I$(MOD_DIR)

# Compiler specific flags
ifeq ($(HDF5_FC), gfortran)

	MODFLAG := -J$(MOD_DIR)
	LDLIBS := -llapack -lblas

	ifeq ($(BUILD), release)
		FFLAGS_OPT := -O2
		FFLAGS_WARN :=
		FFLAGS_DEBUG :=
	else ifeq ($(BUILD), debug)
		FFLAGS_OPT := -Og
		FFLAGS_WARN := -Wall -Wextra
		FFLAGS_DEBUG := -g -fcheck=all -fbacktrace
	else
		$(error Unsupported BUILD='$(BUILD)'. Use debug or release)
	endif

else ifeq ($(HDF5_FC), ifort)

	MODFLAG := -module $(MOD_DIR)
	LDLIBS := -qmkl

	ifeq ($(BUILD), debug)
		FFLAGS_OPT := -O0
		FFLAGS_WARN := -warn all
		FFLAGS_DEBUG := -g -traceback -check all -diag-disable=10448
	else ifeq ($(BUILD), release)
		FFLAGS_OPT := -O2
		FFLAGS_WARN := -warn all
		FFLAGS_DEBUG := -diag-disable=10448
	else
		$(error Unsupported BUILD='$(BUILD)'. Use debug or release)
	endif

else ifeq ($(HDF5_FC), ifx)

	MODFLAG := -module $(MOD_DIR)
	LDLIBS := -qmkl

	ifeq ($(BUILD), debug)
		FFLAGS_OPT := -O0
		FFLAGS_WARN := -warn all
		FFLAGS_DEBUG := -g -traceback -check all
	else ifeq ($(BUILD), release)
		FFLAGS_OPT := -O2
		FFLAGS_WARN := -warn all
		FFLAGS_DEBUG :=
	else
		$(error Unsupported BUILD='$(BUILD)'. Use debug or release)
	endif

else
	$(error Unsupported HDF5_FC='$(HDF5_FC)'. Use gfortran, ifort, or ifx)
endif

# Final flags
FFLAGS := $(FFLAGS_DEBUG) $(FFLAGS_OPT) $(FFLAGS_WARN) $(MODFLAG) $(FFLAGS_COMMON)

# Source files
SRC_FILES := \
	$(SRC_DIR)/precision.f90 \
	$(SRC_DIR)/constants.f90 \
	$(SRC_DIR)/lapack_blas_interface.f90 \
	$(SRC_DIR)/array_io.f90 \
	$(SRC_DIR)/hdf5_io.f90 \
	$(SRC_DIR)/rng_utils.f90 \
	$(SRC_DIR)/matrix_operations.f90 \
	$(SRC_DIR)/peierls_operator.f90 \
	$(SRC_DIR)/lead_green_function.f90 \
	$(SRC_DIR)/transmittance.f90 \
	$(SRC_DIR)/disordered_systems.f90 \
	$(APP_DIR)/cavityaa_rgf.f90

# Object files
OBJ_FILES := $(patsubst %.f90, $(OBJ_DIR)/%.o, $(notdir $(SRC_FILES)))

# Regra principal
all: dirs $(BIN_DIR)/$(TARGET)

# Linkedição
$(BIN_DIR)/$(TARGET): $(OBJ_FILES)
	$(FC) $(FFLAGS) $(OBJ_FILES) -o $@ $(LDLIBS)

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