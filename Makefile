# Disable the default rules
MAKEFLAGS += --no-builtin-rules --no-builtin-variables

# Fortran compiler
COMPILER ?= gfortran

# Compilation profile: debug or release
BUILD ?= release

# User/local configuration
-include config.mk

# Diretórios
SRC_DIR := src
APP_DIR := app
OBJ_DIR := build/obj
MOD_DIR := build/mod
BIN_DIR := build/bin

# Common flags
FFLAGS_COMMON := -I$(MOD_DIR)

# Compiler specific flags
ifeq ($(COMPILER), gfortran)

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

else ifeq ($(COMPILER), ifort)

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

else ifeq ($(COMPILER), ifx)

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
    $(error Unsupported COMPILER='$(COMPILER)'. Use gfortran, ifort, or ifx)
endif

# Final flags
FFLAGS := $(FFLAGS_DEBUG) $(FFLAGS_OPT) $(FFLAGS_WARN) $(MODFLAG) $(FFLAGS_COMMON)

# Modules sources
MOD_FILES := \
    $(SRC_DIR)/precision.f90 \
    $(SRC_DIR)/constants.f90 \
    $(SRC_DIR)/lapack_blas_interface.f90 \
    $(SRC_DIR)/array_io.f90 \
    $(SRC_DIR)/random_number_generator.f90 \
    $(SRC_DIR)/matrix_operations.f90 \
    $(SRC_DIR)/peierls_operator.f90 \
    $(SRC_DIR)/lead_green_function.f90 \
    $(SRC_DIR)/transmittance.f90 \
    $(SRC_DIR)/aubry_andre.f90 \

# Program sources
APP_FILES := $(wildcard $(APP_DIR)/*.f90)

# Object files
OBJ_FILES := $(patsubst %.f90,$(OBJ_DIR)/%.o,$(notdir $(MOD_FILES)))

# Executable files
EXE_FILES := $(patsubst $(APP_DIR)/%.f90,$(BIN_DIR)/%.out,$(APP_FILES))

# Regra principal
all: info dirs $(EXE_FILES)

# Linkedição
$(BIN_DIR)/%.out: $(APP_DIR)/%.f90 $(OBJ_FILES) | dirs
	$(COMPILER) $(FFLAGS) $(OBJ_FILES) $< -o $@ $(LDLIBS)

# Compilação dos arquivos de src
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | dirs
	$(COMPILER) $(FFLAGS) -c $< -o $@

# Compilação do main em app
$(OBJ_DIR)/%.o: $(APP_DIR)/%.f90 | dirs
	$(COMPILER) $(FFLAGS) -c $< -o $@

# Criar diretórios
dirs:
	mkdir -p $(OBJ_DIR) $(MOD_DIR) $(BIN_DIR)

info:
	@echo "Configuration:"
	@echo "  COMPILER = $(COMPILER)"
	@echo "  BUILD    = $(BUILD)"
	@echo "  SRC_DIR  = $(SRC_DIR)"
	@echo "  APP_DIR  = $(APP_DIR)"
	@echo "  OBJ_DIR  = $(OBJ_DIR)"
	@echo "  MOD_DIR  = $(MOD_DIR)"
	@echo "  BIN_DIR  = $(BIN_DIR)"
	@echo "  FFLAGS   = $(FFLAGS)"
	@echo "  LDLIBS   = $(LDLIBS)"

# Limpeza
clean:
	rm -rf build

.PHONY: all clean dirs info
