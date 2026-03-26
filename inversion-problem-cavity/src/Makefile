SHELL := /bin/bash
.ONESHELL:

# Compilador
FC := ifort

# Flags
FFLAGS := -diag-disable=10448 -qmkl -O2

# Todos os fontes .f90 do diretório
SRC := $(wildcard *.f90)

# Executáveis correspondentes
EXE := $(SRC:.f90=.x)

all: $(EXE)

%.x: %.f90
	@module purge
	@module load intel/2024.0.0
	@echo "Compilando $< -> $@"
	@$(FC) $(FFLAGS) $< -o $@
	@if [ $$? -eq 0 ]; then
		echo "✓ Sucesso: $@ criado."
	else
		echo "✗ Erro ao compilar $<."
	fi

clean:
	rm -f $(EXE) *.mod

.PHONY: all clean
