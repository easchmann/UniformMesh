# Build and test one variant:   make check DIM=3d LAYOUT=aos HALO=4
# Test every variant:           make check-all
# Show the tools in use:        make info


DIM    ?= 2d
LAYOUT ?= soa
HALO   ?= 1

# make's built-in defaults for CC and FC are cc and f77, so ?= would never apply
ifeq ($(origin CC),default)
CC = gcc
endif
ifeq ($(origin FC),default)
FC = gfortran
endif
PYTHON ?= python3
RUN    ?=
CFLAGS ?= -std=c11 -Wall -Wextra -O2
FFLAGS ?= -Wall -O2 -fcheck=all


INI = macros/general.ini macros/dim_$(DIM).ini macros/layout_$(LAYOUT).ini \
      $(if $(filter-out 1,$(HALO)),macros/halo_$(HALO).ini)
FDEFS = $(if $(filter soa,$(LAYOUT)),-DMESH_SOA)

.PHONY: all check check-all clean info FORCE

all: test_mesh test_mesh_f

# regenerated on every run, but only rewritten (and rebuilt against) when it changes
mesh_config.h: FORCE
	@$(PYTHON) ini_to_h.py $(INI) -o $@.tmp
	@if cmp -s $@.tmp $@; then rm $@.tmp; else mv $@.tmp $@; fi

mesh.o: mesh.c mesh.h mesh_config.h
	$(CC) $(CFLAGS) -c mesh.c -o $@

mesh_bind.o: mesh_bind.c mesh.h mesh_config.h
	$(CC) $(CFLAGS) -c mesh_bind.c -o $@

mesh_f.o mesh_f.mod: mesh_f.F90 mesh_config.h
	$(FC) $(FFLAGS) -c mesh_f.F90 -o mesh_f.o

test_mesh_f.o: test_mesh_f.F90 mesh_f.mod mesh_config.h
	$(FC) $(FFLAGS) $(FDEFS) -c test_mesh_f.F90 -o $@

test_mesh: mesh.o test_mesh.c
	$(CC) $(CFLAGS) mesh.o test_mesh.c -lm -o $@

test_mesh_f: mesh.o mesh_bind.o mesh_f.o test_mesh_f.o
	$(FC) $(FFLAGS) $^ -o $@

check: test_mesh test_mesh_f
	$(RUN) ./test_mesh out_c.bin
	$(RUN) ./test_mesh_f

# which compilers and python are picked up (useful after `module load`)
info:
	@echo "CC     = $(CC)";     $(CC) --version | head -n 1
	@echo "FC     = $(FC)";     $(FC) --version | head -n 1
	@echo "PYTHON = $(PYTHON)"; $(PYTHON) --version
	@echo "RUN    = $(RUN)"
	@echo "variant: DIM=$(DIM) LAYOUT=$(LAYOUT) HALO=$(HALO)"

check-all:
	@for h in 1 4 8; do for d in 2d 3d; do for l in soa aos; do \
	    $(MAKE) -s clean-obj; \
	    $(MAKE) -s check DIM=$$d LAYOUT=$$l HALO=$$h > check.log 2>&1 \
	        && echo "ok   halo=$$h $$d $$l" || { echo "FAIL halo=$$h $$d $$l"; cat check.log; exit 1; }; \
	done; done; done; rm -f check.log; $(MAKE) -s clean-obj

.PHONY: clean-obj
clean-obj:
	rm -f *.o *.mod test_mesh test_mesh_f

clean: clean-obj
	rm -f mesh_config.h mesh_config.h.tmp *.bin check.log
