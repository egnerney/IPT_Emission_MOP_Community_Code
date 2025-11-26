#==============================================================================
# Makefile for IPT UV Emission Raytracer
#
# Supports macOS and Linux with automatic HDF5 detection
# Run 'make help' for usage information
#==============================================================================

# Compiler
FC = gfortran

# Base flags
FFLAGS_BASE = -cpp -ffree-line-length-none -fimplicit-none

# Optimization flags
FFLAGS_OPT = -O3 -march=native -funroll-loops -ftree-vectorize

# Warning flags
FFLAGS_WARN = -Wall -Wno-unused-dummy-argument

# OpenMP flag
FFLAGS_OMP = -fopenmp

# Debug flags
FFLAGS_DEBUG = -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow

# Profile flags
FFLAGS_PROFILE = -pg -g

#------------------------------------------------------------------------------
# HDF5 Configuration - Auto-detect installation
#------------------------------------------------------------------------------

# Check for user-provided HDF5_DIR
ifdef HDF5_DIR
    HDF5_INC = -I$(HDF5_DIR)/include
    HDF5_LIB = -L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5
    HDF5_RPATH = -Wl,-rpath,$(HDF5_DIR)/lib
else
    # Auto-detect OS
    UNAME_S := $(shell uname -s)
    
    ifeq ($(UNAME_S),Darwin)
        # macOS - check Homebrew locations
        ifneq ($(wildcard /usr/local/opt/hdf5/include),)
            HDF5_DIR_AUTO = /usr/local/opt/hdf5
        else ifneq ($(wildcard /opt/homebrew/opt/hdf5/include),)
            HDF5_DIR_AUTO = /opt/homebrew/opt/hdf5
        else
            $(warning HDF5 not found in standard Homebrew locations)
            $(warning Install with: brew install hdf5)
            HDF5_DIR_AUTO = /usr/local/opt/hdf5
        endif
        HDF5_INC = -I$(HDF5_DIR_AUTO)/include
        HDF5_LIB = -L$(HDF5_DIR_AUTO)/lib -lhdf5_fortran -lhdf5
        HDF5_RPATH = -Wl,-rpath,$(HDF5_DIR_AUTO)/lib
    else ifeq ($(UNAME_S),Linux)
        # Linux - check common locations
        ifneq ($(wildcard /usr/include/hdf5/serial),)
            HDF5_INC = -I/usr/include/hdf5/serial
            HDF5_LIB = -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_fortran -lhdf5
            HDF5_RPATH =
        else ifneq ($(wildcard /usr/include/hdf5.mod),)
            HDF5_INC = -I/usr/include
            HDF5_LIB = -lhdf5_fortran -lhdf5
            HDF5_RPATH =
        else
            $(warning HDF5 not found in standard Linux locations)
            $(warning Install with: sudo apt install libhdf5-dev)
            HDF5_INC = -I/usr/include
            HDF5_LIB = -lhdf5_fortran -lhdf5
            HDF5_RPATH =
        endif
    else
        $(warning Unknown OS: $(UNAME_S))
        HDF5_INC = -I/usr/local/include
        HDF5_LIB = -lhdf5_fortran -lhdf5
        HDF5_RPATH =
    endif
endif

#------------------------------------------------------------------------------
# Combined Flags
#------------------------------------------------------------------------------

# Default: optimized parallel build
FFLAGS = $(FFLAGS_BASE) $(FFLAGS_OPT) $(FFLAGS_OMP) $(FFLAGS_WARN) $(HDF5_INC)
LDFLAGS = $(HDF5_LIB) $(HDF5_RPATH)

#------------------------------------------------------------------------------
# Source Files
#------------------------------------------------------------------------------

MODULE_SRC = IPT_emiss_MOP_community_code.f90
MAIN_SRC = basic_example_uv_integration_emission_model_use_tables.f90

MODULE_OBJ = $(MODULE_SRC:.f90=.o)
MAIN_OBJ = $(MAIN_SRC:.f90=.o)

TARGET = uv_emission_example

#------------------------------------------------------------------------------
# Build Rules
#------------------------------------------------------------------------------

.PHONY: all serial parallel debug profile clean distclean help info test-hdf5

# Default target - parallel optimized build
all: parallel

# Serial build (no OpenMP)
serial: FFLAGS = $(FFLAGS_BASE) $(FFLAGS_OPT) $(FFLAGS_WARN) $(HDF5_INC)
serial: $(TARGET)
	@echo ""
	@echo "Serial build complete: $(TARGET)"

# Parallel build (with OpenMP)
parallel: FFLAGS = $(FFLAGS_BASE) $(FFLAGS_OPT) $(FFLAGS_OMP) $(FFLAGS_WARN) $(HDF5_INC)
parallel: LDFLAGS += $(FFLAGS_OMP)
parallel: $(TARGET)
	@echo ""
	@echo "Build complete: $(TARGET)"
	@echo ""
	@echo "To run the example:"
	@echo "  ./$(TARGET)"
	@echo ""
	@echo "Make sure the data files are in the expected locations:"
	@echo "  ../3D_Torus_Model/jovian_plasma_interpolated_381x381x231.h5"
	@echo "  ../Emiss_tables/CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5"
	@echo "  ../Emiss_tables/CHIANTI_11.0.2_emiss_tables_double_maxwellian_24x10x24x12.h5"

# Debug build
debug: FFLAGS = $(FFLAGS_BASE) $(FFLAGS_DEBUG) $(FFLAGS_WARN) $(HDF5_INC)
debug: $(TARGET)
	@echo ""
	@echo "Debug build complete: $(TARGET)"

# Profile build
profile: FFLAGS = $(FFLAGS_BASE) $(FFLAGS_OPT) $(FFLAGS_PROFILE) $(FFLAGS_WARN) $(HDF5_INC)
profile: LDFLAGS += $(FFLAGS_PROFILE)
profile: $(TARGET)
	@echo ""
	@echo "Profile build complete: $(TARGET)"
	@echo "Run with: ./$(TARGET)"
	@echo "Then analyze with: gprof $(TARGET) gmon.out"

#------------------------------------------------------------------------------
# Compilation Rules
#------------------------------------------------------------------------------

$(TARGET): $(MODULE_OBJ) $(MAIN_OBJ)
	$(FC) $(FFLAGS) -o $@ $(MODULE_OBJ) $(MAIN_OBJ) $(LDFLAGS)

$(MODULE_OBJ): $(MODULE_SRC)
	$(FC) $(FFLAGS) -c $< -o $@

$(MAIN_OBJ): $(MAIN_SRC) $(MODULE_OBJ)
	$(FC) $(FFLAGS) -c $< -o $@

#------------------------------------------------------------------------------
# Utility Targets
#------------------------------------------------------------------------------

clean:
	rm -f $(TARGET) *.o *.mod *.smod
	rm -f *.dat *.gp
	@echo "Cleaned build files"

distclean: clean
	rm -f *.png gmon.out
	@echo "Cleaned all generated files"

info:
	@echo "=============================================="
	@echo "IPT UV Emission Raytracer - Build Information"
	@echo "=============================================="
	@echo ""
	@echo "Compiler: $(FC)"
	@echo "OS: $(shell uname -s)"
	@echo ""
	@echo "HDF5 Configuration:"
ifdef HDF5_DIR
	@echo "  HDF5_DIR (user): $(HDF5_DIR)"
else
	@echo "  HDF5_DIR (auto): $(HDF5_DIR_AUTO)"
endif
	@echo "  Include: $(HDF5_INC)"
	@echo "  Library: $(HDF5_LIB)"
	@echo ""
	@echo "Flags:"
	@echo "  FFLAGS: $(FFLAGS)"
	@echo "  LDFLAGS: $(LDFLAGS)"

test-hdf5:
	@echo "Testing HDF5 installation..."
	@echo 'program test_hdf5; use hdf5; integer :: err; call h5open_f(err); print *, "HDF5 OK"; end' > /tmp/test_hdf5.f90
	@$(FC) $(HDF5_INC) -o /tmp/test_hdf5 /tmp/test_hdf5.f90 $(HDF5_LIB) $(HDF5_RPATH) 2>/dev/null && \
		/tmp/test_hdf5 && rm -f /tmp/test_hdf5 /tmp/test_hdf5.f90 || \
		(echo "HDF5 test failed - check installation"; rm -f /tmp/test_hdf5 /tmp/test_hdf5.f90; exit 1)

help:
	@echo "=============================================="
	@echo "IPT UV Emission Raytracer - Makefile Help"
	@echo "=============================================="
	@echo ""
	@echo "Build targets:"
	@echo "  make          - Build optimized parallel version (default)"
	@echo "  make serial   - Build optimized serial version (no OpenMP)"
	@echo "  make parallel - Build optimized parallel version"
	@echo "  make debug    - Build with debug flags and runtime checks"
	@echo "  make profile  - Build with profiling enabled"
	@echo ""
	@echo "Utility targets:"
	@echo "  make clean     - Remove object files and modules"
	@echo "  make distclean - Remove all generated files"
	@echo "  make info      - Display build configuration"
	@echo "  make test-hdf5 - Test HDF5 installation"
	@echo "  make help      - Display this help message"
	@echo ""
	@echo "Environment variables:"
	@echo "  HDF5_DIR - Override HDF5 installation path"
	@echo "             Example: make HDF5_DIR=/opt/hdf5"
	@echo ""
	@echo "Requirements:"
	@echo "  - gfortran (GCC Fortran compiler)"
	@echo "  - HDF5 with Fortran bindings"
	@echo "  - gnuplot (for plotting)"
	@echo ""
	@echo "macOS installation:"
	@echo "  brew install gcc hdf5 gnuplot"
	@echo ""
	@echo "Linux installation (Debian/Ubuntu):"
	@echo "  sudo apt install gfortran libhdf5-dev gnuplot"
