#==============================================================================
# Makefile for IPT Emission Movie Frame Generator (MPI-Parallelized)
#
# Io Plasma Torus Hybrid Optical+UV Emission Movie Generator
# MOP (Magnetospheres of Outer Planets) Community Code
#
# USAGE:
#   make              - Build the MPI-parallelized frame generator
#   make debug        - Build with debugging enabled
#   make clean        - Remove build artifacts
#   make help         - Show available targets
#
# PREREQUISITES:
#   - MPI Fortran compiler (mpifort)
#   - HDF5 library with Fortran bindings
#   - IPT_emiss_MOP_community_code.f90 module in this directory
#
# PLATFORM SUPPORT:
#   - macOS with Homebrew HDF5 and Open MPI
#   - Linux with system or module HDF5 and MPI
#
# AUTHOR: Edward (Eddie) G. Nerney
# INSTITUTION: Laboratory for Atmospheric and Space Physics, CU Boulder
# DATE: November 2025
#==============================================================================

# MPI Fortran Compiler
FC = mpifort

# Executable name
TARGET = make_emission_movie_frames_mpi

# Source files
MODULE_SRC = IPT_emiss_MOP_community_code.f90
MAIN_SRC = make_emission_movie_frames_mpi.f90

# Object files
MODULE_OBJ = IPT_emiss_MOP_community_code.o
MAIN_OBJ = make_emission_movie_frames_mpi.o

# Module files
MODULES = ipt_emission_constants.mod \
          ipt_emission_utils.mod \
          ipt_emission_types.mod \
          ipt_emission_io.mod \
          ipt_emission_raytracer.mod

#==============================================================================
# HDF5 Configuration - Auto-detection
#==============================================================================

# Check for HDF5_DIR environment variable
ifdef HDF5_DIR
    HDF5_INC = -I$(HDF5_DIR)/include
    HDF5_LIB = -L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5
    HDF5_TYPE = environment
else
    # Check for Homebrew HDF5 on macOS
    BREW_HDF5 = $(shell brew --prefix hdf5 2>/dev/null)
    ifneq ($(BREW_HDF5),)
        HDF5_INC = -I$(BREW_HDF5)/include
        HDF5_LIB = -L$(BREW_HDF5)/lib -lhdf5_fortran -lhdf5
        HDF5_TYPE = homebrew
    else
        # Check common Linux paths
        ifneq ($(wildcard /usr/include/hdf5/serial),)
            # Debian/Ubuntu with serial HDF5
            HDF5_INC = -I/usr/include/hdf5/serial
            HDF5_LIB = -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_fortran -lhdf5
            HDF5_TYPE = system-debian
        else ifneq ($(wildcard /usr/include/hdf5.mod),)
            # Generic system HDF5
            HDF5_INC = -I/usr/include
            HDF5_LIB = -lhdf5_fortran -lhdf5
            HDF5_TYPE = system
        else ifneq ($(wildcard /usr/local/include/hdf5.mod),)
            # /usr/local installation
            HDF5_INC = -I/usr/local/include
            HDF5_LIB = -L/usr/local/lib -lhdf5_fortran -lhdf5
            HDF5_TYPE = local
        else
            # Fallback - try pkg-config
            HDF5_INC = $(shell pkg-config --cflags hdf5 2>/dev/null)
            HDF5_LIB = $(shell pkg-config --libs hdf5 2>/dev/null) -lhdf5_fortran
            HDF5_TYPE = pkg-config
        endif
    endif
endif

#==============================================================================
# Compiler Flags
#==============================================================================

# Common flags
FFLAGS_COMMON = -cpp -ffree-line-length-none -fimplicit-none $(HDF5_INC)

# Release flags (optimized for performance)
FFLAGS_RELEASE = $(FFLAGS_COMMON) -O3 -march=native -funroll-loops -ftree-vectorize

# Debug flags
FFLAGS_DEBUG = $(FFLAGS_COMMON) -g -O0 -Wall -Wextra -fcheck=all \
               -fbacktrace -ffpe-trap=invalid,zero,overflow

# Default to release mode
DEBUG ?= 0
ifeq ($(DEBUG),1)
    FFLAGS = $(FFLAGS_DEBUG)
else
    FFLAGS = $(FFLAGS_RELEASE)
endif

# Linker flags
LDFLAGS = $(HDF5_LIB)

#==============================================================================
# Build Rules
#==============================================================================

.PHONY: all clean debug release help info

# Default target
all: info $(TARGET)
	@echo ""
	@echo "Build successful!"
	@echo "Run with: mpirun -np <num_procs> ./$(TARGET)"
	@echo ""

# Release build
release:
	@$(MAKE) DEBUG=0 all

# Debug build
debug:
	@$(MAKE) DEBUG=1 all

# Link executable
$(TARGET): $(MODULE_OBJ) $(MAIN_OBJ)
	@echo ""
	@echo "Linking $(TARGET)..."
	$(FC) -o $@ $^ $(LDFLAGS)

# Compile module file (must be first)
$(MODULE_OBJ): $(MODULE_SRC)
	@echo "Compiling module: $<"
	$(FC) $(FFLAGS) -c $< -o $@

# Compile main program (depends on module)
$(MAIN_OBJ): $(MAIN_SRC) $(MODULE_OBJ)
	@echo "Compiling main program: $<"
	$(FC) $(FFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	rm -f $(TARGET) $(MODULE_OBJ) $(MAIN_OBJ) $(MODULES) *.mod
	rm -f output_data/*.bin
	@echo "Done."

# Print build configuration
info:
	@echo "=============================================================="
	@echo "IPT Emission Movie Generator (MPI) - Build Configuration"
	@echo "=============================================================="
	@echo "Compiler:   $(FC)"
	@echo "HDF5:       $(HDF5_TYPE)"
	@echo "Debug:      $(DEBUG)"
	@echo "=============================================================="
	@echo ""

# Help target
help:
	@echo ""
	@echo "IPT Emission Movie Frame Generator (MPI) - Makefile Help"
	@echo "========================================================="
	@echo ""
	@echo "Available targets:"
	@echo "  make              - Build optimized executable"
	@echo "  make debug        - Build with debugging enabled"
	@echo "  make release      - Build optimized (same as 'make')"
	@echo "  make clean        - Remove build artifacts"
	@echo "  make info         - Show build configuration"
	@echo "  make help         - Show this help message"
	@echo ""
	@echo "Environment variables:"
	@echo "  HDF5_DIR          - Set custom HDF5 installation path"
	@echo "  DEBUG=1           - Enable debug mode"
	@echo ""
	@echo "Example usage:"
	@echo "  make                                    # Build optimized"
	@echo "  mpirun -np 8 ./$(TARGET)  # Run with 8 processes"
	@echo ""
	@echo "Output:"
	@echo "  Binary data files in output_data/ directory"
	@echo "  Run generate_frames_from_binary.py to create PNG frames"
	@echo "  Run make_movie_ffmpeg.py to create MP4 movie"
	@echo ""
