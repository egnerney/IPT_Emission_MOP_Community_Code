#==============================================================================
# Makefile for IPT UV/Optical Emission Raytracer
#
# Io Plasma Torus UV/Optical Emission Line-of-Sight Integration Code
# MOP (Magnetospheres of Outer Planets) Community Code
#
# USAGE:
#   make              - Build both UV and optical examples (release mode)
#   make uv           - Build UV emission example only
#   make optical      - Build optical emission example only
#   make debug        - Build with debugging enabled
#   make clean        - Remove build artifacts
#   make help         - Show available targets
#
# PREREQUISITES:
#   - gfortran (or compatible Fortran compiler)
#   - HDF5 library with Fortran bindings
#   - gnuplot (for generating PNG plots)
#
# PLATFORM SUPPORT:
#   - macOS with Homebrew HDF5
#   - Linux with system or module HDF5
#
# AUTHOR: Edward (Eddie) G. Nerney
# INSTITUTION: Laboratory for Atmospheric and Space Physics, CU Boulder
# DATE: November 2025
#==============================================================================

# Compiler
FC = gfortran

# Executable names
TARGET_UV = uv_emission_example
TARGET_OPTICAL = optical_emission_example

# Source files
MODULE_SRC = IPT_emiss_MOP_community_code.f90
UV_SRC = basic_example_uv_integration_emission_model_use_tables.f90
OPTICAL_SRC = basic_example_optical_integration_emission_model_use_tables.f90

# Object files
MODULE_OBJ = IPT_emiss_MOP_community_code.o
UV_OBJ = basic_example_uv_integration_emission_model_use_tables.o
OPTICAL_OBJ = basic_example_optical_integration_emission_model_use_tables.o

# Module files
MODULES = ipt_emission_constants.mod \
          ipt_emission_utils.mod \
          ipt_emission_types.mod \
          ipt_emission_io.mod \
          ipt_emission_raytracer.mod

#==============================================================================
# HDF5 Configuration - Auto-detection
#==============================================================================

# Try to detect HDF5 installation
# Priority: 1. Environment variables, 2. Homebrew, 3. System paths

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

# Release flags (optimized)
FFLAGS_RELEASE = $(FFLAGS_COMMON) -O3 -march=native -funroll-loops

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
LDFLAGS = -fopenmp $(HDF5_LIB)

#==============================================================================
# Build Rules
#==============================================================================

.PHONY: all clean debug release help info uv optical both

# Default target: build both executables
all: info both

# Build both UV and optical examples
both: $(TARGET_UV) $(TARGET_OPTICAL)

# Build UV example only
uv: info $(TARGET_UV)
	@echo ""
	@echo "UV emission example built successfully!"
	@echo "Run with: ./$(TARGET_UV)"
	@echo ""

# Build optical example only
optical: info $(TARGET_OPTICAL)
	@echo ""
	@echo "Optical emission example built successfully!"
	@echo "Run with: ./$(TARGET_OPTICAL)"
	@echo ""

# Release build
release:
	@$(MAKE) DEBUG=0 all

# Debug build
debug:
	@$(MAKE) DEBUG=1 all

# Link UV executable
$(TARGET_UV): $(MODULE_OBJ) $(UV_OBJ)
	@echo ""
	@echo "Linking $(TARGET_UV)..."
	$(FC) -o $@ $^ $(LDFLAGS)

# Link optical executable
$(TARGET_OPTICAL): $(MODULE_OBJ) $(OPTICAL_OBJ)
	@echo ""
	@echo "Linking $(TARGET_OPTICAL)..."
	$(FC) -o $@ $^ $(LDFLAGS)

# Compile module file (must be first)
$(MODULE_OBJ): $(MODULE_SRC)
	@echo "Compiling module: $<"
	$(FC) $(FFLAGS) -c $< -o $@

# Compile UV example program (depends on module)
$(UV_OBJ): $(UV_SRC) $(MODULE_OBJ)
	@echo "Compiling UV example: $<"
	$(FC) $(FFLAGS) -c $< -o $@

# Compile optical example program (depends on module)
$(OPTICAL_OBJ): $(OPTICAL_SRC) $(MODULE_OBJ)
	@echo "Compiling optical example: $<"
	$(FC) $(FFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	rm -f $(TARGET_UV) $(TARGET_OPTICAL) $(MODULE_OBJ) $(UV_OBJ) $(OPTICAL_OBJ) $(MODULES) *.mod
	rm -f *.dat *.gp *.png
	@echo "Done."

# Print build configuration
info:
	@echo "=============================================================="
	@echo "IPT UV/Optical Emission Code - Build Configuration"
	@echo "=============================================================="
	@echo "Compiler:   $(FC)"
	@echo "HDF5:       $(HDF5_TYPE)"
	@echo "Debug:      $(DEBUG)"
	@echo "=============================================================="
	@echo ""

# Help target
help:
	@echo ""
	@echo "IPT UV/Optical Emission Raytracer - Makefile Help"
	@echo "=================================================="
	@echo ""
	@echo "Available targets:"
	@echo "  make              - Build both UV and optical executables"
	@echo "  make uv           - Build UV emission example only"
	@echo "  make optical      - Build optical emission example only"
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
	@echo "Example:"
	@echo "  HDF5_DIR=/opt/hdf5 make"
	@echo "  make optical      # Build only optical example"
	@echo "  make uv           # Build only UV example"
	@echo ""
	@echo "Prerequisites:"
	@echo "  1. Run 'python convert_hdf5_strings.py' first"
	@echo "  2. Ensure HDF5 files are in correct directories"
	@echo ""
	@echo "Wavelength Ranges:"
	@echo "  UV:      550-2100 Angstroms (space-based spectrographs)"
	@echo "  Optical: 3000-10000 Angstroms (ground-based telescopes)"
	@echo ""
