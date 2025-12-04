# =============================================================================
# Makefile for IPT UV/Optical Emission Line-of-Sight Integration (C++)
# =============================================================================
#
# Io Plasma Torus (IPT) Emission Model - MOP Community Code
#
# Author: Edward (Eddie) G. Nerney
# Institution: LASP, University of Colorado Boulder
# Date: November 2025
#
# IMPORTANT: Do NOT use -ffast-math as it breaks NaN handling and causes
# numerical instability in floating point comparisons.
# =============================================================================

CXX := g++
CXXSTD := -std=c++17

# Base warning flags
WARN_FLAGS := -Wall -Wextra -Wpedantic

# Optimization flags - NOTE: NO -ffast-math for numerical stability
OPT_FLAGS := -O3 -march=native -mtune=native -funroll-loops -ftree-vectorize

# Debug flags
DEBUG_FLAGS := -g -O0 -DDEBUG

ifdef DEBUG
    CXXFLAGS := $(CXXSTD) $(WARN_FLAGS) $(DEBUG_FLAGS)
else
    CXXFLAGS := $(CXXSTD) $(WARN_FLAGS) $(OPT_FLAGS)
endif

# HDF5 Configuration - Auto-detect platform
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
    # macOS with Homebrew
    HDF5_DIR := $(shell brew --prefix hdf5 2>/dev/null || echo /usr/local/opt/hdf5)
    HDF5_CFLAGS := -I$(HDF5_DIR)/include
    HDF5_LIBS := -L$(HDF5_DIR)/lib -lhdf5_cpp -lhdf5 -lz -ldl -lm
else
    # Linux - try pkg-config first, fall back to common paths
    HDF5_CFLAGS := $(shell pkg-config --cflags hdf5 2>/dev/null || echo -I/usr/include/hdf5/serial)
    HDF5_LIBS := $(shell pkg-config --libs hdf5 2>/dev/null || echo -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5)
    HDF5_LIBS += -lhdf5_cpp
endif

# Source files
LIB_SRC := IPT_emiss_MOP_community_code.cpp
LIB_OBJ := IPT_emiss_MOP_community_code.o
HDRS := IPT_emiss_MOP_community_code.hpp

UV_SRC := basic_example_uv_integration_emission_model_use_tables.cpp
UV_OBJ := basic_example_uv_integration_emission_model_use_tables.o
UV_TARGET := ipt_emission_model

OPTICAL_SRC := basic_example_optical_integration_emission_model_use_tables.cpp
OPTICAL_OBJ := basic_example_optical_integration_emission_model_use_tables.o
OPTICAL_TARGET := ipt_optical_emission_model

# Phony targets
.PHONY: all clean uv optical run run_uv run_optical plot help

# Default: build both UV and optical executables
all: $(UV_TARGET) $(OPTICAL_TARGET)

# Build UV example only
uv: $(UV_TARGET)

# Build optical example only
optical: $(OPTICAL_TARGET)

# Link UV executable
$(UV_TARGET): $(LIB_OBJ) $(UV_OBJ)
	@echo "Linking $(UV_TARGET)..."
	$(CXX) $(CXXFLAGS) -o $@ $^ $(HDF5_LIBS)
	@echo ""
	@echo "UV build complete: $(UV_TARGET)"
	@echo "To run: ./$(UV_TARGET)"
	@echo "To plot results: python3 plot_spectra.py"

# Link optical executable
$(OPTICAL_TARGET): $(LIB_OBJ) $(OPTICAL_OBJ)
	@echo "Linking $(OPTICAL_TARGET)..."
	$(CXX) $(CXXFLAGS) -o $@ $^ $(HDF5_LIBS)
	@echo ""
	@echo "Optical build complete: $(OPTICAL_TARGET)"
	@echo "To run: ./$(OPTICAL_TARGET)"
	@echo "To plot results: python3 plot_optical_spectra.py"

# Compile source files
%.o: %.cpp $(HDRS)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) $(HDF5_CFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	rm -f $(LIB_OBJ) $(UV_OBJ) $(OPTICAL_OBJ) $(UV_TARGET) $(OPTICAL_TARGET)
	rm -f *.dat *.png plot_spectra.py plot_optical_spectra.py
	@echo "Clean complete."

# Run UV example
run_uv: $(UV_TARGET)
	@echo ""
	./$(UV_TARGET)

# Run optical example
run_optical: $(OPTICAL_TARGET)
	@echo ""
	./$(OPTICAL_TARGET)

# Run both
run: run_uv run_optical

# Generate plots
plot:
	@if [ -f plot_spectra.py ]; then python3 plot_spectra.py; fi
	@if [ -f plot_optical_spectra.py ]; then python3 plot_optical_spectra.py; fi

# Help target
help:
	@echo "IPT UV/Optical Emission Model - Build Targets"
	@echo "=============================================="
	@echo ""
	@echo "Targets:"
	@echo "  all          - Build both UV and optical executables (default)"
	@echo "  uv           - Build UV example only"
	@echo "  optical      - Build optical example only"
	@echo "  clean        - Remove build artifacts"
	@echo "  run_uv       - Build and run UV example"
	@echo "  run_optical  - Build and run optical example"
	@echo "  run          - Build and run both examples"
	@echo "  plot         - Generate plots from output data"
	@echo "  help         - Show this help message"
	@echo ""
	@echo "Options:"
	@echo "  DEBUG=1      - Build with debug symbols"
	@echo ""
	@echo "Example usage:"
	@echo "  make               # Build both executables"
	@echo "  make optical       # Build optical example only"
	@echo "  make run_optical   # Build and run optical example"
	@echo "  make DEBUG=1       # Debug build"
