# =============================================================================
# Makefile for MPI-Parallelized Io Plasma Torus Emission Movie Generator
# =============================================================================
#
# This Makefile compiles the MPI C++ program for generating emission movie
# frames in parallel. The program uses the JovianUVEmissionRaytracer class
# from the adjacent Cpp_Code directory.
#
# USAGE:
#   make              Build the MPI executable
#   make run          Run with default 4 MPI processes
#   make run NP=8     Run with 8 MPI processes
#   make frames       Generate PNG frames from binary data
#   make movie        Create MP4 movie from frames
#   make all_steps    Run everything: compute, frames, movie
#   make clean        Remove compiled files
#   make clean_all    Remove compiled files and output directories
#   make help         Show this help
#
# AUTHOR: Edward (Eddie) G. Nerney
# INSTITUTION: LASP, University of Colorado Boulder
# DATE: November 2025
# =============================================================================

# Compiler settings - use MPI C++ compiler wrapper
CXX = mpic++
CXXFLAGS = -std=c++17 -O3 -march=native -mtune=native -funroll-loops -ftree-vectorize

# DO NOT use -ffast-math (breaks NaN handling required for proper physics)

# Executable name
TARGET = make_emission_movie_frames_mpi

# Source files
SRC = make_emission_movie_frames_mpi.cpp
RAYTRACER_SRC = ../Cpp_Code/IPT_emiss_MOP_community_code.cpp

# Object files
OBJ = make_emission_movie_frames_mpi.o IPT_emiss_MOP_community_code.o

# HDF5 configuration - auto-detect based on platform
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
    # macOS with Homebrew
    HDF5_DIR ?= $(shell brew --prefix hdf5 2>/dev/null || echo /usr/local)
    HDF5_INCLUDE = -I$(HDF5_DIR)/include
    HDF5_LIBS = -L$(HDF5_DIR)/lib -lhdf5 -lhdf5_cpp
else
    # Linux
    HDF5_INCLUDE = $(shell pkg-config --cflags hdf5 2>/dev/null || echo -I/usr/include/hdf5/serial)
    HDF5_LIBS = $(shell pkg-config --libs hdf5 2>/dev/null || echo -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -lhdf5_cpp)
    HDF5_LIBS += -lhdf5_cpp
endif

# Include paths
INCLUDES = $(HDF5_INCLUDE) -I../Cpp_Code

# Library paths
LIBS = $(HDF5_LIBS)

# Default number of MPI processes
NP ?= 4

# Output directories
MOVIE_DATA_DIR = movie_data
FRAMES_DIR = movie_frames/frames_stacked

# =============================================================================
# BUILD TARGETS
# =============================================================================

.PHONY: all clean clean_all run frames movie all_steps help

# Default target
all: $(TARGET)

# Link executable
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)
	@echo ""
	@echo "Build successful!"
	@echo "Run with: mpirun -np $(NP) ./$(TARGET)"

# Compile main source
make_emission_movie_frames_mpi.o: make_emission_movie_frames_mpi.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Compile raytracer library
IPT_emiss_MOP_community_code.o: $(RAYTRACER_SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# =============================================================================
# RUN TARGETS
# =============================================================================

# Run MPI program
run: $(TARGET)
	@echo "Running with $(NP) MPI processes..."
	mpirun -np $(NP) ./$(TARGET)

# Generate PNG frames from binary data
frames:
	@if [ -d "$(MOVIE_DATA_DIR)" ]; then \
		echo "Generating PNG frames..."; \
		python3 generate_frames_from_binary.py; \
	else \
		echo "ERROR: $(MOVIE_DATA_DIR) not found. Run 'make run' first."; \
		exit 1; \
	fi

# Create movie from frames
movie:
	@if [ -d "$(FRAMES_DIR)" ]; then \
		echo "Creating MP4 movie..."; \
		python3 make_movies_ffmpeg.py; \
	else \
		echo "ERROR: $(FRAMES_DIR) not found. Run 'make frames' first."; \
		exit 1; \
	fi

# Run all steps: compute, frames, movie
all_steps: run frames movie
	@echo ""
	@echo "All steps complete!"
	@echo "Movie: hybrid_emission.mp4"

# =============================================================================
# CLEAN TARGETS
# =============================================================================

# Clean compiled files only
clean:
	rm -f $(OBJ) $(TARGET)
	@echo "Cleaned compiled files."

# Clean everything including output directories
clean_all: clean
	rm -rf $(MOVIE_DATA_DIR) movie_frames hybrid_emission.mp4
	@echo "Cleaned all output files and directories."

# =============================================================================
# HELP TARGET
# =============================================================================

help:
	@echo "========================================================================"
	@echo "MPI Io Plasma Torus Emission Movie Generator - Makefile Help"
	@echo "========================================================================"
	@echo ""
	@echo "BUILD TARGETS:"
	@echo "  make              Build the MPI executable"
	@echo ""
	@echo "RUN TARGETS:"
	@echo "  make run          Run with default $(NP) MPI processes"
	@echo "  make run NP=8     Run with 8 MPI processes"
	@echo "  make frames       Generate PNG frames from binary data"
	@echo "  make movie        Create MP4 movie from frames"
	@echo "  make all_steps    Run everything: compute, frames, movie"
	@echo ""
	@echo "CLEAN TARGETS:"
	@echo "  make clean        Remove compiled files"
	@echo "  make clean_all    Remove compiled files and output directories"
	@echo ""
	@echo "WORKFLOW:"
	@echo "  1. make              # Build executable"
	@echo "  2. make run NP=8     # Run with 8 cores (adjust as needed)"
	@echo "  3. make frames       # Generate PNG images"
	@echo "  4. make movie        # Create MP4 video"
	@echo ""
	@echo "  Or simply: make all_steps NP=8"
	@echo ""
	@echo "OUTPUT FILES:"
	@echo "  movie_data/              Binary emission data files"
	@echo "  movie_frames/            PNG frame images"
	@echo "  hybrid_emission.mp4      Final movie"
	@echo ""
	@echo "REQUIREMENTS:"
	@echo "  - MPI (mpic++ compiler)"
	@echo "  - HDF5 library with C++ bindings"
	@echo "  - Python 3 with numpy, matplotlib"
	@echo "  - FFmpeg for movie creation"
	@echo "========================================================================"
