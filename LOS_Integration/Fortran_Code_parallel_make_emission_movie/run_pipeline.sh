#!/bin/bash
#==============================================================================
# run_pipeline.sh
#
# Complete Pipeline for IPT Emission Movie Generation
# ====================================================
#
# This script runs the complete pipeline:
#   1. Build the Fortran/MPI executable
#   2. Run the MPI computation to generate binary data
#   3. Generate PNG frames from binary data
#   4. Create MP4 movie from frames
#
# USAGE:
#   ./run_pipeline.sh [num_procs]
#
# ARGUMENTS:
#   num_procs  Number of MPI processes (default: 4)
#
# EXAMPLE:
#   ./run_pipeline.sh 8    # Use 8 MPI processes
#
# AUTHOR: Edward (Eddie) G. Nerney
# INSTITUTION: Laboratory for Atmospheric and Space Physics, CU Boulder
# DATE: November 2025
#==============================================================================

set -e  # Exit on error

# Default number of MPI processes
NPROCS=${1:-4}

echo "======================================================================"
echo "Io Plasma Torus Emission Movie Generation Pipeline"
echo "======================================================================"
echo ""
echo "Configuration:"
echo "  MPI processes: ${NPROCS}"
echo ""

# Step 1: Build
echo "======================================================================"
echo "Step 1: Building Fortran/MPI executable..."
echo "======================================================================"
make clean
make
echo ""

# Step 2: Run MPI computation
echo "======================================================================"
echo "Step 2: Running MPI computation..."
echo "======================================================================"
mpirun -np ${NPROCS} ./make_emission_movie_frames_mpi
echo ""

# Step 3: Generate PNG frames
echo "======================================================================"
echo "Step 3: Generating PNG frames..."
echo "======================================================================"
python3 generate_frames_from_binary.py
echo ""

# Step 4: Create movie
echo "======================================================================"
echo "Step 4: Creating MP4 movie..."
echo "======================================================================"
python3 make_movie_ffmpeg.py
echo ""

echo "======================================================================"
echo "Pipeline Complete!"
echo "======================================================================"
echo ""
echo "Output movie: hybrid_emission.mp4"
echo ""
