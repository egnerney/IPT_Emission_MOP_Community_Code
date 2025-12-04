# IPT Emission Movie Generator (Fortran/MPI)

MPI-parallelized Fortran implementation for generating emission movie frames showing the torus at different Central Meridian Longitudes (CML).

## Requirements

- gfortran 9.0+ with MPI (mpifort)
- HDF5 with Fortran bindings
- Python 3.8+, NumPy, Matplotlib
- FFmpeg

**macOS:**
```bash
brew install gcc open-mpi hdf5 ffmpeg
pip3 install numpy matplotlib h5py
```

**Linux:**
```bash
sudo apt install gfortran openmpi-bin libopenmpi-dev libhdf5-dev ffmpeg
pip3 install numpy matplotlib h5py
```

## Files

| File | Description |
|------|-------------|
| `IPT_emiss_MOP_community_code.f90` | Core raytracer module |
| `make_emission_movie_frames_mpi.f90` | MPI main program |
| `generate_frames_from_binary.py` | PNG frame generator |
| `make_movie_ffmpeg.py` | Movie creator |
| `run_pipeline.sh` | Complete workflow script |

## Movie Output

Three-panel stacked view (S⁺ 6731 Å, S⁺⁺ 680 Å, O⁺ 833 Å) with 75 frames covering 360° rotation at 15 fps.

## Quick Start

```bash
cd LOS_Integration/Fortran_Code_parallel_make_emission_movie

# Complete pipeline with 8 MPI processes
./run_pipeline.sh 8
```

## Step-by-Step

```bash
# 1. Build
make clean && make

# 2. Generate binary data (Fortran/MPI)
mpirun -np 6 ./make_emission_movie_frames_mpi

# 3. Generate PNG frames (Python)
python3 generate_frames_from_binary.py

# 4. Create movie (FFmpeg)
python3 make_movie_ffmpeg.py
```

## Performance

| MPI Processes | Time (75 frames) |
|---------------|------------------|
| 1 | ~3–4 hours |
| 4 | ~45–60 min |
| 6 | ~30–40 min |
| 8 | ~25–35 min |

## Customization

Edit `make_emission_movie_frames_mpi.f90`:

```fortran
! Frame rate and duration
integer(int32), parameter :: FRAME_RATE = 15
integer(int32), parameter :: DURATION = 5

! UV colorbar limits
real(real64), parameter :: S2P_680_VMAX = 185.0_real64
real(real64), parameter :: OP_833_VMAX = 175.0_real64
```

Rebuild after changes: `make clean && make`
