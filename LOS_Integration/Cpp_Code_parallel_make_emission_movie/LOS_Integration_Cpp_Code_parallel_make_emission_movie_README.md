# IPT Emission Movie Generator (C++/MPI)

MPI-parallelized C++ implementation for generating emission movie frames showing the torus at different Central Meridian Longitudes (CML).

## Requirements

- C++17 compiler with MPI (mpic++)
- HDF5 with C++ bindings
- Python 3.6+, NumPy, Matplotlib
- FFmpeg

**macOS:**
```bash
brew install open-mpi hdf5 ffmpeg
pip3 install numpy matplotlib
```

**Linux:**
```bash
sudo apt install build-essential libopenmpi-dev openmpi-bin libhdf5-dev libhdf5-cpp-103 ffmpeg
pip3 install numpy matplotlib
```

## Files

| File | Description |
|------|-------------|
| `make_emission_movie_frames_mpi.cpp` | MPI main program |
| `generate_frames_from_binary.py` | PNG frame generator |
| `make_movies_ffmpeg.py` | Movie creator |
| `Makefile` | Build system |

## Movie Output

Three-panel stacked view (S⁺ 6731 Å, S⁺⁺ 680 Å, O⁺ 833 Å) with 75 frames covering 360° rotation at 15 fps.

## Quick Start

```bash
cd LOS_Integration/Cpp_Code_parallel_make_emission_movie

# Build
make clean && make

# Run with 6 MPI processes
mpirun -np 6 ./make_emission_movie_frames_mpi

# Generate frames and movie
python3 generate_frames_from_binary.py
python3 make_movies_ffmpeg.py
```

## Performance

| MPI Ranks | Time (75 frames) |
|-----------|------------------|
| 1 | ~112 min |
| 4 | ~30 min |
| 6 | ~22 min |
| 8 | ~19 min |

Memory: ~1.5 GB per MPI rank

## Customization

Edit `make_emission_movie_frames_mpi.cpp`:

```cpp
namespace config {
    constexpr int FRAME_RATE = 15;
    constexpr int DURATION = 5;
    constexpr double R_OBS = 20.0;      // Observer distance [R_J]
    constexpr double DS = 0.1;          // Integration step [R_J]
    constexpr double S2P_680_VMAX = 185.0;  // S++ 680 Å colorbar max [R]
    constexpr double OP_833_VMAX = 175.0;   // O+ 833 Å colorbar max [R]
}
```

Rebuild after changes: `make clean && make`

Also update frame rate in `make_movies_ffmpeg.py` if changed.
