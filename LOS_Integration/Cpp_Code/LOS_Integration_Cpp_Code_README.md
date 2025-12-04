# IPT Line-of-Sight Integration (C++)

C++ implementation for UV and optical emission raytracing.

## Requirements

- g++ 8.0+ (C++17) or clang++ 7.0+
- HDF5 with C++ bindings

**macOS:**
```bash
brew install gcc hdf5
```

**Linux:**
```bash
sudo apt install g++ libhdf5-dev libhdf5-cpp-dev
```

## Files

| File | Description |
|------|-------------|
| `IPT_emiss_MOP_community_code.hpp` | Header file |
| `IPT_emiss_MOP_community_code.cpp` | Implementation |
| `basic_example_uv_integration_emission_model_use_tables.cpp` | UV example |
| `basic_example_optical_integration_emission_model_use_tables.cpp` | Optical example |
| `Makefile` | Build system |

## Building

```bash
cd LOS_Integration/Cpp_Code
make           # Build both UV and optical
make uv        # UV only
make optical   # Optical only
make DEBUG=1   # Debug build
```

## Running

```bash
./ipt_emission_model           # UV
./ipt_optical_emission_model   # Optical
make run                       # Both
```

## API

```cpp
#include "IPT_emiss_MOP_community_code.hpp"

// Initialize raytracer (auto-loads data files)
ipt::JovianUVEmissionRaytracer raytracer;

// Define line of sight
std::array<double, 3> start = {6.0, -20.0, 0.0};
std::array<double, 3> dir = {0.0, 1.0, 0.0};

// Calculate UV spectrum
auto result = raytracer.calculate_spectrum_single(
    start, dir,
    {550.0, 2100.0},  // wavelength range [Å]
    1.0,              // bin width [Å]
    6.0,              // FWHM [Å]
    0.1               // step size [R_J]
);

// Access results
std::cout << "Total: " << result.total_brightness << " R\n";
```

## Performance

| Configuration | UV | Optical |
|---------------|-----|---------|
| Single Maxwellian | ~80–100 ms | ~150–200 ms |
| Double Maxwellian | ~80–100 ms | ~150–200 ms |
