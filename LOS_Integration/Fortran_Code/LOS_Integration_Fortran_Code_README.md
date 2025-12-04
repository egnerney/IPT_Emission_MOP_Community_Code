# IPT Line-of-Sight Integration (Fortran)

Production-quality Fortran implementation for UV and optical emission raytracing.

## Requirements

- gfortran 9.0+ (or Intel ifort)
- HDF5 with Fortran bindings
- gnuplot (for plotting)

**macOS:**
```bash
brew install gcc hdf5 gnuplot
```

**Linux:**
```bash
sudo apt install gfortran libhdf5-dev gnuplot
```

## Files

| File | Description |
|------|-------------|
| `IPT_emiss_MOP_community_code.f90` | Main raytracer module |
| `basic_example_uv_integration_emission_model_use_tables.f90` | UV example |
| `basic_example_optical_integration_emission_model_use_tables.f90` | Optical example |
| `Makefile` | Build system |

## Building

```bash
cd LOS_Integration/Fortran_Code
make clean
make
```

For debug build:
```bash
make debug
```

If HDF5 not auto-detected:
```bash
HDF5_DIR=/path/to/hdf5 make
```

## Running

```bash
./uv_emission_example
./optical_emission_example
```

## API

```fortran
use ipt_emission_raytracer

type(jovian_uv_raytracer) :: raytracer

! Initialize
call raytracer%initialize(ierr=ierr)

! Calculate spectrum
call raytracer%calculate_spectrum_single( &
    slit_pos, direction, &
    wave_range=[550.0_real64, 2100.0_real64], &
    bin_width=1.0_real64, &
    fwhm=6.0_real64, &
    ds=0.01_real64, &
    wave_bins=wave_bins, &
    spectrum=spectrum, &
    lines=lines)

! Cleanup
call raytracer%cleanup()
```

## Performance

Typical times on 6-core Intel i7:

| Calculation | Time |
|-------------|------|
| UV spectrum (single Maxwellian) | ~2 s |
| UV spectrum (double Maxwellian) | ~3 s |
| Optical spectrum | ~3–5 s |
