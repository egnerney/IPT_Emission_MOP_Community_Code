# IPT Line-of-Sight Integration (Python)

Ray tracing through 3D plasma models for synthetic UV and optical emission spectra. Includes parallelized movie generation.

## Requirements

- Python 3.10+
- NumPy, SciPy, Matplotlib, h5py
- FFmpeg (for movie generation)

```bash
pip install numpy scipy matplotlib h5py
```

## Files

| File | Description |
|------|-------------|
| `IPT_emiss_MOP_community_code.py` | Main raytracer library |
| `basic_example_uv_integration_emission_model_use_tables.py` | UV example |
| `basic_example_optical_integration_emission_model_use_tables.py` | Optical example |
| `make_emission_frame_UV_example_low_res_parallelized.py` | Parallelized UV map generation |

## Quick Start

```python
from IPT_emiss_MOP_community_code import JovianUVEmissionRaytracer

# Initialize (loads plasma model and emission tables)
raytracer = JovianUVEmissionRaytracer()

# Define observation geometry
slit_pos = [6.0, -20.0, 0.0]   # [R_J]
direction = [0.0, 1.0, 0.0]    # unit vector

# Calculate UV spectrum (single Maxwellian)
wave, spectrum, lines = raytracer.calculate_spectrum_single(
    slit_pos, direction,
    wavelength_range=(550, 2100),
    bin_width=1.0,
    fwhm=6.0,
    ds=0.01
)

# Calculate optical spectrum (double Maxwellian)
wave, spectrum, lines = raytracer.calculate_spectrum_double(
    slit_pos, direction,
    wavelength_range=(3000, 10000),
    bin_width=2.0,
    fwhm=3.0,
    ds=0.01
)
```

## Running Examples

```bash
cd LOS_Integration/Python_Code
python basic_example_uv_integration_emission_model_use_tables.py
python basic_example_optical_integration_emission_model_use_tables.py
```

## Key Methods

- `JovianUVEmissionRaytracer()` — Initialize raytracer (auto-loads data)
- `calculate_spectrum_single()` — Single Maxwellian LOS integration
- `calculate_spectrum_double()` — Double Maxwellian LOS integration
- `trace_ray()` — Generate ray sample points

## Custom Models

```python
raytracer = JovianUVEmissionRaytracer(
    plasma_file='/path/to/custom_plasma.h5',
    emission_file_single='/path/to/custom_single.h5',
    emission_file_double='/path/to/custom_double.h5'
)
```

## Units

| Quantity | Unit |
|----------|------|
| Position | R_J (Jupiter radii) |
| Temperature | eV |
| Density | cm⁻³ |
| Wavelength | Å |
| Brightness | R (Rayleighs) |
