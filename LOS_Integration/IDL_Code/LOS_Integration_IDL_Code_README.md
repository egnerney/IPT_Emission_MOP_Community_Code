# IPT Line-of-Sight Integration (IDL)

Ray tracing through 3D plasma models for synthetic UV and optical emission spectra.

## Requirements

- IDL 8.5+ with HDF5 support

## Files

| File | Description |
|------|-------------|
| `IPT_emiss_MOP_community_code.pro` | Main raytracer library |
| `basic_example_uv_integration_emission_model_use_tables.pro` | UV example |
| `basic_example_optical_integration_emission_model_use_tables.pro` | Optical example |

## Quick Start

```idl
; Compile library
.COMPILE IPT_emiss_MOP_community_code

; Initialize raytracer (loads all data files)
raytracer = IPT_INIT_RAYTRACER()

; Define observation geometry
slit_pos = [6.0d, -20.0d, 0.0d]  ; [R_J]
direction = [0.0d, 1.0d, 0.0d]   ; unit vector

; Calculate UV spectrum (single Maxwellian)
spectrum = IPT_CALCULATE_SPECTRUM_SINGLE(raytracer, slit_pos, direction, $
    WAVELENGTH_RANGE=[550.0d, 2100.0d], BIN_WIDTH=1.0d, FWHM=6.0d, DS=0.01d, $
    WAVE_BINS=wave_bins, LINE_LIST=lines)

; Calculate optical spectrum (double Maxwellian)
spectrum = IPT_CALCULATE_SPECTRUM_DOUBLE(raytracer, slit_pos, direction, $
    WAVELENGTH_RANGE=[3000.0d, 10000.0d], BIN_WIDTH=2.0d, FWHM=3.0d, DS=0.01d, $
    WAVE_BINS=wave_bins, LINE_LIST=lines)
```

## Running Examples

```idl
cd, 'LOS_Integration/IDL_Code'
.RUN basic_example_uv_integration_emission_model_use_tables
.RUN basic_example_optical_integration_emission_model_use_tables
```

## Key Procedures

- `IPT_INIT_RAYTRACER` — Initialize and load data files
- `IPT_CALCULATE_SPECTRUM_SINGLE` — Single Maxwellian LOS integration
- `IPT_CALCULATE_SPECTRUM_DOUBLE` — Double Maxwellian LOS integration
- `IPT_TRACE_RAY` — Generate ray sample points

## Output

- `spectrum` — Convolved spectrum [R/Å]
- `wave_bins` — Wavelength bin centers [Å]
- `line_list` — Structure array with `.wavelength`, `.brightness`, `.species`
