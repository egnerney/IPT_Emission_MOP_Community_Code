# IPT Line-of-Sight Integration (MATLAB)

MATLAB implementation for UV and optical emission raytracing.

## Requirements

- MATLAB R2019b+ (no additional toolboxes required)

## Files

| File | Description |
|------|-------------|
| `IPT_emiss_MOP_community_code.m` | Main raytracer class |
| `basic_example_uv_integration_emission_model_use_tables.m` | UV example |
| `basic_example_optical_integration_emission_model_use_tables.m` | Optical example |

## Quick Start

```matlab
% Navigate to directory
cd('LOS_Integration/MATLAB_Code')

% Initialize raytracer (loads data files)
raytracer = IPT_emiss_MOP_community_code();

% Define observation geometry
slit_pos = [6.0, -20.0, 0.0];  % [R_J]
direction = [0.0, 1.0, 0.0];   % unit vector

% Calculate UV spectrum (single Maxwellian)
[wavelength, spectrum, lines] = raytracer.calculate_spectrum_single(...
    slit_pos, direction, [550, 2100], 1.0, 6.0, 0.01);

% Calculate optical spectrum (double Maxwellian)
[wavelength, spectrum, lines] = raytracer.calculate_spectrum_double(...
    slit_pos, direction, [3000, 10000], 2.0, 3.0, 0.01);

% Plot
fig = raytracer.plot_spectrum(wavelength, spectrum, lines);
```

## Running Examples

```matlab
results = basic_example_uv_integration_emission_model_use_tables();
results = basic_example_optical_integration_emission_model_use_tables();
```

## Key Methods

- `IPT_emiss_MOP_community_code()` — Constructor (loads data)
- `calculate_spectrum_single()` — Single Maxwellian LOS integration
- `calculate_spectrum_double()` — Double Maxwellian LOS integration
- `trace_ray()` — Generate ray sample points
- `plot_spectrum()` — Create publication-quality plots

## Output

- `wavelength` — Bin centers [Å]
- `spectrum` — Brightness [R/Å]
- `lines` — Cell array: `{wavelength, brightness, species}`
