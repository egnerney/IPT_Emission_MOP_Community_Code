# IPT cm³ Emission Model (IDL)

Volume emission rate calculations for the Io Plasma Torus. Supports both pre-calculated tables and direct CHIANTI database calls for custom table generation.

## Requirements

- IDL 8.0+ with HDF5 support
- CHIANTI 11.0.2 (required for direct CHIANTI calls and table generation)

## CHIANTI Database Setup

Download `chianti.zip` from [Zenodo (DOI: 10.5281/zenodo.17808352)](https://doi.org/10.5281/zenodo.17808352) and extract in the `cm3/` directory:

```bash
cd cm3
unzip chianti.zip
```

This creates `cm3/chianti/` with `dbase/` and `idl/` subdirectories.

Set the environment variable before starting IDL:
```bash
export CHIANTI_DATA=/path/to/IPT_Emission_MOP_Community_Code/cm3/chianti/dbase
```

## Files

| File | Description |
|------|-------------|
| `ipt_cm3_emission_model.pro` | Main library |
| `basic_example_uv_cm3_emission_model_use_tables.pro` | UV example using tables |
| `basic_example_uv_cm3_emission_model_use.pro` | UV example using direct CHIANTI |
| `basic_example_optical_cm3_emission_model_use_tables.pro` | Optical example |
| `make_emission_tables_chianti11_single_max.pro` | Generate single Maxwellian tables |
| `make_emission_tables_chianti11_double_max.pro` | Generate double Maxwellian tables |
| `run_make_emission_tables_chianti11_single_max.pro` | Runner for single Maxwellian table generation |
| `run_make_emission_tables_chianti11_double_max.pro` | Runner for double Maxwellian table generation |

## Quick Start

```idl
; Navigate to directory
cd, 'cm3/IDL_Code'

; Run UV example using pre-calculated tables
basic_example_uv_cm3_emission_model_use_tables

; Run optical example
basic_example_optical_cm3_emission_model_use_tables

; Run UV example using direct CHIANTI calls (requires CHIANTI setup)
basic_example_uv_cm3_emission_model_use
```

## Generating Custom Emission Tables

Requires CHIANTI database installation (see above).

```idl
cd, 'cm3/IDL_Code'

; Generate single Maxwellian tables
run_make_emission_tables_chianti11_single_max

; Generate double Maxwellian tables
run_make_emission_tables_chianti11_double_max
```

## Programmatic Usage

```idl
; Load emission tables
tables = IPT_LOAD_EMISSION_TABLES('../../Emiss_tables/CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5')

; Plasma parameters
Te = 5.0d  ; [eV]
ne = 2200.0d  ; [cm⁻³]
col_sp = 1.2d13  ; S+ column density [cm⁻²]
col_s2p = 4.2d13
col_s3p = 5.9d12
col_op = 5.2d13
col_o2p = 5.9d12

; Calculate emission
brightnesses = IPT_CALCULATE_EMISSION_SINGLE(tables, Te, ne, col_sp, col_s2p, col_s3p, col_op, col_o2p, WAVELENGTHS=wavelengths)

; Convolve spectrum
wave_grid = dindgen(1551) + 550.0d
spectrum = IPT_SIMULATE_SPECTRUM_ERF(wave_grid, 1.0d, wavelengths, brightnesses, FWHM=6.0d)
```

## Key Procedures

- `IPT_LOAD_EMISSION_TABLES` — Load HDF5 emission tables
- `IPT_CALCULATE_EMISSION_SINGLE` — Single Maxwellian emission
- `IPT_CALCULATE_EMISSION_DOUBLE` — Double Maxwellian emission
- `IPT_SIMULATE_SPECTRUM_ERF` — ERF-based Gaussian convolution
