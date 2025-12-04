# IPT cm³ Emission Model (Python)

Volume emission rate calculations for the Io Plasma Torus using pre-calculated CHIANTI 11.0.2 emission tables.

## Requirements

- Python 3.9+
- NumPy, SciPy, Matplotlib, h5py

```bash
pip install numpy scipy matplotlib h5py
```

## Files

| File | Description |
|------|-------------|
| `IPT_cm3_emission_model.py` | Main library module |
| `basic_example_uv_cm3_emission_model_use_tables.py` | UV example (550–2100 Å) |
| `basic_example_optical_cm3_emission_model_use_tables.py` | Optical example (3000–10000 Å) |
| `convert_idl_emiss_tables_sav_files_to_h5.py` | IDL .sav to HDF5 converter |

## Quick Start

```python
from IPT_cm3_emission_model import EmissionTables, calculate_ipt_emiss_tables_single, simulate_ipt_spectrum_rayleighs_erf_form
import numpy as np

# Load emission tables
tables = EmissionTables()
tables.load_single_maxwellian_tables('../../Emiss_tables/CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5')

# Plasma parameters
Te, ne = 5.0, 2200.0  # [eV], [cm⁻³]
column_densities = {'S+': 1.2e13, 'S++': 4.2e13, 'S+++': 5.9e12, 'O+': 5.2e13, 'O++': 5.9e12}

# Calculate emission lines
wavelengths, brightnesses = calculate_ipt_emiss_tables_single(tables, Te, ne, column_densities)

# Convolve with instrument response
wave_grid = np.linspace(550, 2100, 1551)
spectrum = simulate_ipt_spectrum_rayleighs_erf_form(wave_grid, 1.0, wavelengths, brightnesses, fwhm=6.0)
```

## Running Examples

```bash
cd cm3/Python_Code
python basic_example_uv_cm3_emission_model_use_tables.py
python basic_example_optical_cm3_emission_model_use_tables.py
```

## Key Functions

- `EmissionTables.load_single_maxwellian_tables()` — Load single Maxwellian tables
- `EmissionTables.load_double_maxwellian_tables()` — Load double Maxwellian tables
- `calculate_ipt_emiss_tables_single()` — Calculate emission (single Maxwellian)
- `calculate_ipt_emiss_tables_double()` — Calculate emission (double Maxwellian)
- `simulate_ipt_spectrum_rayleighs_erf_form()` — ERF-based Gaussian convolution

## Units

| Quantity | Unit |
|----------|------|
| Temperature | eV |
| Density | cm⁻³ |
| Column density | cm⁻² |
| Wavelength | Å |
| Brightness | R (Rayleighs) |
