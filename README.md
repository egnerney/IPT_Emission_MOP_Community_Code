# IPT_Emission_MOP_Community_Code

**Io Plasma Torus UV and Optical Emission Modeling for the Magnetospheres of Outer Planets (MOP) Community**

[![License: Open Source](https://img.shields.io/badge/License-Open%20Source-green.svg)](LICENSE)
[![CHIANTI: v11.0.2](https://img.shields.io/badge/CHIANTI-v11.0.2-blue.svg)](https://www.chiantidatabase.org/)

## Overview

This repository provides a comprehensive toolkit for modeling UV and optical emission spectra from the Io Plasma Torus (IPT). The code supports two complementary approaches:

1. **Volume emission rates (cm⁻³ approach)**: Calculate emission spectra given plasma parameters and column densities
2. **Line-of-sight (LOS) integration**: Ray tracing through 3D plasma models with trilinear interpolation

The emission models use pre-calculated tables from the [CHIANTI 11.0.2 atomic database](https://www.chiantidatabase.org/) and support both single and double Maxwellian electron velocity distributions, enabling accurate modeling of both thermal and suprathermal (hot) electron populations observed in the Jovian magnetosphere.

Designed for the Magnetospheres of Outer Planets (MOP) research community for interpreting observations from missions including Juno, Europa-UVS, JUICE-UVS, and HST.

## Features

- **Multi-language support**: Python, IDL/GDL, and Fortran implementations
- **Emission table generation**: Create custom emission tables from CHIANTI atomic database
- **Single & double Maxwellian distributions**: Model both thermal core and hot electron populations
- **Line-of-sight integration**: Ray tracing through 3D plasma models with trilinear interpolation
- **UV and optical wavelength coverage**: 550–2100 Å (UV) and 3000–10000 Å (optical)
- **Instrument convolution**: Gaussian PSF convolution with ERF-based integration
- **Parallelized computation**: Multi-core support for emission map generation
- **Publication-quality plotting**: Matplotlib (Python) and gnuplot (Fortran) output
- **Cross-platform**: Tested on macOS and Linux

## Repository Structure

```
IPT_Emission_MOP_Community_Code/
│
├── cm3/                                    # Volume emission rate calculations
│   │
│   ├── Python/
│   │   ├── IPT_cm3_emission_model.py                         # Main Python library
│   │   ├── basic_example_uv_cm3_emission_model_use_tables.py # UV example
│   │   ├── basic_example_optical_cm3_emission_model_use_tables.py  # Optical example
│   │   └── convert_idl_emiss_tables_sav_files_to_h5.py       # IDL→HDF5 converter
│   │
│   ├── IDL/
│   │   ├── ipt_cm3_emission_model.pro                        # Main IDL library
│   │   ├── basic_example_uv_cm3_emission_model_use.pro       # UV example (direct CHIANTI)
│   │   ├── basic_example_uv_cm3_emission_model_use_tables.pro    # UV example (tables)
│   │   ├── basic_example_optical_cm3_emission_model_use_tables.pro   # Optical example
│   │   ├── run_basic_example_uv_cm3_emission_model_use.pro   # Runner script
│   │   ├── make_emission_tables_chianti11_single_max.pro     # Single Maxwellian table generator
│   │   ├── make_emission_tables_chianti11_double_max.pro     # Double Maxwellian table generator
│   │   ├── run_make_emission_tables_chianti11_single_max.pro # Runner script
│   │   └── run_make_emission_tables_chianti11_double_max.pro # Runner script
│   │
│   └── chianti/                            # CHIANTI database (download required)
│       ├── dbase/                          # Atomic database files
│       └── idl/                            # CHIANTI IDL procedures
│
├── LOS_Integration/                        # Line-of-sight integration
│   │
│   ├── Python/
│   │   ├── IPT_emiss_MOP_community_code.py                   # Main raytracer library
│   │   ├── basic_example_uv_integration_emission_model_use_tables.py     # UV example
│   │   ├── basic_example_uv_integration_emission_model_use_tables_aligned.py  # Aligned geometry
│   │   ├── basic_example_optical_integration_emission_model_use_tables.py    # Optical example
│   │   ├── make_emission_frame_UV_example_low_res.py         # UV emission maps (serial)
│   │   └── make_emission_frame_UV_example_low_res_parallelized.py  # UV maps (parallel)
│   │
│   ├── IDL/
│   │   ├── IPT_emiss_MOP_community_code.pro                  # Main raytracer library
│   │   ├── basic_example_uv_integration_emission_model_use_tables.pro    # UV example
│   │   └── basic_example_optical_integration_emission_model_use_tables.pro   # Optical example
│   │
│   └── Fortran/
│       ├── IPT_emiss_MOP_community_code.f90                  # Main Fortran module
│       ├── basic_example_uv_integration_emission_model_use_tables.f90    # UV example
│       └── makefile                                          # Build system (macOS/Linux)
│
├── Emiss_tables/                           # Pre-calculated emission tables (download required)
│   ├── CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5
│   └── CHIANTI_11.0.2_emiss_tables_double_maxwellian_24x10x24x12.h5
│
├── 3D_Torus_Model/                         # 3D plasma model (download required)
│   └── jovian_plasma_interpolated_381x381x231.h5
│
└── README.md
```

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/egnerney/IPT_Emission_MOP_Community_Code.git
cd IPT_Emission_MOP_Community_Code
```

### 2. Download Required Data Files

Due to file size limitations on GitHub, several data files must be downloaded separately from Google Drive. All links are publicly accessible (no sign-in required).

#### Emission Tables (place in `Emiss_tables/`)

| File | Description | Download Link |
|------|-------------|---------------|
| `CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5` | Single Maxwellian (50×50 Te-ne grid) | [Download](https://drive.google.com/file/d/1w1Pw1h2otMhqbR4d9V6XogovK1xaWJe6/view?usp=sharing) |
| `CHIANTI_11.0.2_emiss_tables_double_maxwellian_24x10x24x12.h5` | Double Maxwellian (24×10×24×12 grid) | [Download](https://drive.google.com/file/d/1WvRH7Zyk_-cb03L4CsesYm20LWjvea3a/view?usp=sharing) |

#### 3D Plasma Model (place in `3D_Torus_Model/`)

| File | Description | Download Link |
|------|-------------|---------------|
| `jovian_plasma_interpolated_381x381x231.h5` | 3D torus model (381×381×231 grid) | [Download](https://drive.google.com/file/d/1eqFOvSBFJwJvGc4Yd6O8cBSr7QwfOoGC/view?usp=sharing) |

#### CHIANTI Database (place in `cm3/` and unzip)

| File | Description | Download Link |
|------|-------------|---------------|
| `chianti.zip` | CHIANTI 11.0.2 database + IDL procedures | [Download](https://drive.google.com/file/d/1hdDlgyza4p88vOlViORC0cbZq35JuZxS/view?usp=sharing) |

> **Note**: The CHIANTI database is originally from the [CHIANTI Database website](https://www.chiantidatabase.org/chianti_download.html). The provided archive includes both the atomic database files and required IDL procedures.

After downloading, organize the files:

```bash
# Move emission tables to Emiss_tables/
mv CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5 Emiss_tables/
mv CHIANTI_11.0.2_emiss_tables_double_maxwellian_24x10x24x12.h5 Emiss_tables/

# Move 3D plasma model to 3D_Torus_Model/
mv jovian_plasma_interpolated_381x381x231.h5 3D_Torus_Model/

# Extract CHIANTI database in cm3/
cd cm3
unzip chianti.zip
cd ..
```

### 3. Install Dependencies

#### Python

```bash
pip install numpy scipy matplotlib h5py
```

#### IDL/GDL

Requires IDL 8.0+ or GDL with HDF5 support. CHIANTI IDL procedures are included in the download.

Set the CHIANTI environment variable:
```bash
export CHIANTI_DATA=/path/to/IPT_Emission_MOP_Community_Code/cm3/chianti/dbase
```

#### Fortran

**macOS (Homebrew):**
```bash
brew install gcc hdf5 gnuplot
```

**Linux (Debian/Ubuntu):**
```bash
sudo apt install gfortran libhdf5-dev gnuplot
```

## Quick Start

### Python: Volume Emission Calculation (cm³ approach)

```python
from cm3.IPT_cm3_emission_model import (
    EmissionTables,
    calculate_ipt_emiss_tables_single,
    simulate_ipt_spectrum_rayleighs_erf_form
)
import numpy as np

# Load emission tables
tables = EmissionTables()
tables.load_single_maxwellian_tables('Emiss_tables/CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5')

# Define plasma parameters
Te = 5.0      # Electron temperature [eV]
ne = 2200.0   # Electron density [cm⁻³]

# Column densities [cm⁻²] for typical IPT conditions
column_densities = {
    'S+': 1.2e13, 'S++': 4.2e13, 'S+++': 5.92e12,
    'S++++': 6.0e11, 'O+': 5.2e13, 'O++': 5.92e12
}

# Calculate discrete emission lines
xwavi, yptsi = calculate_ipt_emiss_tables_single(
    tables, Te, ne, column_densities, min_wav=550.0, max_wav=2100.0
)

# Convolve with instrument response (FWHM = 6 Å for Europa-UVS/JUICE-UVS)
xwav = np.linspace(550, 2100, 1551)
spectrum = simulate_ipt_spectrum_rayleighs_erf_form(xwav, 1.0, xwavi, yptsi, fwhm=6.0)
```

### Python: Line-of-Sight Integration

```python
from LOS_Integration.IPT_emiss_MOP_community_code import JovianUVEmissionRaytracer
import numpy as np

# Initialize raytracer (loads plasma model and emission tables)
raytracer = JovianUVEmissionRaytracer()

# Define observation geometry
slit_pos = np.array([6.0, -20.0, 0.0])  # Start at x=6 R_J, equator
direction = np.array([0.0, 1.0, 0.0])   # Look in +y direction

# Calculate UV spectrum
wave_bins, spectrum, line_list = raytracer.calculate_spectrum_single(
    slit_pos, direction,
    wavelength_range=(550, 2100),
    bin_width=1.0,
    fwhm=6.0,
    ds=0.01  # Integration step size [R_J]
)
```

### Fortran: Line-of-Sight Integration

```bash
cd LOS_Integration/Fortran
make              # Build with OpenMP (parallel)
./uv_emission_example
```

The Fortran example traces rays through the 3D plasma model and calculates UV spectra for both equatorial and off-equatorial lines of sight.

### IDL: Volume Emission Calculation

```idl
; Run from cm3/IDL/ directory
@run_basic_example_uv_cm3_emission_model_use
```

### IDL: Line-of-Sight Integration

```idl
; Compile the library first
.compile IPT_emiss_MOP_community_code.pro

; Run the example
@basic_example_uv_integration_emission_model_use_tables
```

## Physical Background

### The Io Plasma Torus

The Io Plasma Torus is a donut-shaped region of plasma in Jupiter's magnetosphere, fed by volcanic activity on Io (~1 ton/s of SO₂). Key characteristics:

| Parameter | Value | Description |
|-----------|-------|-------------|
| Location | ~5–10 R_J | Cylindrical radius from Jupiter's spin axis |
| Peak density | ~2000 cm⁻³ | Near 5.9 R_J (Io's orbit) |
| Scale height | 0.5–1 R_J | Perpendicular to centrifugal equator |
| Core temperature | ~5 eV | Thermal electron population |
| Hot temperature | ~50–500 eV | Suprathermal electrons |
| Hot fraction | ~0.25% | Fraction of hot electrons |

### Electron Distributions

The code supports two electron velocity distributions:

**Single Maxwellian**: Thermal equilibrium distribution characterized by a single temperature T_e. Appropriate for most torus conditions.

**Double Maxwellian**: Superposition of cold (core) and hot populations:
- Core: T_ec ~ 5 eV, fraction f_ec ~ 99.75%
- Hot: T_eh ~ 270 eV, fraction f_eh ~ 0.25%

Hot electrons significantly enhance high-excitation emission lines (e.g., S IV, O III) relative to a single Maxwellian.

### Key Diagnostic Emission Lines

| Ion | Wavelength (Å) | Transition | Diagnostic Use |
|-----|----------------|------------|----------------|
| S II | 765.0 | 4S–4P | Density diagnostic |
| S III | 680.4 | 3P–3D | Temperature diagnostic |
| S III | 729.5 | 3P–3S | Abundance diagnostic |
| O II | 833.3 | 4S–4P | O/S ratio |
| S IV | 657.3 | 2P–2D | Hot electron diagnostic |

### Coordinate System and Units

| Quantity | Units | Description |
|----------|-------|-------------|
| Position | R_J | Jupiter radii (1 R_J = 71,492 km) |
| Temperature | eV | Electron volts |
| Density | cm⁻³ | Particles per cubic centimeter |
| Wavelength | Å | Angstroms |
| Emission rate | photons s⁻¹ ion⁻¹ | Per-ion photon emission rate |
| Brightness | R | Rayleighs (10⁶ photons s⁻¹ cm⁻² (4π sr)⁻¹) |

## Generating Custom Emission Tables

To generate emission tables with different parameter grids using IDL:

```idl
; Single Maxwellian tables
@run_make_emission_tables_chianti11_single_max

; Double Maxwellian tables  
@run_make_emission_tables_chianti11_double_max
```

Convert IDL `.sav` files to HDF5 for use with Python/Fortran:

```bash
python cm3/convert_idl_emiss_tables_sav_files_to_h5.py input.sav output.h5
```

## API Reference

### Python cm³ Module (`IPT_cm3_emission_model.py`)

```python
class EmissionTables:
    """Container for CHIANTI emission rate tables."""
    def load_single_maxwellian_tables(filepath: str) -> None
    def load_double_maxwellian_tables(filepath: str) -> None

def calculate_ipt_emiss_tables_single(
    tables: EmissionTables,
    Te: float,              # Electron temperature [eV]
    ne: float,              # Electron density [cm⁻³]
    column_densities: dict, # Ion column densities [cm⁻²]
    min_wav: float = 550.0,
    max_wav: float = 2100.0
) -> Tuple[np.ndarray, np.ndarray]

def calculate_ipt_emiss_tables_double(
    tables: EmissionTables,
    Tec: float,             # Core temperature [eV]
    Teh: float,             # Hot temperature [eV]
    ne_total: float,        # Total electron density [cm⁻³]
    feh: float,             # Hot electron fraction
    column_densities: dict,
    min_wav: float = 550.0,
    max_wav: float = 2100.0
) -> Tuple[np.ndarray, np.ndarray]

def simulate_ipt_spectrum_rayleighs_erf_form(
    wavelength_grid: np.ndarray,
    bin_width: float,
    line_wavelengths: np.ndarray,
    line_brightnesses: np.ndarray,
    fwhm: float = 6.0
) -> np.ndarray
```

### Python LOS Module (`IPT_emiss_MOP_community_code.py`)

```python
class JovianUVEmissionRaytracer:
    """Main class for LOS-integrated emission calculations."""
    
    def __init__(
        self,
        plasma_file: str = None,
        emission_file_single: str = None,
        emission_file_double: str = None
    )
    
    def calculate_spectrum_single(
        self,
        slit_pos_vec: np.ndarray,   # Starting position [R_J]
        norm_vec: np.ndarray,        # Direction vector (normalized)
        wavelength_range: tuple = (550, 2100),
        bin_width: float = 1.0,
        fwhm: float = 6.0,
        ds: float = 0.01             # Step size [R_J]
    ) -> Tuple[np.ndarray, np.ndarray, list]
    
    def calculate_spectrum_double(
        self,
        slit_pos_vec: np.ndarray,
        norm_vec: np.ndarray,
        wavelength_range: tuple = (550, 2100),
        bin_width: float = 1.0,
        fwhm: float = 6.0,
        ds: float = 0.01
    ) -> Tuple[np.ndarray, np.ndarray, list]
    
    def trace_ray(
        self,
        start_pos: np.ndarray,
        direction: np.ndarray,
        ds: float = 0.01,
        max_distance: float = 40.0
    ) -> Tuple[np.ndarray, np.ndarray]
```

### Fortran Module (`IPT_emiss_MOP_community_code.f90`)

```fortran
module ipt_emission_raytracer
    type :: jovian_uv_raytracer
        procedure :: initialize
        procedure :: calculate_spectrum_single
        procedure :: calculate_spectrum_double
        procedure :: trace_ray
        procedure :: interpolate_plasma_trilinear
        procedure :: cleanup
    end type
end module
```

## References

If you use this code in your research, please cite the following:

### Software Citation
- Nerney, E. G. (2025). IPT_Emission_MOP_Community_Code. GitHub repository. https://github.com/egnerney/IPT_Emission_MOP_Community_Code

### Scientific References
- Nerney, E. G., Bagenal, F., & Steffl, A. J. (2017). Io plasma torus ion composition: Voyager, Galileo, and Cassini. *Journal of Geophysical Research: Space Physics*, 122, 727-744.
- Nerney, E. G., & Bagenal, F. (2020). Combining UV spectra and physical chemistry to constrain Io plasma torus composition. *Journal of Geophysical Research: Space Physics*, 125, e2019JA027458.
- Steffl, A. J., Stewart, A. I. F., & Bagenal, F. (2004). Cassini UVIS observations of the Io plasma torus. I. Initial results. *Icarus*, 172, 78-90.
- Thomas, N., Bagenal, F., Hill, T. W., & Wilson, J. K. (2004). The Io neutral clouds and plasma torus. In *Jupiter: The Planet, Satellites and Magnetosphere* (pp. 561-591). Cambridge University Press.
- Del Zanna, G., Dere, K. P., Young, P. R., & Landi, E. (2021). CHIANTI—An atomic database for emission lines. XVI. Version 10. *The Astrophysical Journal*, 909, 38.

## Author

**Edward (Eddie) G. Nerney**  
Laboratory for Atmospheric and Space Physics  
University of Colorado Boulder  
edward.nerney@colorado.edu

## License

This code is open source for academic and research use. See [LICENSE](LICENSE) for details.

The CHIANTI atomic database is subject to its own license. Please see the [CHIANTI website](https://www.chiantidatabase.org/) for terms of use.

## Contributing

Contributions are welcome! Please open an issue or submit a pull request for:

- Bug fixes
- New emission line species
- Performance improvements
- Documentation updates
- Additional example scripts

## Acknowledgments

- The CHIANTI team for the atomic database
- NASA Juno, Europa Clipper, and ESA JUICE mission teams
- The Magnetospheres of Outer Planets (MOP) community for feedback and testing
- NASA Solar System Workings program for funding support

---

**CHIANTI Acknowledgment**: CHIANTI is a collaborative project involving George Mason University, the University of Michigan (USA), University of Cambridge (UK), and NASA Goddard Space Flight Center (USA).
