# IPT Emission MOP Community Code

**Io Plasma Torus UV and Optical Emission Modeling for the Magnetospheres of Outer Planets (MOP) Community**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17808825.svg)](https://doi.org/10.5281/zenodo.17808825)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![CHIANTI: v11.0.2](https://img.shields.io/badge/CHIANTI-v11.0.2-orange.svg)](https://www.chiantidatabase.org/)

## Overview

A comprehensive toolkit for modeling UV (550–2100 Å) and optical (3000–10000 Å) emission spectra from the Io Plasma Torus (IPT). Supports volume emission rate calculations and line-of-sight integration through 3D plasma models using pre-calculated CHIANTI 11.0.2 atomic database tables.

Designed for interpreting observations from Juno-UVS, JUICE-UVS, Europa-UVS, HST/STIS, and ground-based telescopes.

## Features

- **Two modeling approaches**: Volume emission rates (cm³) and line-of-sight integration
- **Multi-language**: Python, IDL/GDL, Fortran, C++, MATLAB
- **Electron distributions**: Single and double Maxwellian (core + hot populations)
- **Parallelized workflows**: MPI-based emission map and movie generation
- **Ion species**: S II, S III, S IV, O II, O III

## Repository Structure

```
IPT_Emission_MOP_Community_Code/
├── cm3/                          # Volume emission rate calculations
│   ├── Python_Code/              # Python implementation
│   ├── IDL_Code/                 # IDL implementation (+ CHIANTI table generation)
│   └── chianti/                  # CHIANTI database (download from Zenodo)
│       ├── dbase/                # Atomic database files
│       └── idl/                  # CHIANTI IDL procedures
├── LOS_Integration/              # Line-of-sight integration
│   ├── Python_Code/              # Python (includes parallelized movie generation)
│   ├── IDL_Code/                 # IDL
│   ├── Fortran_Code/             # Fortran
│   ├── Cpp_Code/                 # C++
│   ├── MATLAB_Code/              # MATLAB
│   ├── Fortran_Code_parallel_make_emission_movie/   # MPI Fortran movie generator
│   └── Cpp_Code_parallel_make_emission_movie/       # MPI C++ movie generator
├── Emiss_tables/                 # CHIANTI emission tables (download from Zenodo)
├── 3D_Torus_Model/               # 3D plasma model (download from Zenodo)
├── CITATION.cff
└── README.md
```

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/egnerney/IPT_Emission_MOP_Community_Code.git
cd IPT_Emission_MOP_Community_Code
```

### 2. Download Data Files from Zenodo

Download the companion dataset from [Zenodo (DOI: 10.5281/zenodo.17808352)](https://doi.org/10.5281/zenodo.17808352) and place files in the appropriate directories:

| File | Location | Required For |
|------|----------|--------------|
| `CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5` | `Emiss_tables/` | All implementations |
| `CHIANTI_11.0.2_emiss_tables_double_maxwellian_16x8x12x10.h5` | `Emiss_tables/` | All implementations |
| `jovian_plasma_interpolated_381x381x231.h5` | `3D_Torus_Model/` | LOS integration |
| `chianti.zip` | `cm3/` (unzip here) | IDL direct CHIANTI calls & table generation |

### 3. Extract CHIANTI Database (IDL users)

Required for direct CHIANTI calculations or generating custom emission tables:

```bash
cd cm3
unzip chianti.zip
```

This creates `cm3/chianti/` containing `dbase/` and `idl/` subdirectories.

### 4. Install Dependencies

**Python:**
```bash
pip install numpy scipy matplotlib h5py
```

**IDL/GDL:** IDL 8.0+ with HDF5 support. For direct CHIANTI calls, set the environment variable:
```bash
export CHIANTI_DATA=/path/to/IPT_Emission_MOP_Community_Code/cm3/chianti/dbase
```

**Fortran:**
```bash
# macOS
brew install gcc hdf5 gnuplot

# Linux
sudo apt install gfortran libhdf5-dev gnuplot
```

**C++:**
```bash
# macOS
brew install gcc hdf5

# Linux
sudo apt install g++ libhdf5-dev libhdf5-cpp-dev
```

## Quick Start

See the README in each subdirectory for language-specific usage. Basic Python example:

```python
from LOS_Integration.Python_Code.IPT_emiss_MOP_community_code import JovianUVEmissionRaytracer

raytracer = JovianUVEmissionRaytracer()
wave, spectrum, lines = raytracer.calculate_spectrum_single(
    [6.0, -20.0, 0.0],  # position [R_J]
    [0.0, 1.0, 0.0],    # direction
    (550, 2100),        # wavelength range [Å]
    1.0, 6.0, 0.01      # bin_width, FWHM, step_size
)
```

## Documentation

Each implementation has its own README with detailed usage instructions:

| Directory | Description |
|-----------|-------------|
| `cm3/Python_Code/` | Python volume emission model |
| `cm3/IDL_Code/` | IDL volume emission model + table generation |
| `LOS_Integration/Python_Code/` | Python LOS raytracer |
| `LOS_Integration/IDL_Code/` | IDL LOS raytracer |
| `LOS_Integration/Fortran_Code/` | Fortran LOS raytracer |
| `LOS_Integration/Cpp_Code/` | C++ LOS raytracer |
| `LOS_Integration/MATLAB_Code/` | MATLAB LOS raytracer |
| `LOS_Integration/Fortran_Code_parallel_make_emission_movie/` | MPI Fortran movie generator |
| `LOS_Integration/Cpp_Code_parallel_make_emission_movie/` | MPI C++ movie generator |

## Citation

If you use this code in your research, please cite both the software and the companion dataset:

### Software

Nerney, E. G. (2025). IPT Emission MOP Community Code (Version 1.0.0) [Software]. Zenodo. https://doi.org/10.5281/zenodo.17808825

```bibtex
@software{nerney_2025_ipt_emission_code,
  author       = {Nerney, Edward G.},
  title        = {{IPT Emission MOP Community Code}},
  month        = dec,
  year         = 2025,
  publisher    = {Zenodo},
  version      = {1.0.0},
  doi          = {10.5281/zenodo.17808825},
  url          = {https://doi.org/10.5281/zenodo.17808825}
}
```

### Companion Dataset

Nerney, E. G. (2025). IPT Emission Model Data Files for MOP Community Code [Data set]. Zenodo. https://doi.org/10.5281/zenodo.17808352

```bibtex
@dataset{nerney_2025_ipt_emission_data,
  author       = {Nerney, Edward G.},
  title        = {{IPT Emission Model Data Files for MOP Community Code}},
  year         = 2025,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.17808352},
  url          = {https://doi.org/10.5281/zenodo.17808352}
}
```

### Related Publications

- Nerney, E. G., Bagenal, F., & Schmidt, C. (2025). Simulations of optical emissions in Io's plasma torus. *Journal of Geophysical Research: Space Physics*, 130, e2024JA033232. https://doi.org/10.1029/2024JA033232

- Nerney, E. G. (2025). Diffusive equilibrium: Modeling anisotropic Maxwellian and kappa field line distributions in Io's plasma torus using multi-fluid and kinetic approaches. *Journal of Geophysical Research: Space Physics*. https://doi.org/10.1029/2024JA033582

- Nerney, E. G., Bagenal, F., & Steffl, A. J. (2017). Io plasma torus ion composition: Voyager, Galileo, and Cassini. *Journal of Geophysical Research: Space Physics*, 122, 727–744. https://doi.org/10.1002/2016JA023306

- Nerney, E. G., & Bagenal, F. (2020). Combining UV spectra and physical chemistry to constrain the hot electron fraction in the Io plasma torus. *Journal of Geophysical Research: Space Physics*, 125, e2019JA027458. https://doi.org/10.1029/2019JA027458

### CHIANTI Database

- Del Zanna, G., Dere, K. P., Young, P. R., & Landi, E. (2021). CHIANTI—An atomic database for emission lines. XVI. Version 10. *The Astrophysical Journal*, 909, 38. https://doi.org/10.3847/1538-4357/abdaf9

## Author

**Edward (Eddie) G. Nerney**  
egnerney@gmail.com

## License

MIT License. See [LICENSE](LICENSE) for details.

## Acknowledgments

- The CHIANTI team for the atomic database
- NASA Juno, Europa Clipper, and ESA JUICE mission teams
- The Magnetospheres of Outer Planets (MOP) community

**CHIANTI Acknowledgment**: CHIANTI is a collaborative project involving George Mason University, the University of Michigan (USA), University of Cambridge (UK), and NASA Goddard Space Flight Center (USA).
