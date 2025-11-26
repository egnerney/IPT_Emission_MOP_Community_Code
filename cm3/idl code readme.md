IPT CM3 EMISSION MODEL IDL IMPLEMENTATION
==========================================
MOP Community Code v1.0
November 2025

FILE STRUCTURE
--------------
1. IPT_cm3_emission_model.pro
   - Main module containing all functions
   - Load emission tables from HDF5 files
   - Calculate emission line brightnesses
   - Perform 2D (single Maxwellian) and 4D (double Maxwellian) interpolation
   - Convolve discrete lines with instrument response
   - Analyze key diagnostic emission lines

2. basic_example_uv_cm3_emission_model_use_tables.pro
   - Example for UV emission (550-2100 Angstroms)
   - Demonstrates single and double Maxwellian calculations
   - Creates publication-quality plots
   - Analyzes key UV diagnostic lines

3. basic_example_optical_cm3_emission_model_use_tables.pro
   - Example for optical emission (3000-10000 Angstroms)
   - Demonstrates single and double Maxwellian calculations
   - Creates publication-quality plots with zoom regions
   - Analyzes key optical diagnostic lines

USAGE
-----
1. Ensure emission table HDF5 files are in ../Emiss_tables/ directory:
   - CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5
   - CHIANTI_11.0.2_emiss_tables_double_maxwellian_15x10x20x10.h5

2. Run UV example:
   IDL> .compile basic_example_uv_cm3_emission_model_use_tables
   IDL> basic_example_uv_cm3_emission_model_use_tables

3. Run optical example:
   IDL> .compile basic_example_optical_cm3_emission_model_use_tables
   IDL> basic_example_optical_cm3_emission_model_use_tables

KEY FEATURES
------------
- Species-by-species interpolation for accurate emission calculations
- Support for single and double Maxwellian electron distributions
- Error function convolution for realistic instrument response
- Matplotlib-compatible color scheme for consistent plotting
- Comprehensive line analysis and enhancement factor calculations
- Publication-ready plot generation

PLASMA PARAMETERS
-----------------
Single Maxwellian:
- Te: Electron temperature [eV]
- ne: Electron density [cm^-3]

Double Maxwellian:
- Tec: Core electron temperature [eV]
- Teh: Hot electron temperature [eV]
- ne_total: Total electron density [cm^-3]
- feh: Hot electron fraction

Column Densities:
- S+, S++, S+++, S++++: Sulfur ions [cm^-2]
- O+, O++: Oxygen ions [cm^-2]

OUTPUT
------
- Wavelength arrays with emission line positions
- Brightness arrays in Rayleighs
- Convolved spectra in R/Angstrom
- PNG plots of emission spectra
- Line analysis with enhancement factors

INTERPOLATION METHOD
--------------------
Single Maxwellian (2D):
- Bilinear interpolation in log(Te) and log(ne) space
- Species-by-species processing

Double Maxwellian (4D):
- Sequential linear interpolation through each dimension
- Order: Tec -> Teh -> ne -> feh
- Species-by-species processing

NOTES
-----
- All functions are fully self-contained with no external dependencies
- Code optimized for the Magnetospheres of Outer Planets (MOP) community
- Based on CHIANTI 11.0.2 atomic database calculations
- Emission coefficients in units of erg cm^3/s
- Brightness output in Rayleighs (10^6 photons/cm^2/s/sr)