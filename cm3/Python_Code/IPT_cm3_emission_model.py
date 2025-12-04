#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IPT_cm3_emission_model.py

Io Plasma Torus (IPT) UV/Optical Emission Model Library
========================================================

This library provides comprehensive tools for calculating emission spectra from the Io Plasma 
Torus in Jupiter's magnetosphere using pre-calculated CHIANTI atomic database emission tables. 
The code supports both single and double Maxwellian electron distributions and includes 
realistic instrument response modeling for UV and optical wavelength ranges.

PHYSICAL MODEL:
- Pre-calculated volume emission rates from CHIANTI 11.0.2 (photons/s)
- Temperature and density dependent emissivities for single Maxwellian distributions
- Four-dimensional interpolation for double Maxwellian distributions  
- Optically thin plasma approximation (valid for IPT)
- Electron impact excitation with proper atomic physics
- Line-of-sight integration via column densities

APPLICABLE WAVELENGTH RANGES:
- UV instruments: 550-2100 Å (Europa-UVS, JUICE-UVS, HST/STIS)
- Optical instruments: 3000-10000 Å (ground-based telescopes)
- Supports any wavelength range where CHIANTI data is available

COORDINATE SYSTEMS AND UNITS:
- Temperatures: electron volts [eV]
- Densities: particles per cubic centimeter [cm^-3]
- Column densities: particles per square centimeter [cm^-2]
- Wavelengths: Angstroms [Å]
- Emissivities: photons per second per cubic centimeter [photons s^-1 cm^-3]
- Brightnesses: Rayleighs [R], where 1 R = 10^6 photons s^-1 cm^-2 (4π sr)^-1

REFERENCES:
- CHIANTI database: Dere et al. 1997; Del Zanna et al. 2020; Dufresne et al. 2024
- IPT observations: Steffl et al. 2004a,b; Thomas et al. 2004; Bagenal & Delamere 2011
- Emission modeling: Nerney et al. 2017, 2020, 2022, 2025a, 2025b
- Electron distributions: Meyer-Vernet & Moncuquet 1989; Moncuquet et al. 2002

AUTHOR: Edward (Eddie) G. Nerney
INSTITUTION: Laboratory for Atmospheric and Space Physics, University of Colorado Boulder
LICENSE: Open source for academic and research use
VERSION: 1.0
DATE: November 2025

CHIANTI ACKNOWLEDGMENT:
CHIANTI is a collaborative project involving George Mason University, the University 
of Michigan (USA), University of Cambridge (UK), and NASA Goddard Space Flight Center (USA).
"""

import numpy as np
import h5py
from scipy import interpolate
from scipy.special import erf
from scipy.integrate import simpson
import matplotlib.pyplot as plt
from typing import Dict, Tuple, Optional, List
import warnings
from pathlib import Path


class EmissionTables:
    """
    Container for pre-calculated emission tables from CHIANTI atomic database.
    
    This class handles loading and organizing emission data for both single and
    double Maxwellian electron distributions. The tables contain volume emission
    rates (emissivities) as functions of temperature, density, and (for double
    Maxwellian) hot electron fraction.
    
    Attributes
    ----------
    single_maxwellian_loaded : bool
        Flag indicating if single Maxwellian tables are loaded
    double_maxwellian_loaded : bool
        Flag indicating if double Maxwellian tables are loaded
    temp_arr : ndarray
        Temperature grid for single Maxwellian [eV]
    dens_arr : ndarray
        Density grid for single Maxwellian [cm^-3]
    tec_arr : ndarray
        Core temperature grid for double Maxwellian [eV]
    teh_arr : ndarray
        Hot temperature grid for double Maxwellian [eV]
    ne_arr : ndarray
        Density grid for double Maxwellian [cm^-3]
    feh_arr : ndarray
        Hot electron fraction grid for double Maxwellian
    wavelengths_single : dict
        Wavelengths organized by species for single Maxwellian
    emissivities_single : dict
        Volume emission rates organized by species for single Maxwellian
    wavelengths_double : dict
        Wavelengths organized by species for double Maxwellian
    emissivities_double : dict
        Volume emission rates organized by species for double Maxwellian
    
    Notes
    -----
    Species use spectroscopic notation:
    - 'SP': S II (S+)
    - 'S2P': S III (S++)
    - 'S3P': S IV (S+++)
    - 'S4P': S V (S++++)
    - 'OP': O II (O+)
    - 'O2P': O III (O++)
    - 'NAP': Na II (Na+) - trace species in IPT
    """
    
    def __init__(self):
        """Initialize empty emission tables container."""
        self.single_maxwellian_loaded = False
        self.double_maxwellian_loaded = False
        
        # Species names mapping
        self.ion_names = {
            'S': 'S I',      # Neutral sulfur (negligible in IPT)
            'SP': 'S II',    # S+ (singly ionized)
            'S2P': 'S III',  # S++ (doubly ionized) - dominant S ion
            'S3P': 'S IV',   # S+++ (triply ionized)
            'S4P': 'S V',    # S++++ (quadruply ionized)
            'O': 'O I',      # Neutral oxygen (negligible in IPT)
            'OP': 'O II',    # O+ (singly ionized) - dominant ion overall
            'O2P': 'O III',  # O++ (doubly ionized)
            'NAP': 'Na II'   # Na+ (trace species from Io)
        }
    
    def load_single_maxwellian_tables(self, filename: str):
        """
        Load single Maxwellian emission tables from HDF5 file.
        
        Parameters
        ----------
        filename : str
            Path to HDF5 file containing single Maxwellian emission tables
            
        Table Structure
        ---------------
        The HDF5 file contains:
        - T : temperature grid [eV], shape (n_T,)
        - n : density grid [cm^-3], shape (n_n,)
        - emissivity : volume emission rates [photons cm^-3 s^-1], shape (n_T, n_n, n_lines)
        - wavelength : emission line wavelengths [Angstroms], shape (n_lines,)
        - species : species labels, shape (n_lines,)
        
        Notes
        -----
        Emissivities are calculated using CHIANTI atomic database assuming:
        - Coronal equilibrium ionization balance (appropriate for IPT)
        - Electron impact excitation cross sections from CHIANTI
        - Radiative cascade and metastable populations included
        - Photospheric abundances (solar with Io volcanic modifications)
        """
        print("Loading single Maxwellian emission tables...")
        
        with h5py.File(filename, 'r') as f:
            # Load parameter grids
            self.temp_arr = f['T'][:].astype(np.float64)  # Electron Temperature of Single Maxwellian [eV]
            self.dens_arr = f['n'][:].astype(np.float64)  # Electron Number Density of Single Maxwellian [cm^-3]
            
            # Load combined per ion or atom photon emission rate array for each emission line for each species in Photons/s/ion or atom 
            # For Single Maxwellian electron impact excitation from CHIANTI 11.0.2
            # Dimensions:(n_T, n_n, n_lines)
            emissivity_all = f['emiss'][:]
            
            # Load wavelengths and species labels
            wavelength_all = f['wavelength'][:]
            
            # Handle species array (may be bytes or strings depending on HDF5 version)
            species_all = f['species'][:]
            if species_all.dtype.kind == 'S' or species_all.dtype.kind == 'O':
                species_all = np.array([s.decode('utf-8') if isinstance(s, bytes) else str(s) 
                                       for s in species_all])
            
            # Organize by species
            self.wavelengths_single = {}
            self.emissivities_single = {}
            
            # Get unique species
            unique_species = []
            for s in species_all:
                if s not in unique_species:
                    unique_species.append(s)
            
            # Split combined arrays into per-species dictionaries
            for species_key in unique_species:
                # Find indices for this species
                indices = np.where(species_all == species_key)[0]
                
                # Extract wavelengths for this species
                self.wavelengths_single[species_key] = wavelength_all[indices]
                
                # Extract emissivities for this species
                # emissivity_all is (n_T, n_n, n_lines)
                # We reorganize to (n_lines_species, n_T, n_n) for efficient interpolation
                emiss_species = emissivity_all[:, :, indices]  # (n_T, n_n, n_lines_species)
                self.emissivities_single[species_key] = np.transpose(emiss_species, (2, 0, 1))
        
        # Create log-space grids for interpolation
        # Emissivities vary as power laws in T and n, so log interpolation is more accurate
        self.log_temp = np.log10(self.temp_arr)
        self.log_dens = np.log10(self.dens_arr)
        
        self.single_maxwellian_loaded = True
        
        print("Single Maxwellian tables loaded successfully:")
        print(f"  Temperature range: {self.temp_arr.min():.3f} - {self.temp_arr.max():.1f} eV")
        print(f"  Density range: {self.dens_arr.min():.1f} - {self.dens_arr.max():.0f} cm^-3")
        print(f"  Grid size: {len(self.temp_arr)} x {len(self.dens_arr)}")
        print(f"  Species: {', '.join(sorted(self.wavelengths_single.keys()))}")
    
    def load_double_maxwellian_tables(self, filename: str):
            """
            Load double Maxwellian emission tables from HDF5 file.
            
            Parameters
            ----------
            filename : str
                Path to HDF5 file containing double Maxwellian emission tables
                
            Table Structure
            ---------------
            The HDF5 file contains:
            - T_cold : core temperature grid [eV], shape (n_Tc,)
            - T_hot : hot temperature grid [eV], shape (n_Th,)
            - n : total density grid [cm^-3], shape (n_n,)
            - feh : hot electron fraction grid, shape (n_feh,)
            - emissivity : volume emission rates [photons cm^-3 s^-1], 
                          shape (n_Tc, n_Th, n_n, n_feh, n_lines)
            - wavelength : emission line wavelengths [Angstroms], shape (n_lines,)
            - species : species labels, shape (n_lines,)
            
            Physical Model
            --------------
            The double Maxwellian electron distribution function is:
            
            f(v) = (1 - feh) * f_Maxwell(v, Tec) + feh * f_Maxwell(v, Teh)
            
            where:
            - Tec is the core electron temperature
            - Teh is the hot electron temperature
            - feh is the hot electron fraction
            
            This distribution arises from:
            - Wave-particle interactions (electron cyclotron heating)
            - Pickup ion acceleration
            - Magnetic reconnection events
            
            Notes
            -----
            Emission from double Maxwellian is highly nonlinear and cannot be
            approximated by linear superposition of single Maxwellians. High-energy
            transitions are enhanced by hot electrons due to their larger collision
            cross sections at high energies.
            
            Performance
            -----------
            Tables are pre-transposed to (n_Tec, n_Teh, n_ne, n_feh, n_lines) order
            during loading for optimal interpolation performance.
            """
            print("Loading double Maxwellian emission tables...")
            
            with h5py.File(filename, 'r') as f:
                # Load parameter grids
                self.tec_arr = f['T_cold'][:].astype(np.float64)
                self.teh_arr = f['T_hot'][:].astype(np.float64)
                self.ne_arr = f['n'][:].astype(np.float64)
                self.feh_arr = f['feh'][:].astype(np.float64)
                
                # Load combined emission array
                # Original shape: (n_Tc, n_Th, n_n, n_feh, n_lines)
                emissivity_all = f['emiss'][:]
                
                # Load wavelengths and species labels
                wavelength_all = f['wavelength'][:]
                
                # Handle species array
                species_all = f['species'][:]
                if species_all.dtype.kind == 'S' or species_all.dtype.kind == 'O':
                    species_all = np.array([s.decode('utf-8') if isinstance(s, bytes) else str(s) 
                                           for s in species_all])
                
                # Organize by species
                self.wavelengths_double = {}
                self.emissivities_double = {}
                
                # Get unique species preserving order
                unique_species = []
                for s in species_all:
                    if s not in unique_species:
                        unique_species.append(s)
                
                # Split combined arrays into per-species dictionaries
                # Pre-transpose to interpolation-optimal order: (n_lines, n_Tec, n_Teh, n_ne, n_feh)
                for species_key in unique_species:
                    indices = np.where(species_all == species_key)[0]
                    self.wavelengths_double[species_key] = wavelength_all[indices]
                    
                    # Extract and transpose once at load time for optimal memory access
                    # Original: (n_Tc, n_Th, n_n, n_feh, n_lines_species)
                    # Target:   (n_lines_species, n_Tec, n_Teh, n_ne, n_feh)
                    emiss_species = emissivity_all[:, :, :, :, indices]
                    self.emissivities_double[species_key] = np.ascontiguousarray(
                        np.transpose(emiss_species, (4, 0, 1, 2, 3))
                    )
            
            # Create log-space grids for interpolation
            self.log_tec = np.log10(self.tec_arr)
            self.log_teh = np.log10(self.teh_arr)
            self.log_ne = np.log10(self.ne_arr)
            self.log_feh = np.log10(self.feh_arr)
            
            self.double_maxwellian_loaded = True
            
            print("Double Maxwellian tables loaded successfully:")
            print(f"  Core temperature range: {self.tec_arr.min():.3f} - {self.tec_arr.max():.1f} eV")
            print(f"  Hot temperature range: {self.teh_arr.min():.1f} - {self.teh_arr.max():.0f} eV")
            print(f"  Density range: {self.ne_arr.min():.1f} - {self.ne_arr.max():.0f} cm^-3")
            print(f"  Hot fraction range: {self.feh_arr.min():.6f} - {self.feh_arr.max():.4f}")
            print(f"  Grid size: {len(self.tec_arr)} x {len(self.teh_arr)} x {len(self.ne_arr)} x {len(self.feh_arr)}")
            print(f"  Species: {', '.join(sorted(self.wavelengths_double.keys()))}")


def interpolate_emissivity_2d(tables: EmissionTables, 
                              temperature: float, 
                              density: float,
                              species_list: Optional[List[str]] = None) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Interpolate volume emission rates from 2D single Maxwellian tables.
    
    Uses bilinear interpolation in log10(T) - log10(n) parameter space with
    log-space interpolation of emissivity values. This approach accurately
    captures the exponential temperature dependence of atomic emission rates
    near excitation thresholds.
    
    Parameters
    ----------
    tables : EmissionTables
        Object containing loaded emission tables with attributes:
        - temp_arr: 1D array of temperatures [eV]
        - dens_arr: 1D array of densities [cm^-3]
        - log_temp: log10(temp_arr)
        - log_dens: log10(dens_arr)
        - wavelengths_single: dict mapping species to wavelength arrays
        - emissivities_single: dict mapping species to emission rate arrays
                               with shape (n_lines, n_T, n_n)
    temperature : float
        Electron temperature [eV]
    density : float
        Electron density [cm^-3]
    species_list : list of str, optional
        List of species to include. Default: ['SP', 'S2P', 'S3P', 'S4P', 'OP', 'O2P']
        which excludes neutral atoms (negligible in IPT) and trace species.
        
    Returns
    -------
    wavelengths : ndarray
        Emission line wavelengths [Angstroms], sorted in ascending order
    emissivities : ndarray
        Interpolated per-ion emission rates [photons s^-1 per ion]
        (multiply by ion number density to get volume emission rate)
    species_tags : list of str
        Species tag for each emission line, maintaining wavelength association
        
    Physics Notes
    -------------
    Atomic emission rates scale approximately as:
    
        ε ∝ n_e × <σv> ∝ n_e × T^(-1/2) × exp(-ΔE/kT)
    
    where ΔE is the excitation energy. Near excitation thresholds, the
    exponential factor dominates, causing emission to vary by orders of
    magnitude over small temperature ranges.
    
    Log-space interpolation computes the geometric mean of corner values
    rather than the arithmetic mean, which is mathematically correct for
    exponentially-varying functions:
    
        Linear:  f_interp = w*f1 + (1-w)*f0           (arithmetic mean)
        Log:     f_interp = exp(w*ln(f1) + (1-w)*ln(f0))  (geometric mean)
    
    For emission rates that change by factor R between grid points:
        - Linear interpolation error at midpoint: O(|R-1|/(R+1))
        - Log interpolation error at midpoint: ~0 (exact for exponentials)
    
    Implementation
    --------------
    - Grid lookup performed in log10(T), log10(n) space
    - Emissivity interpolation performed in ln(emissivity) space
    - Handles zero/negative values with linear fallback
    - Fully vectorized across all emission lines
    
    References
    ----------
    CHIANTI database: Dere et al. 1997; Del Zanna et al. 2020
    IPT emission modeling: Nerney et al. 2017, 2020, 2022
    
    Raises
    ------
    ValueError
        If tables not loaded
    """
    if not tables.single_maxwellian_loaded:
        raise ValueError("Single Maxwellian tables not loaded. Call load_single_maxwellian_tables() first.")
    
    # Default species list (exclude neutrals and trace species)
    if species_list is None:
        species_list = ['SP', 'S2P', 'S3P', 'S4P', 'OP', 'O2P']
    
    # Convert query point to log space for grid lookup
    log_T = np.log10(temperature)
    log_n = np.log10(density)
    
    # Get grid arrays in log space
    log_T_grid = tables.log_temp
    log_n_grid = tables.log_dens
    
    # Check bounds and issue warnings if extrapolating
    if temperature < tables.temp_arr.min() or temperature > tables.temp_arr.max():
        warnings.warn(f"Temperature {temperature} eV outside table range "
                     f"[{tables.temp_arr.min():.1f}, {tables.temp_arr.max():.1f}] eV. "
                     f"Extrapolation may be inaccurate.")
    
    if density < tables.dens_arr.min() or density > tables.dens_arr.max():
        warnings.warn(f"Density {density} cm^-3 outside table range "
                     f"[{tables.dens_arr.min():.0f}, {tables.dens_arr.max():.0f}] cm^-3. "
                     f"Extrapolation may be inaccurate.")
    
    # Find bracketing indices for temperature in log space
    i_T = np.searchsorted(log_T_grid, log_T)
    i_T = np.clip(i_T, 1, len(log_T_grid) - 1)
    i_T0 = i_T - 1  # Lower temperature index
    i_T1 = i_T      # Upper temperature index
    
    # Find bracketing indices for density in log space
    i_n = np.searchsorted(log_n_grid, log_n)
    i_n = np.clip(i_n, 1, len(log_n_grid) - 1)
    i_n0 = i_n - 1  # Lower density index
    i_n1 = i_n      # Upper density index
    
    # Get bracketing grid values in log space
    log_T0 = log_T_grid[i_T0]
    log_T1 = log_T_grid[i_T1]
    log_n0 = log_n_grid[i_n0]
    log_n1 = log_n_grid[i_n1]
    
    # Compute normalized interpolation weights in [0, 1]
    if log_T1 != log_T0:
        w_T = (log_T - log_T0) / (log_T1 - log_T0)
    else:
        w_T = 0.0
    
    if log_n1 != log_n0:
        w_n = (log_n - log_n0) / (log_n1 - log_n0)
    else:
        w_n = 0.0
    
    # Pre-compute bilinear weight products for the 4 corners
    w00 = (1.0 - w_T) * (1.0 - w_n)  # low T, low n
    w01 = (1.0 - w_T) * w_n          # low T, high n
    w10 = w_T * (1.0 - w_n)          # high T, low n
    w11 = w_T * w_n                  # high T, high n
    
    corner_weights = np.array([w00, w01, w10, w11], dtype=np.float64)
    
    # Minimum emissivity threshold for log interpolation
    # Values below this use linear interpolation to avoid log(0)
    min_emiss = 1.0e-50
    
    # Collect results from all requested species
    wavelengths_list = []
    emissivities_list = []
    species_tags_list = []
    
    for species_key in species_list:
        if species_key not in tables.wavelengths_single:
            warnings.warn(f"Species {species_key} not found in tables. Skipping.")
            continue
        
        wav_species = tables.wavelengths_single[species_key]
        emiss_table = tables.emissivities_single[species_key]
        n_lines = len(wav_species)
        
        if n_lines == 0:
            continue
        
        # Extract emissivities at the four corners of the interpolation cell
        # emiss_table shape: (n_lines, n_T, n_n)
        E_00 = emiss_table[:, i_T0, i_n0]  # (n_lines,) - low T, low n
        E_01 = emiss_table[:, i_T0, i_n1]  # (n_lines,) - low T, high n
        E_10 = emiss_table[:, i_T1, i_n0]  # (n_lines,) - high T, low n
        E_11 = emiss_table[:, i_T1, i_n1]  # (n_lines,) - high T, high n
        
        # Stack corners for vectorized operations: shape (4, n_lines)
        corners = np.vstack([E_00, E_01, E_10, E_11])
        
        # Determine which lines can use log interpolation
        # (all four corners must be positive)
        min_corners = np.min(corners, axis=0)
        use_log = min_corners > min_emiss
        
        # Initialize output array
        emiss_interpolated = np.zeros(n_lines, dtype=np.float64)
        
        # Log-space interpolation for lines with all positive corners
        # Computes geometric mean weighted by corner weights
        if np.any(use_log):
            log_corners = np.log(corners[:, use_log])
            log_interp = np.dot(corner_weights, log_corners)
            emiss_interpolated[use_log] = np.exp(log_interp)
        
        # Linear interpolation fallback for lines with zero/negative corners
        if np.any(~use_log):
            linear_interp = np.dot(corner_weights, corners[:, ~use_log])
            emiss_interpolated[~use_log] = linear_interp
        
        wavelengths_list.append(wav_species)
        emissivities_list.append(emiss_interpolated)
        species_tags_list.extend([species_key] * n_lines)
    
    # Handle case where no valid species were found
    if len(wavelengths_list) == 0:
        return np.array([]), np.array([]), []
    
    # Concatenate all species data
    wavelengths = np.concatenate(wavelengths_list)
    emissivities = np.concatenate(emissivities_list)
    
    # Sort all arrays by wavelength while maintaining species associations
    sort_idx = np.argsort(wavelengths)
    wavelengths_sorted = wavelengths[sort_idx]
    emissivities_sorted = emissivities[sort_idx]
    species_sorted = [species_tags_list[i] for i in sort_idx]
    
    return wavelengths_sorted, emissivities_sorted, species_sorted


def interpolate_emissivity_4d(tables: EmissionTables, 
                              core_temp: float, 
                              hot_temp: float,
                              total_density: float,
                              hot_fraction: float,
                              species_list: Optional[List[str]] = None) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Interpolate volume emission rates from 4D double Maxwellian tables.
    
    Uses quadrilinear interpolation in the 4D parameter space
    (log10(Tec), log10(Teh), log10(ne), log10(feh)) with log-space 
    interpolation of emissivity values. This approach accurately captures
    the exponential temperature dependence of atomic emission rates and
    the nonlinear enhancement from hot electron populations.
    
    Parameters
    ----------
    tables : EmissionTables
        Object containing loaded emission tables with attributes:
        - tec_arr: 1D array of core temperatures [eV]
        - teh_arr: 1D array of hot temperatures [eV]
        - ne_arr: 1D array of total electron densities [cm^-3]
        - feh_arr: 1D array of hot electron fractions (dimensionless)
        - log_tec, log_teh, log_ne, log_feh: log10 of the above arrays
        - wavelengths_double: dict mapping species to wavelength arrays
        - emissivities_double: dict mapping species to emission rate arrays
                               with shape (n_lines, n_feh, n_ne, n_Teh, n_Tec)
    core_temp : float
        Core electron temperature [eV]
    hot_temp : float
        Hot electron temperature [eV]
    total_density : float
        Total electron density [cm^-3]
    hot_fraction : float
        Fraction of electrons in hot population (0 < feh < 1)
    species_list : list of str, optional
        List of species to include. Default: ['SP', 'S2P', 'S3P', 'S4P', 'OP', 'O2P']
        
    Returns
    -------
    wavelengths : ndarray
        Emission line wavelengths [Angstroms], sorted in ascending order
    emissivities : ndarray
        Interpolated per-ion emission rates [photons s^-1 per ion]
    species_tags : list of str
        Species tag for each emission line
        
    Physics Notes
    -------------
    The double Maxwellian electron velocity distribution is:
    
        f(v) = (1 - f_eh) × f_Maxwell(v, T_ec) + f_eh × f_Maxwell(v, T_eh)
    
    where T_ec is the core temperature, T_eh is the hot temperature, and
    f_eh is the hot electron number fraction. This distribution arises from:
    - Wave-particle interactions (electron cyclotron heating)
    - Pickup ion acceleration
    - Magnetic reconnection events
    
    Emission from double Maxwellian distributions is highly nonlinear and
    CANNOT be approximated by linear superposition of two single Maxwellians.
    Hot electrons enhance high-excitation transitions disproportionately
    because collision cross sections increase with energy above threshold.
    
    Atomic emission rates scale approximately as:
    
        ε ∝ n_e × <σv> ∝ n_e × T^(-1/2) × exp(-ΔE/kT)
    
    The exponential factor causes emission to vary by orders of magnitude
    over small temperature ranges, especially near excitation thresholds.
    
    Log-space interpolation computes the geometric mean of corner values:
    
        Linear:  f_interp = Σ w_i × f_i           (arithmetic mean)
        Log:     f_interp = exp(Σ w_i × ln(f_i))  (geometric mean)
    
    This is mathematically correct for exponentially-varying functions and
    reduces interpolation errors from ~10-15% to ~2% for typical IPT grids.
    
    Implementation
    --------------
    - Grid lookup performed in log10 space for all 4 parameters
    - Emissivity interpolation performed in ln(emissivity) space
    - Handles zero/negative values with linear fallback
    - Fully vectorized 16-point quadrilinear interpolation
    - All emission lines processed simultaneously
    
    References
    ----------
    Double Maxwellian distributions: Meyer-Vernet & Moncuquet 1989
    IPT electron distributions: Sittler & Strobel 1987; Moncuquet et al. 2002
    CHIANTI database: Dere et al. 1997; Del Zanna et al. 2020
    IPT emission modeling: Nerney et al. 2017, 2020, 2022, 2025
    
    Raises
    ------
    ValueError
        If tables not loaded
    """
    if not tables.double_maxwellian_loaded:
        raise ValueError("Double Maxwellian tables not loaded. Call load_double_maxwellian_tables() first.")
    
    if species_list is None:
        species_list = ['SP', 'S2P', 'S3P', 'S4P', 'OP', 'O2P']
    
    # Convert parameters to log space for grid lookup
    log_Tec = np.log10(core_temp)
    log_Teh = np.log10(hot_temp)
    log_ne = np.log10(total_density)
    log_feh = np.log10(hot_fraction)
    
    # Get grid arrays in log space
    log_Tec_grid = tables.log_tec
    log_Teh_grid = tables.log_teh
    log_ne_grid = tables.log_ne
    log_feh_grid = tables.log_feh
    
    # Check bounds and issue warnings if extrapolating
    if core_temp < tables.tec_arr[0] or core_temp > tables.tec_arr[-1]:
        warnings.warn(f"Core temperature {core_temp} eV outside table range "
                     f"[{tables.tec_arr[0]:.1f}, {tables.tec_arr[-1]:.1f}] eV.")
    if hot_temp < tables.teh_arr[0] or hot_temp > tables.teh_arr[-1]:
        warnings.warn(f"Hot temperature {hot_temp} eV outside table range "
                     f"[{tables.teh_arr[0]:.0f}, {tables.teh_arr[-1]:.0f}] eV.")
    if total_density < tables.ne_arr[0] or total_density > tables.ne_arr[-1]:
        warnings.warn(f"Density {total_density} cm^-3 outside table range "
                     f"[{tables.ne_arr[0]:.0f}, {tables.ne_arr[-1]:.0f}] cm^-3.")
    if hot_fraction < tables.feh_arr[0] or hot_fraction > tables.feh_arr[-1]:
        warnings.warn(f"Hot fraction {hot_fraction} outside table range "
                     f"[{tables.feh_arr[0]:.6f}, {tables.feh_arr[-1]:.4f}].")
    
    # Find bracketing indices using searchsorted
    # searchsorted returns insertion point; clip ensures valid bracketing
    i_Tec = np.searchsorted(log_Tec_grid, log_Tec).clip(1, len(log_Tec_grid) - 1)
    i_Teh = np.searchsorted(log_Teh_grid, log_Teh).clip(1, len(log_Teh_grid) - 1)
    i_ne = np.searchsorted(log_ne_grid, log_ne).clip(1, len(log_ne_grid) - 1)
    i_feh = np.searchsorted(log_feh_grid, log_feh).clip(1, len(log_feh_grid) - 1)
    
    # Lower and upper indices for each dimension
    i0_Tec, i1_Tec = i_Tec - 1, i_Tec
    i0_Teh, i1_Teh = i_Teh - 1, i_Teh
    i0_ne, i1_ne = i_ne - 1, i_ne
    i0_feh, i1_feh = i_feh - 1, i_feh
    
    # Compute interpolation weights in log space
    # Weight = fractional distance from lower grid point to query point
    dTec = log_Tec_grid[i1_Tec] - log_Tec_grid[i0_Tec]
    dTeh = log_Teh_grid[i1_Teh] - log_Teh_grid[i0_Teh]
    dne = log_ne_grid[i1_ne] - log_ne_grid[i0_ne]
    dfeh = log_feh_grid[i1_feh] - log_feh_grid[i0_feh]
    
    w_Tec = (log_Tec - log_Tec_grid[i0_Tec]) / dTec if dTec != 0 else 0.0
    w_Teh = (log_Teh - log_Teh_grid[i0_Teh]) / dTeh if dTeh != 0 else 0.0
    w_ne = (log_ne - log_ne_grid[i0_ne]) / dne if dne != 0 else 0.0
    w_feh = (log_feh - log_feh_grid[i0_feh]) / dfeh if dfeh != 0 else 0.0
    
    # Pre-compute weight products for all 16 corners of the 4D hypercube
    # Binary indexing convention: bit 0 = Tec, bit 1 = Teh, bit 2 = ne, bit 3 = feh
    # Corner 0000 = (low Tec, low Teh, low ne, low feh)
    # Corner 1111 = (high Tec, high Teh, high ne, high feh)
    w0_Tec, w1_Tec = 1.0 - w_Tec, w_Tec
    w0_Teh, w1_Teh = 1.0 - w_Teh, w_Teh
    w0_ne, w1_ne = 1.0 - w_ne, w_ne
    w0_feh, w1_feh = 1.0 - w_feh, w_feh
    
    corner_weights = np.array([
        w0_Tec * w0_Teh * w0_ne * w0_feh,  # 0000
        w1_Tec * w0_Teh * w0_ne * w0_feh,  # 1000
        w0_Tec * w1_Teh * w0_ne * w0_feh,  # 0100
        w1_Tec * w1_Teh * w0_ne * w0_feh,  # 1100
        w0_Tec * w0_Teh * w1_ne * w0_feh,  # 0010
        w1_Tec * w0_Teh * w1_ne * w0_feh,  # 1010
        w0_Tec * w1_Teh * w1_ne * w0_feh,  # 0110
        w1_Tec * w1_Teh * w1_ne * w0_feh,  # 1110
        w0_Tec * w0_Teh * w0_ne * w1_feh,  # 0001
        w1_Tec * w0_Teh * w0_ne * w1_feh,  # 1001
        w0_Tec * w1_Teh * w0_ne * w1_feh,  # 0101
        w1_Tec * w1_Teh * w0_ne * w1_feh,  # 1101
        w0_Tec * w0_Teh * w1_ne * w1_feh,  # 0011
        w1_Tec * w0_Teh * w1_ne * w1_feh,  # 1011
        w0_Tec * w1_Teh * w1_ne * w1_feh,  # 0111
        w1_Tec * w1_Teh * w1_ne * w1_feh,  # 1111
    ], dtype=np.float64)
    
    # Index arrays for extracting all 16 corners simultaneously
    # These arrays map corner index (0-15) to grid indices
    idx_Tec = np.array([i0_Tec, i1_Tec, i0_Tec, i1_Tec, i0_Tec, i1_Tec, i0_Tec, i1_Tec,
                        i0_Tec, i1_Tec, i0_Tec, i1_Tec, i0_Tec, i1_Tec, i0_Tec, i1_Tec])
    idx_Teh = np.array([i0_Teh, i0_Teh, i1_Teh, i1_Teh, i0_Teh, i0_Teh, i1_Teh, i1_Teh,
                        i0_Teh, i0_Teh, i1_Teh, i1_Teh, i0_Teh, i0_Teh, i1_Teh, i1_Teh])
    idx_ne = np.array([i0_ne, i0_ne, i0_ne, i0_ne, i1_ne, i1_ne, i1_ne, i1_ne,
                       i0_ne, i0_ne, i0_ne, i0_ne, i1_ne, i1_ne, i1_ne, i1_ne])
    idx_feh = np.array([i0_feh, i0_feh, i0_feh, i0_feh, i0_feh, i0_feh, i0_feh, i0_feh,
                        i1_feh, i1_feh, i1_feh, i1_feh, i1_feh, i1_feh, i1_feh, i1_feh])
    
    # Minimum emissivity threshold for log interpolation
    # Values below this use linear interpolation to avoid log(0)
    min_emiss = 1.0e-50
    
    # Process all species with vectorized operations
    wavelengths_list = []
    emissivities_list = []
    species_tags_list = []
    
    for species_key in species_list:
        if species_key not in tables.wavelengths_double:
            warnings.warn(f"Species {species_key} not found in double Maxwellian tables. Skipping.")
            continue
        
        wav_species = tables.wavelengths_double[species_key]
        emiss_table = tables.emissivities_double[species_key]
        n_lines = len(wav_species)
        
        if n_lines == 0:
            continue
        
        # Extract all 16 corners at once using advanced indexing
        # emiss_table shape: (n_lines, n_feh, n_ne, n_Teh, n_Tec)
        # Advanced indexing extracts values at all 16 corner combinations
        # Result shape after transpose: (16, n_lines)
        corners = emiss_table[:, idx_Tec, idx_Teh, idx_ne, idx_feh].T
        
        # Determine which lines can use log interpolation
        # Requires all 16 corners to be positive (non-zero emission)
        min_corners = np.min(corners, axis=0)
        use_log = min_corners > min_emiss
        
        # Initialize output array
        emiss_interpolated = np.zeros(n_lines, dtype=np.float64)
        
        # Log-space interpolation for lines with all positive corners
        # This computes the weighted geometric mean: exp(Σ w_i × ln(f_i))
        if np.any(use_log):
            log_corners = np.log(corners[:, use_log])
            log_interp = np.dot(corner_weights, log_corners)
            emiss_interpolated[use_log] = np.exp(log_interp)
        
        # Linear interpolation fallback for lines with zero/negative corners
        # Uses standard weighted arithmetic mean: Σ w_i × f_i
        if np.any(~use_log):
            linear_interp = np.dot(corner_weights, corners[:, ~use_log])
            emiss_interpolated[~use_log] = linear_interp
        
        wavelengths_list.append(wav_species)
        emissivities_list.append(emiss_interpolated)
        species_tags_list.extend([species_key] * n_lines)
    
    if len(wavelengths_list) == 0:
        return np.array([]), np.array([]), []
    
    # Concatenate and sort by wavelength
    wavelengths = np.concatenate(wavelengths_list)
    emissivities = np.concatenate(emissivities_list)
    
    sort_idx = np.argsort(wavelengths)
    
    return wavelengths[sort_idx], emissivities[sort_idx], [species_tags_list[i] for i in sort_idx]

def calculate_ipt_emiss_tables_single(tables: EmissionTables,
                                     temperature: float,
                                     density: float,
                                     column_densities: Dict[str, float],
                                     min_wav: float = 550.0,
                                     max_wav: float = 2100.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate IPT emission line brightnesses for single Maxwellian electron distribution.
    
    Converts volume emission rates (emissivities) to observable line brightnesses
    by multiplying by ion column densities along the line of sight.
    
    Parameters
    ----------
    tables : EmissionTables
        Object containing loaded emission tables
    temperature : float
        Electron temperature [eV]
    density : float
        Electron density [cm^-3]
    column_densities : dict
        Ion column densities [cm^-2] with keys:
        - 'S+': S II (singly ionized sulfur)
        - 'S++': S III (doubly ionized sulfur)
        - 'S+++': S IV (triply ionized sulfur)
        - 'S++++': S V (quadruply ionized sulfur)
        - 'O+': O II (singly ionized oxygen)
        - 'O++': O III (doubly ionized oxygen)
    min_wav : float, optional
        Minimum wavelength to include [Angstroms], default=550
    max_wav : float, optional
        Maximum wavelength to include [Angstroms], default=2100
        
    Returns
    -------
    wavelengths : ndarray
        Emission line wavelengths within specified range [Angstroms]
    brightnesses : ndarray
        Line brightnesses [Rayleighs]
        
    Physics
    -------
    The brightness of an emission line in Rayleighs is:
    
    I = ε(Te, ne) * N_ion * 10^-6  [Rayleighs]
    
    where:
    - ε(Te, ne) is the volume emissivity [photons s^-1 cm^-3]
    - N_ion is the ion column density [cm^-2]
    - 10^-6 converts to Rayleighs (1 R = 10^6 photons s^-1 cm^-2 (4π sr)^-1)
    
    This assumes:
    - Optically thin plasma (no self-absorption)
    - Uniform temperature and density along line of sight
    - Ion column density approximates line-of-sight integral
    
    Notes
    -----
    Typical IPT column densities for ~6 R_J path length:
    - O+: 5 × 10^13 cm^-2 (dominant ion)
    - S++: 4 × 10^13 cm^-2 (dominant sulfur ion)
    - S+: 1 × 10^13 cm^-2
    - O++, S+++: ~6 × 10^12 cm^-2
    - S++++: ~6 × 10^11 cm^-2
    
    See Steffl et al. (2004b), Nerney et al. (2017) for typical values.
    """
    print("Calculating single Maxwellian emission using tables...")
    
    # Interpolate emissivities at specified temperature and density
    # Returns wavelengths, emissivities, and species tags
    wavelengths, emissivities, species_tags = interpolate_emissivity_2d(
        tables, temperature, density
    )
    
    # Map user-friendly notation to internal species keys
    species_columns = {
        'SP': column_densities.get('S+', 0.0),
        'S2P': column_densities.get('S++', 0.0),
        'S3P': column_densities.get('S+++', 0.0),
        'S4P': column_densities.get('S++++', 0.0),
        'OP': column_densities.get('O+', 0.0),
        'O2P': column_densities.get('O++', 0.0)
    }
    
    # Calculate brightnesses by multiplying emissivities by column densities
    # Species tags ensure each line is multiplied by correct column density
    brightnesses = np.zeros_like(emissivities)
    for i, species_key in enumerate(species_tags):
        if species_key in species_columns:
            # Convert from photons s^-1 cm^-3 to Rayleighs
            # 1 Rayleigh = 10^6 photons s^-1 cm^-2 (4π sr)^-1
            brightnesses[i] = emissivities[i] * species_columns[species_key] * 1e-6
    
    # Filter to requested wavelength range
    mask = (wavelengths >= min_wav) & (wavelengths <= max_wav)
    
    return wavelengths[mask], brightnesses[mask]


def calculate_ipt_emiss_tables_double(tables: EmissionTables,
                                     core_temp: float,
                                     hot_temp: float,
                                     density: float,
                                     hot_fraction: float,
                                     column_densities: Dict[str, float],
                                     min_wav: float = 550.0,
                                     max_wav: float = 2100.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate IPT emission line brightnesses for double Maxwellian electron distribution.
    
    Accounts for suprathermal electron population from wave heating, which significantly
    enhances high-excitation transitions.
    
    Parameters
    ----------
    tables : EmissionTables
        Object containing loaded emission tables
    core_temp : float
        Core electron temperature [eV]
    hot_temp : float
        Hot electron temperature [eV]
    density : float
        Total electron density [cm^-3]
    hot_fraction : float
        Fraction of electrons in hot population (0 < feh < 1)
    column_densities : dict
        Ion column densities [cm^-2] with keys:
        - 'S+': S II, 'S++': S III, 'S+++': S IV, 'S++++': S V
        - 'O+': O II, 'O++': O III
    min_wav : float, optional
        Minimum wavelength [Angstroms], default=550
    max_wav : float, optional
        Maximum wavelength [Angstroms], default=2100
        
    Returns
    -------
    wavelengths : ndarray
        Emission line wavelengths within specified range [Angstroms]
    brightnesses : ndarray
        Line brightnesses [Rayleighs]
        
    Physics
    -------
    The double Maxwellian distribution represents:
    
    f(v) = (1 - feh) * f_Maxwell(v, Tec) + feh * f_Maxwell(v, Teh)
    
    High-excitation lines (e.g., S IV, S V) are strongly enhanced by even
    small fractions of hot electrons due to their exponentially larger
    collision cross sections at high energies.
    """
    print("Calculating double Maxwellian emission using 4D tables...")
    
    # Interpolate emissivities using 4D interpolation
    wavelengths, emissivities, species_tags = interpolate_emissivity_4d(
        tables, core_temp, hot_temp, density, hot_fraction
    )
    
    # Build column density lookup array for vectorized brightness calculation
    species_to_column = {
        'SP': column_densities.get('S+', 0.0),
        'S2P': column_densities.get('S++', 0.0),
        'S3P': column_densities.get('S+++', 0.0),
        'S4P': column_densities.get('S++++', 0.0),
        'OP': column_densities.get('O+', 0.0),
        'O2P': column_densities.get('O++', 0.0)
    }
    
    # Vectorized brightness calculation using numpy array operations
    column_array = np.array([species_to_column.get(s, 0.0) for s in species_tags], 
                            dtype=np.float64)
    brightnesses = emissivities * column_array * 1e-6
    
    # Filter to wavelength range using boolean indexing
    mask = (wavelengths >= min_wav) & (wavelengths <= max_wav)
    
    return wavelengths[mask], brightnesses[mask]


def simulate_ipt_spectrum_rayleighs_erf_form(wavelength_grid: np.ndarray,
                                            bin_width: float,
                                            line_wavelengths: np.ndarray,
                                            line_brightnesses: np.ndarray,
                                            fwhm: float = 6.0) -> np.ndarray:
    """
    Convolve discrete emission lines with instrument response function.
    
    Uses Error Function (ERF) formulation for exact integration of Gaussian
    line profiles over finite wavelength bins. This is more accurate than
    simple Gaussian evaluation at bin centers, especially for wide bins or
    narrow lines.
    
    Parameters
    ----------
    wavelength_grid : ndarray
        Output wavelength grid (bin centers) [Angstroms]
    bin_width : float
        Width of each wavelength bin [Angstroms]
    line_wavelengths : ndarray
        Discrete emission line wavelengths [Angstroms]
    line_brightnesses : ndarray
        Discrete emission line brightnesses [Rayleighs]
    fwhm : float, optional
        Full width at half maximum of instrument response [Angstroms], default=6.0
        (appropriate for HST/STIS low resolution and JUICE/Europa-UVS)
        
    Returns
    -------
    spectrum : ndarray
        Convolved spectrum [Rayleighs/Angstrom] at each wavelength grid point
        
    Method
    ------
    For each emission line, calculates the integral of a Gaussian line profile
    over each wavelength bin using the error function:
    
    I_bin = I_line * (1/2) * [erf((λ_max - λ_line)/(σ√2)) - erf((λ_min - λ_line)/(σ√2))]
    
    where:
    - σ = FWHM / (2√(2 ln 2)) ≈ FWHM / 2.355 is the Gaussian standard deviation
    - λ_min, λ_max are the bin edges (center ± bin_width/2)
    - erf(x) is the error function: erf(x) = (2/√π) ∫[0,x] exp(-t²) dt
    - The factor of 1/2 comes from the normalization of the Gaussian CDF
    
    This exactly integrates the Gaussian profile rather than sampling at the bin center,
    which is critical when:
    - Bin width is comparable to or larger than FWHM
    - Line falls near bin edge rather than center
    - High accuracy is needed for photometric measurements
    
    Physical Interpretation
    -----------------------
    The erf terms represent the cumulative distribution function (CDF) of the Gaussian:
    - erf((x - μ)/(σ√2)) gives the integral of the normalized Gaussian from -∞ to x
    - The difference [erf(b) - erf(a)] gives the fraction of total line brightness
      that falls within the wavelength interval [a, b]
    - Dividing by bin_width converts from integrated brightness [Rayleighs] 
      to brightness per unit wavelength [Rayleighs/Angstrom]
    
    Instrument Response
    -------------------
    Typical FWHM values for various spectrographs:
    - HST/STIS G140M: ~0.6 Å (~40 km/s at 1200 Å, R~1000)
    - HST/STIS G230M: ~1.0 Å (~130 km/s at 2300 Å, R~2000)
    - HST/STIS G430L: ~6 Å (low resolution, R~500)
    - Ground-based echelle: ~0.05-0.1 Å (R~50,000-100,000)
    - Europa-UVS/JUICE-UVS: ~6 Å (best case, goal for extended source)
    - New Horizons Alice: ~8-12 Å (far-UV)
    
    Implementation Notes
    --------------------
    - Fully vectorized using NumPy broadcasting for optimal performance
    - All emission lines processed simultaneously without Python loops
    - Memory efficient: only creates arrays of shape (n_lines, n_bins)
    - The Gaussian approximation is appropriate for most UV/optical spectrographs
    - For high-resolution work, actual instrument line spread functions (LSF)
      may be asymmetric and require more sophisticated convolution
    
    Examples
    --------
    >>> # Create wavelength grid
    >>> wavelengths = np.linspace(1000, 2000, 1000)
    >>> bin_width = 1.0  # 1 Angstrom bins
    >>> 
    >>> # Define emission lines
    >>> line_wavs = np.array([1215.67, 1304.35, 1355.6])  # Ly-alpha, O I, O I
    >>> line_brights = np.array([100.0, 50.0, 25.0])  # Rayleighs
    >>> 
    >>> # Convolve with 6 Angstrom FWHM instrument response
    >>> spectrum = simulate_ipt_spectrum_rayleighs_erf_form(
    ...     wavelengths, bin_width, line_wavs, line_brights, fwhm=6.0)
    """
    # Convert FWHM to Gaussian standard deviation
    # Relationship: FWHM = 2σ√(2 ln 2) ≈ 2.355σ
    # Therefore: σ = FWHM / (2√(2 ln 2))
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    
    # Pre-compute normalization factor for ERF argument
    # This is 1/(σ√2), used in the ERF argument: (x - μ)/(σ√2)
    sigma_sqrt2 = sigma * np.sqrt(2.0)
    
    # Calculate bin edges from bin centers
    # Each bin extends from (center - width/2) to (center + width/2)
    bin_edges_low = wavelength_grid - bin_width / 2.0
    bin_edges_high = wavelength_grid + bin_width / 2.0
    
    # Check for empty line lists
    if len(line_wavelengths) == 0 or len(line_brightnesses) == 0:
        return np.zeros(len(wavelength_grid))
    
    # Reshape line arrays for broadcasting
    # Shape: (n_lines, 1) so they broadcast with bin_edges shape (n_bins,)
    # This creates arrays of shape (n_lines, n_bins) without explicit loops
    line_wav = line_wavelengths[:, np.newaxis]  # Column vector
    line_int = line_brightnesses[:, np.newaxis]  # Column vector
    
    # Calculate ERF arguments for all line-bin combinations simultaneously
    # For the lower bin edge:
    # Shape: (n_lines, n_bins) via broadcasting
    # Each element is (bin_edge_low[j] - line_wavelength[i]) / (σ√2)
    erf_arg_low = (bin_edges_low - line_wav) / sigma_sqrt2
    
    # For the upper bin edge:
    # Shape: (n_lines, n_bins)
    # Each element is (bin_edge_high[j] - line_wavelength[i]) / (σ√2)
    erf_arg_high = (bin_edges_high - line_wav) / sigma_sqrt2
    
    # Evaluate error function at all points
    # erf(x) is the integral of normalized Gaussian from 0 to x (scaled)
    # The difference [erf(high) - erf(low)] gives the fraction of the
    # Gaussian that falls within [low, high]
    erf_low = erf(erf_arg_low)
    erf_high = erf(erf_arg_high)
    
    # Calculate the fraction of each line's brightness in each bin
    # The factor of 1/2 comes from the normalization:
    # erf represents integral from -∞ to x, scaled to [-1, 1] range
    # The 1/2 converts this to the [0, 1] probability range
    # Shape: (n_lines, n_bins)
    bin_contributions = 0.5 * (erf_high - erf_low)
    
    # Multiply by line brightnesses to get absolute contributions
    # Shape: (n_lines, n_bins)
    weighted_contributions = line_int * bin_contributions
    
    # Sum contributions from all lines for each bin
    # Shape: (n_bins,)
    spectrum_integrated = np.sum(weighted_contributions, axis=0)
    
    # Convert from integrated brightness [Rayleighs] to brightness per wavelength
    # [Rayleighs/Angstrom] by dividing by bin width
    spectrum = spectrum_integrated / bin_width
    
    return spectrum


def analyze_emission_line(wavelength_array, spectrum_array, target_wavelength, 
                         fwhm_instrument):
    """
    Analyze a specific emission line and calculate integrated brightness.
    
    This function extracts a spectral region around a target wavelength,
    identifies the peak emission, and integrates the line brightness over
    ±3σ from the peak (approximately the FWHM region).
    
    Parameters
    ----------
    wavelength_array : np.ndarray
        Wavelength grid in Angstroms
    spectrum_array : np.ndarray
        Spectral brightness in Rayleighs/Angstrom
    target_wavelength : float
        Target wavelength for the emission line in Angstroms
    fwhm_instrument : float
        Instrumental FWHM in Angstroms for convolved spectrum
    
    Returns
    -------
    dict
        Dictionary containing:
        - wavelength: wavelength array for extracted region (centered on peak)
        - spectrum: spectrum array for extracted region
        - peak_wavelength: wavelength of peak brightness
        - peak_brightness: maximum brightness value in R/Å
        - fwhm: instrumental FWHM in Å
        - sigma: Standard deviation of Gaussian at given FWHM (FWHM/2.35482) in Å
        - integration_bounds: tuple of (lower, upper) wavelength bounds for integration
          These are the actual min/max wavelengths of the extracted region
        - integrated_brightness: integrated line brightness in Rayleighs
    
    Notes
    -----
    Integration is performed using Simpson's rule over approximately ±3σ from 
    the peak, which captures ~99.7% of a Gaussian line profile. This approach 
    is standard for UV/optical spectroscopy of planetary magnetospheres.
    
    The spectral extraction window is centered on the identified peak wavelength
    (not the nominal target wavelength) to ensure consistency between the plotted
    region and integration bounds. Peak identification prioritizes finding the
    maximum brightness within ±3σ of the target wavelength to avoid confusion
    with nearby emission lines in crowded spectral regions.
    
    The integration bounds reported are the actual wavelength range of the
    extracted discrete data points, ensuring exact correspondence between
    plotted spectral region and integration limits for visualization consistency.
    """
    # Calculate sigma from instrumental FWHM
    # For Gaussian profiles: FWHM = 2.35482 * sigma
    sigma = fwhm_instrument / 2.35482
    
    # Define search window to locate peak emission
    # Use ±6σ from target wavelength to be generous with initial search
    search_half_width = 6.0 * sigma
    search_mask = (wavelength_array >= target_wavelength - search_half_width) & \
                  (wavelength_array <= target_wavelength + search_half_width)
    
    wl_search = wavelength_array[search_mask]
    spec_search = spectrum_array[search_mask]
    
    if len(wl_search) == 0:
        raise ValueError(f"No data found near {target_wavelength} Å")
    
    # Find peak that is both high in brightness AND close to target wavelength
    # Prioritize peaks within ±3σ of target to avoid confusion with nearby lines
    # This is critical for crowded spectral regions (e.g., [S II] 4068/4076 doublet)
    proximity_limit = 3.0 * sigma
    
    # Create mask for points within ±3σ of target wavelength
    proximity_mask = np.abs(wl_search - target_wavelength) <= proximity_limit
    
    if np.any(proximity_mask):
        # Find the maximum brightness within the proximity region
        # This ensures we select the correct line in crowded spectra
        wl_proximity = wl_search[proximity_mask]
        spec_proximity = spec_search[proximity_mask]
        peak_idx_proximity = np.argmax(spec_proximity)
        peak_wl = wl_proximity[peak_idx_proximity]
        peak_brightness = spec_proximity[peak_idx_proximity]
    else:
        # Fallback: if no data within ±3σ (rare edge case with very coarse binning)
        # Use the point closest to target wavelength
        closest_idx = np.argmin(np.abs(wl_search - target_wavelength))
        peak_wl = wl_search[closest_idx]
        peak_brightness = spec_search[closest_idx]
    
    # Define theoretical integration bounds: ±3σ from peak
    # This captures ~99.7% of Gaussian line flux
    theoretical_wl_min = peak_wl - 3.0 * sigma
    theoretical_wl_max = peak_wl + 3.0 * sigma
    
    # Extract spectral region centered on peak wavelength with ±3σ extent
    # This ensures plotted region matches integration bounds exactly
    region_mask = (wavelength_array >= theoretical_wl_min) & \
                  (wavelength_array <= theoretical_wl_max)
    
    wl_region = wavelength_array[region_mask]
    spec_region = spectrum_array[region_mask]
    
    if len(wl_region) == 0:
        raise ValueError(f"No data found in integration region [{theoretical_wl_min:.3f}, {theoretical_wl_max:.3f}] Å")
    
    # Integrate using Simpson's rule over the entire extracted region
    # Since region is already ±3σ from peak, integrate all points
    # Converts from R/Å × Å = R (total line brightness)
    integrated_brightness = simpson(spec_region, x=wl_region)
    
    # Report actual integration bounds as the min/max of extracted wavelengths
    # This ensures exact correspondence between plotted data and shaded region
    # for visualization consistency across different wavelength grid spacings
    actual_wl_min = wl_region[0]
    actual_wl_max = wl_region[-1]
    
    # Package results
    results = {
        'wavelength': wl_region,
        'spectrum': spec_region,
        'peak_wavelength': peak_wl,
        'peak_brightness': peak_brightness,
        'fwhm': fwhm_instrument,
        'sigma': sigma,
        'integration_bounds': (actual_wl_min, actual_wl_max),
        'integrated_brightness': integrated_brightness
    }
    
    return results


def plot_emission_line_detail(line_data_single, line_data_double, line_name,
                              ion_species, wavelength_nominal, filename):
    """
    Create detailed publication-quality plot of a specific emission line.
    
    Generates a comparison plot showing single vs double Maxwellian emission
    for a specific spectral line, with integration region highlighted and
    quantitative results annotated.
    
    Parameters
    ----------
    line_data_single : dict
        Analysis results from analyze_emission_line() for single Maxwellian
    line_data_double : dict
        Analysis results from analyze_emission_line() for double Maxwellian
    line_name : str
        Descriptive name for the emission line (e.g., "S++ 680Å")
    ion_species : str
        Ion species notation (e.g., "S++", "O+")
    wavelength_nominal : float
        Nominal/rest wavelength of the transition in Angstroms
    filename : str
        Output filename for saved figure
    
    Returns
    -------
    None
        Saves figure to disk
    
    Notes
    -----
    Plot includes:
    - Single and double Maxwellian spectra
    - Peak wavelength indicator
    - Integration region (±3σ) shading
    - Quantitative summary of integrated brightnesses
    - Brightness enhancement factor (double/single)
    """
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    # Plot both spectra
    ax.plot(line_data_single['wavelength'], line_data_single['spectrum'], 
            "C0",linewidth=2, label='Single Maxwellian', alpha=0.7)
    ax.plot(line_data_double['wavelength'], line_data_double['spectrum'], 
            "C1",linewidth=2, label='Double Maxwellian', alpha=0.7)
    
    # Mark peak wavelength
    ax.axvline(line_data_double['peak_wavelength'], color='gray', 
               linestyle='--', linewidth=1.5, alpha=0.6, 
               label=f'Peak: {line_data_double["peak_wavelength"]:.3f} Å')
    
    # Highlight integration region (±3σ)
    ax.axvspan(line_data_double['integration_bounds'][0], 
               line_data_double['integration_bounds'][1], 
               alpha=0.15, color='green', label='Integration region (±3σ)')
    
    # Create annotation with quantitative results
    brightness_ratio = (line_data_double['integrated_brightness'] / 
                       line_data_single['integrated_brightness'])
    
    annotation_text = (
        f'{line_name}\n'
        f'Nominal λ: {wavelength_nominal:.2f} Å\n'
        f'Peak λ: {line_data_double["peak_wavelength"]:.3f} Å\n'
        f'Peak brightness: {line_data_double["peak_brightness"]:.2f} R/Å\n'
        f'FWHM: {line_data_double["fwhm"]:.3f} Å (σ = {line_data_double["sigma"]:.3f} Å)\n'
        f'\n'
        f'Integrated Brightness (±3σ):\n'
        f'  Single Maxwellian: {line_data_single["integrated_brightness"]:.1f} R\n'
        f'  Double Maxwellian: {line_data_double["integrated_brightness"]:.1f} R\n'
        f'  Enhancement Factor: {brightness_ratio:.3f}'
    )
    
    ax.text(0.98, 0.98, annotation_text, transform=ax.transAxes,
            fontsize=10, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8, 
                     edgecolor='black', linewidth=1))
    
    # Labels and formatting
    ax.set_xlabel('Wavelength [Å]', fontsize=12, fontweight='bold')
    ax.set_ylabel('Brightness [R/Å]', fontsize=12, fontweight='bold')
    ax.set_title(f'Io Plasma Torus Emission: {ion_species} at {wavelength_nominal:.1f} Å', 
                 fontsize=13, fontweight='bold')
    ax.legend(loc='upper left', fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f'  Saved: {filename}')
    plt.show()


def plot_emission_lines_summary(line_results, filename):
    """
    Create multi-panel summary plot of analyzed emission lines.
    
    Generates a publication-quality figure with separate panels for each
    analyzed emission line, comparing single and double Maxwellian models.
    
    Parameters
    ----------
    line_results : dict
        Dictionary with line names as keys and nested dicts containing
        'single' and 'double' analysis results from analyze_emission_line()
    filename : str
        Output filename for saved figure
    
    Returns
    -------
    None
        Saves figure to disk
    
    Notes
    -----
    This function creates a stacked subplot layout with consistent formatting
    across all panels for easy comparison between different emission lines.
    Commonly used for comparing sulfur and oxygen ion emissions in the IPT.
    """
    n_lines = len(line_results)
    fig, axes = plt.subplots(n_lines, 1, figsize=(12, 5*n_lines))
    
    # Handle single subplot case
    if n_lines == 1:
        axes = [axes]
    
    fig.suptitle('Io Plasma Torus Emission: Detailed Line Analysis', 
                 fontsize=14, fontweight='bold')
    
    for idx, (line_name, line_data) in enumerate(line_results.items()):
        ax = axes[idx]
        
        single_data = line_data['single']
        double_data = line_data['double']
        
        # Plot spectra
        ax.plot(single_data['wavelength'], single_data['spectrum'],"C0", linewidth=1.5, label='Single Maxwellian', alpha=0.7)
        ax.plot(double_data['wavelength'], double_data['spectrum'],"C1", linewidth=1.5, label='Double Maxwellian', alpha=0.7)
        
        # Mark peak and integration region
        ax.axvline(double_data['peak_wavelength'], color='gray', 
                   linestyle='--', linewidth=1, alpha=0.5)
        ax.axvspan(double_data['integration_bounds'][0], 
                   double_data['integration_bounds'][1], 
                   alpha=0.2, color='green', label='Integration (±3σ)')
        
        # Annotation with key results
        brightness_ratio = (double_data['integrated_brightness'] / 
                           single_data['integrated_brightness'])
        annotation_text = (
            f'Integrated Brightness:\n'
            f'Single: {single_data["integrated_brightness"]:.1f} R\n'
            f'Double: {double_data["integrated_brightness"]:.1f} R\n'
            f'Ratio: {brightness_ratio:.3f}\n'
            f'FWHM: {double_data["fwhm"]:.3f} Å'
        )
        ax.text(0.02, 0.98, annotation_text, transform=ax.transAxes,
                fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
        
        ax.set_xlabel('Wavelength [Å]', fontsize=11)
        ax.set_ylabel('Brightness [R/Å]', fontsize=11)
        ax.set_title(f'{line_name}', fontsize=12, fontweight='bold')
        ax.legend(loc='upper right', fontsize=9)
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f'  Saved: {filename}')
    plt.show()


def analyze_key_uv_emission_lines(wavelength_output, spectrum_single, spectrum_double,
                               fwhm_instrument=6.0):
    """
    Analyze key UV/optical emission lines in the Io Plasma Torus spectrum.
    
    This function performs detailed analysis of scientifically important
    emission lines commonly observed in IPT spectra, including sulfur
    and oxygen ion transitions. Calculates integrated brightnesses and
    generates publication-quality diagnostic plots.
    
    Parameters
    ----------
    wavelength_output : np.ndarray
        Output wavelength grid in Angstroms
    spectrum_single : np.ndarray
        Single Maxwellian spectrum in R/Å
    spectrum_double : np.ndarray
        Double Maxwellian spectrum in R/Å
    fwhm_instrument : float, optional
        Instrumental FWHM in Angstroms after convolution (default: 0.6 Å)
        Typical for space-based UV spectrographs like HST/STIS, MAVEN/IUVS
    
    Returns
    -------
    dict
        Nested dictionary with line names as keys, containing 'single' and
        'double' sub-dictionaries with analysis results from analyze_emission_line()
    
    Notes
    -----
    Lines analyzed:
    - S++ 680.9 Å: Important diagnostic for S++ density and temperature
    - O+ 833.7 Å: Key tracer of neutral oxygen ionization in the torus
    
    These lines are frequently used in studies of:
    - Io's volcanic activity and atmospheric escape
    - Plasma source rates and composition
    - Energy deposition from magnetospheric electrons
    - Comparison with in-situ plasma measurements
    
    References
    ----------
    Thomas, N. et al. (2004), The Io plasma torus, In Jupiter: The Planet,
        Satellites and Magnetosphere, Cambridge University Press
    Steffl, A. J. et al. (2004), Cassini UVIS observations of the Io plasma
        torus, Icarus, 172, 91-103
    Nerney, E.G. et al. (2017), (2020), (2025a), & (2025b)
    """
    print('\n' + '='*64)
    print('DETAILED EMISSION LINE ANALYSIS')
    print('='*64)
    print('Analyzing key diagnostic lines for IPT composition and energetics')
    print(f'Instrumental FWHM: {fwhm_instrument:.3f} Å')
    print('Integration method: Simpson\'s rule over ±3σ (99.7% of line flux)')
    
    # Define emission lines of scientific interest
    # Wavelengths are vacuum wavelengths from NIST/CHIANTI atomic database
    emission_lines = {
        'S++ 680.9 Å': {
            'wavelength': 680.9,
            'ion': 'S++'
        },
        'O+ 833.7 Å': {
            'wavelength': 833.7,
            'ion': 'O+'
        }
    }
    
    # Store results for all lines
    line_results = {}
    
    for line_name, line_info in emission_lines.items():
        print(f'\nAnalyzing {line_name}...')
        
        # Analyze single Maxwellian spectrum
        single_result = analyze_emission_line(
            wavelength_output, spectrum_single,
            line_info['wavelength'], fwhm_instrument
        )
        
        # Analyze double Maxwellian spectrum
        double_result = analyze_emission_line(
            wavelength_output, spectrum_double,
            line_info['wavelength'], fwhm_instrument
        )
        
        # Calculate brightness enhancement from hot electron population
        enhancement = double_result['integrated_brightness'] / single_result['integrated_brightness']
        
        # Store results
        line_results[line_name] = {
            'single': single_result,
            'double': double_result,
            'nominal_wavelength': line_info['wavelength'],
            'ion': line_info['ion']
        }
        
        # Print quantitative results
        print(f'  Peak wavelength: {double_result["peak_wavelength"]:.3f} Å')
        print(f'  Peak brightness (double): {double_result["peak_brightness"]:.2f} R/Å')
        print(f'  Integration bounds: [{double_result["integration_bounds"][0]:.3f}, '
              f'{double_result["integration_bounds"][1]:.3f}] Å')
        print(f'  Integrated brightness (single): {single_result["integrated_brightness"]:.1f} R')
        print(f'  Integrated brightness (double): {double_result["integrated_brightness"]:.1f} R')
        print(f'  Hot electron enhancement factor: {enhancement:.3f}')
    
    # Generate plots
    print('\nGenerating emission line diagnostic plots...')
    
    # Create summary plot with all lines
    plot_emission_lines_summary(line_results, 'ipt_uv_emission_lines_summary.png')
    
    # Create individual detailed plots for each line
    for line_name, line_data in line_results.items():
        # Create safe filename
        safe_filename = line_name.replace(' ', '_').replace('+', 'p').replace('.', 'p')
        filename = f'ipt_uv_emission_{safe_filename}.png'
        
        plot_emission_line_detail(
            line_data['single'], line_data['double'],
            line_name, line_data['ion'],
            line_data['nominal_wavelength'], filename
        )
    
    # Print summary table
    print('\n' + '='*64)
    print('EMISSION LINE ANALYSIS SUMMARY')
    print('='*64)
    print(f'{"Line":<20} {"Peak [Å]":<12} {"Single [R]":<15} {"Double [R]":<15} {"Ratio":<8}')
    print('-'*64)
    for line_name, line_data in line_results.items():
        single_int = line_data['single']['integrated_brightness']
        double_int = line_data['double']['integrated_brightness']
        ratio = double_int / single_int
        peak = line_data['double']['peak_wavelength']
        print(f'{line_name:<20} {peak:<12.3f} {single_int:<15.1f} {double_int:<15.1f} {ratio:<8.3f}')
    print('='*64)
    
    return line_results

def analyze_key_optical_emission_lines(wavelength_output, spectrum_single, spectrum_double,
                                       fwhm_instrument=3.0):
    """
    Analyze key optical emission lines in the Io Plasma Torus spectrum.
    
    This function performs detailed analysis of scientifically important
    optical emission lines commonly observed in ground-based IPT spectra,
    including forbidden transitions from sulfur and oxygen ions. Calculates
    integrated brightnesses and generates publication-quality diagnostic plots.
    
    Parameters
    ----------
    wavelength_output : np.ndarray
        Output wavelength grid in Angstroms
    spectrum_single : np.ndarray
        Single Maxwellian spectrum in R/Å
    spectrum_double : np.ndarray
        Double Maxwellian spectrum in R/Å
    fwhm_instrument : float, optional
        Instrumental FWHM in Angstroms after convolution (default: 3.0 Å)
        Typical for ground-based optical spectrographs with R~3000
    
    Returns
    -------
    dict
        Nested dictionary with line names as keys, containing 'single' and
        'double' sub-dictionaries with analysis results from analyze_emission_line()
    
    Notes
    -----
    Lines analyzed include major diagnostic transitions for IPT composition:
    
    Singly Ionized Sulfur [S II]:
    - 6716.4 Å: Forbidden transition, density diagnostic doublet component
    - 6730.8 Å: Forbidden transition, forms density-sensitive doublet with 6716
    - 4068.6 Å: Auroral transition in blue spectral region
    - 4076.3 Å: Auroral transition, forms doublet with 4069
    
    Doubly Ionized Sulfur [S III]:
    - 6312.1 Å: Forbidden transition, temperature diagnostic
    - 9068.6 Å: Near-infrared forbidden transition
    
    Singly Ionized Oxygen [O II]:
    - 3726.0 Å: Forbidden doublet component, density diagnostic
    - 3728.8 Å: Forbidden doublet component, forms ratio with 3726
    
    Doubly Ionized Oxygen [O III]:
    - 4363.2 Å: Auroral line, temperature diagnostic
    - 4958.9 Å: Nebular doublet component, strong emission
    - 5006.8 Å: Nebular doublet component, brightest [O III] line
    
    These optical lines are frequently used in studies of:
    - Torus electron density and temperature diagnostics
    - Ion composition and mixing ratios
    - Ground-based monitoring of volcanic activity at Io
    - Comparison with space-based UV observations
    - Energy budget and radiative cooling processes
    
    The [S II] 6716/6731 ratio is density-sensitive and widely used
    for IPT electron density measurements from ground-based telescopes.
    The [O II] 3726/3729 doublet provides complementary density diagnostics.
    The [O III] 4363/(4959+5007) ratio is temperature-sensitive.
    
    References
    ----------
    Brown, M. E. (1994), Observation of mass loading in the Io plasma torus,
        GRL, 21, 847-850
    Brown, M. E. (1995), Periodicities in the Io plasma torus, JGR, 100, 21683
    Küppers, M. & Schneider, N. M. (2000), Discovery of chlorine in the
        Io torus, GRL, 27, 513
    Thomas, N. et al. (2004), The Io plasma torus, In Jupiter: The Planet,
        Satellites and Magnetosphere, Cambridge University Press
    Steffl, A. J. et al. (2004), Cassini UVIS observations of the Io plasma
        torus, Icarus, 172, 91-103
    Oliversen, R. J. et al. (2001), Sunlit Io atmospheric [O I] 6300 Å
        emission and the plasma torus, JGR, 106, 26183
    Roth, L. et al. (2014), Orbital apocenter is not a sufficient condition
        for HST/STIS detection of Europa's water vapor aurora, PNAS, 111, E5123
    """
    print('\n' + '='*64)
    print('DETAILED OPTICAL EMISSION LINE ANALYSIS')
    print('='*64)
    print('Analyzing key diagnostic lines for IPT optical spectroscopy')
    print(f'Instrumental FWHM: {fwhm_instrument:.3f} Å')
    print('Integration method: Simpson\'s rule over ±3σ (99.7% of line flux)')
    
    # Define optical emission lines of scientific interest for IPT diagnostics
    # Wavelengths are vacuum wavelengths from NIST/CHIANTI atomic database
    # Lines are ordered by wavelength for systematic analysis
    emission_lines = {
        # Singly Ionized Oxygen [O II] - density diagnostic doublet
        # These are the classic nebular lines from O+ (dominant ion in IPT)
        '[O II] 3726.0 Å': {
            'wavelength': 3726.03,
            'ion': 'O+'
        },
        '[O II] 3728.8 Å': {
            'wavelength': 3728.82,
            'ion': 'O+'
        },
        # Singly Ionized Sulfur [S II] - auroral lines in blue
        '[S II] 4068.6 Å': {
            'wavelength': 4068.6,
            'ion': 'S+'
        },
        '[S II] 4076.3 Å': {
            'wavelength': 4076.3,
            'ion': 'S+'
        },
        # Doubly Ionized Oxygen [O III] - temperature and density diagnostics
        # 4363 Å is the auroral line (temperature diagnostic)
        '[O III] 4363.2 Å': {
            'wavelength': 4363.21,
            'ion': 'O++'
        },
        # [O III] nebular doublet - strong forbidden lines
        '[O III] 4958.9 Å': {
            'wavelength': 4958.91,
            'ion': 'O++'
        },
        '[O III] 5006.8 Å': {
            'wavelength': 5006.84,
            'ion': 'O++'
        },
        # Doubly Ionized Sulfur [S III] - temperature diagnostic
        '[S III] 6312.1 Å': {
            'wavelength': 6312.1,
            'ion': 'S++'
        },
        # Singly Ionized Sulfur [S II] - density diagnostic doublet (red lines)
        '[S II] 6716.4 Å': {
            'wavelength': 6716.4,
            'ion': 'S+'
        },
        '[S II] 6730.8 Å': {
            'wavelength': 6730.8,
            'ion': 'S+'
        },
        # Doubly Ionized Sulfur [S III] - near-infrared forbidden line
        '[S III] 9068.6 Å': {
            'wavelength': 9068.6,
            'ion': 'S++'
        }
    }
    
    # Store results for all lines
    line_results = {}
    
    for line_name, line_info in emission_lines.items():
        print(f'\nAnalyzing {line_name}...')
        
        try:
            # Analyze single Maxwellian spectrum
            single_result = analyze_emission_line(
                wavelength_output, spectrum_single,
                line_info['wavelength'], fwhm_instrument
            )
            
            # Analyze double Maxwellian spectrum
            double_result = analyze_emission_line(
                wavelength_output, spectrum_double,
                line_info['wavelength'], fwhm_instrument
            )
            
            # Calculate brightness enhancement from hot electron population
            # This ratio indicates sensitivity to suprathermal electrons
            enhancement = double_result['integrated_brightness'] / single_result['integrated_brightness']
            
            # Store results in structured dictionary for further analysis
            line_results[line_name] = {
                'single': single_result,
                'double': double_result,
                'nominal_wavelength': line_info['wavelength'],
                'ion': line_info['ion']
            }
            
            # Print quantitative results for each line
            print(f'  Peak wavelength: {double_result["peak_wavelength"]:.3f} Å')
            print(f'  Peak brightness (double): {double_result["peak_brightness"]:.2f} R/Å')
            print(f'  Integration bounds: [{double_result["integration_bounds"][0]:.3f}, '
                  f'{double_result["integration_bounds"][1]:.3f}] Å')
            print(f'  Integrated brightness (single): {single_result["integrated_brightness"]:.1f} R')
            print(f'  Integrated brightness (double): {double_result["integrated_brightness"]:.1f} R')
            print(f'  Hot electron enhancement factor: {enhancement:.3f}')
            
        except ValueError as e:
            # Handle cases where emission line is too weak or outside wavelength range
            print(f'  Warning: Could not analyze {line_name} - {str(e)}')
            continue
    
    # Generate diagnostic plots if we have results
    if line_results:
        print('\nGenerating optical emission line diagnostic plots...')
        
        # Create summary plot showing all analyzed lines in one figure
        # Useful for comparing relative brightness and enhancement factors
        plot_emission_lines_summary(line_results, 'ipt_optical_emission_lines_summary.png')
        
        # Create individual detailed plots for each emission line
        # These are publication-quality figures with full quantitative annotations
        for line_name, line_data in line_results.items():
            # Create safe filename by replacing special characters
            # Removes brackets, spaces, and converts + and . to alphanumeric
            safe_filename = line_name.replace(' ', '_').replace('+', 'p').replace('.', 'p').replace('[', '').replace(']', '')
            filename = f'ipt_optical_emission_{safe_filename}.png'
            
            plot_emission_line_detail(
                line_data['single'], line_data['double'],
                line_name, line_data['ion'],
                line_data['nominal_wavelength'], filename
            )
        
        # Print summary table for quick reference
        # Shows integrated brightnesses and hot electron enhancement ratios
        print('\n' + '='*64)
        print('OPTICAL EMISSION LINE ANALYSIS SUMMARY')
        print('='*64)
        print(f'{"Line":<20} {"Peak [Å]":<12} {"Single [R]":<15} {"Double [R]":<15} {"Ratio":<8}')
        print('-'*64)
        for line_name, line_data in line_results.items():
            single_int = line_data['single']['integrated_brightness']
            double_int = line_data['double']['integrated_brightness']
            ratio = double_int / single_int
            peak = line_data['double']['peak_wavelength']
            print(f'{line_name:<20} {peak:<12.3f} {single_int:<15.1f} {double_int:<15.1f} {ratio:<8.3f}')
        print('='*64)
    
    return line_results