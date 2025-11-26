#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IPT_emiss_MOP_community_code.py

Io Plasma Torus (IPT) UV/Optical Emission Line-of-Sight Integration Module
==========================================================================

This module provides comprehensive functionality for calculating UV and optical emission
spectra from the Io Plasma Torus (IPT) using ray tracing through a 3D plasma model with
CHIANTI atomic emission tables. Supports both single and double Maxwellian electron
distributions with proper species-by-species interpolation and vectorized computations.

PHYSICAL MODEL:
- Plasma torus extends from ~5-10 R_J in cylindrical radius
- Peak emission near 6 R_J at centrifugal equator
- Scale height ~0.5-1 R_J
- Based on IPT Isotropic Model interpolated to rectilinear grid
- Per-ion photon emission rates from CHIANTI 11.0.2 atomic database
- Line-of-sight integration using Simpson's rule
- Optically thin plasma approximation (valid for IPT)

COMPUTATIONAL APPROACH:
- Ray tracing through 3D plasma model with trilinear interpolation
- Bilinear interpolation of single Maxwellian emission tables (2D: Te, ne)
- Quadrilinear interpolation of double Maxwellian emission tables (4D: Tec, Teh, ne, feh)
- Vectorized emission rate interpolation for all lines simultaneously
- Per-ion photon emission rates multiplied by ion density and 1e-6 Rayleigh factor
- Simpson's rule integration over line of sight
- Analytic ERF-based Gaussian convolution for instrument response

COORDINATE SYSTEMS AND UNITS:
- Positions: Jupiter radii [R_J]
- Temperatures: electron volts [eV]
- Densities: particles per cubic centimeter [cm^-3]
- Wavelengths: Angstroms [Å]
- Emission rates: photons per second per ion [photons s^-1 ion^-1]
- Brightnesses: Rayleighs [R], where 1 R = 10^6 photons s^-1 cm^-2 (4π sr)^-1

APPLICABLE WAVELENGTH RANGES:
- UV instruments: 550-2100 Å (JUICE-UVS, Europa-UVS, HST/STIS)
- Optical instruments: 3000-10000 Å (ground-based telescopes, HST optical)

REFERENCES:
- CHIANTI database: Dere et al. 1997; Del Zanna et al. 2020; Dufresne et al. 2024
- IPT observations: Steffl et al. 2004a,b; Thomas et al. 2004; Bagenal & Delamere 2011
- Emission modeling: Nerney et al. 2017, 2020, 2022, 2025a, 2025b
- Electron distributions: Meyer-Vernet & Moncuquet 1989; Moncuquet et al. 2002

AUTHOR: Edward (Eddie) G. Nerney
INSTITUTION: Laboratory for Atmospheric and Space Physics, University of Colorado Boulder
LICENSE: Open source for academic and research use
VERSION: 2.0
DATE: November 2025

CHIANTI ACKNOWLEDGMENT:
CHIANTI is a collaborative project involving George Mason University, the University 
of Michigan (USA), University of Cambridge (UK), and NASA Goddard Space Flight Center (USA).
"""

import numpy as np
import h5py
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import simpson
from scipy.special import erf
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
from typing import Dict, Tuple, Optional, List

# ==============================================================================
# PHYSICAL CONSTANTS
# ==============================================================================

R_J = 71492.0  # Jupiter radius in km
R_J_CM = 7.1492e9  # Jupiter radius in cm
RAYLEIGH_FACTOR = 1e-6 * R_J_CM  # Conversion factor for integration to Rayleighs


# ==============================================================================
# MAIN RAYTRACER CLASS
# ==============================================================================

class JovianUVEmissionRaytracer:
    """
    Main class for calculating UV and optical emission through the Jovian plasma torus.
    
    This class provides a complete workflow for line-of-sight integrated emission
    calculations using ray tracing through a 3D plasma model combined with CHIANTI
    atomic emission tables. Supports both single and double Maxwellian electron
    distributions with optimized vectorized computations.
    
    Attributes
    ----------
    x_axis, y_axis, z_axis : ndarray
        Coordinate axes of the 3D plasma model in Jupiter radii
    nec, neh : ndarray
        Cold and hot electron density fields [cm^-3]
    nsp, ns2p, ns3p, nop, no2p : ndarray
        Ion density fields for S+, S++, S+++, O+, O++ [cm^-3]
    Tec, Teh : ndarray
        Cold and hot electron temperature fields [eV]
    temp_arr, dens_arr : ndarray
        Temperature and density grids for single Maxwellian tables
    tec_arr, teh_arr, ne_arr, feh_arr : ndarray
        Parameter grids for double Maxwellian tables
    wavelengths_single, emissivities_single : dict
        Emission data organized by species for single Maxwellian
    wavelengths_double, emissivities_double : dict
        Emission data organized by species for double Maxwellian
        
    Methods
    -------
    load_plasma_model(filename)
        Load 3D plasma model from HDF5 file
    load_emission_tables_single(filename)
        Load single Maxwellian emission tables
    load_emission_tables_double(filename)
        Load double Maxwellian emission tables
    trace_ray(start_pos, direction, ds, max_distance)
        Trace ray through plasma model
    calculate_spectrum_single(slit_pos_vec, norm_vec, ...)
        Calculate LOS-integrated spectrum for single Maxwellian
    calculate_spectrum_double(slit_pos_vec, norm_vec, ...)
        Calculate LOS-integrated spectrum for double Maxwellian
    """
    
    def __init__(self, plasma_file=None, emission_file_single=None, emission_file_double=None):
        """
        Initialize the raytracer with plasma model and emission tables.
        
        Parameters
        ----------
        plasma_file : str or Path, optional
            Path to the plasma model HDF5 file. If None, uses default location.
        emission_file_single : str or Path, optional
            Path to single Maxwellian emission tables HDF5 file. If None, uses default.
        emission_file_double : str or Path, optional
            Path to double Maxwellian emission tables HDF5 file. If None, uses default.
        """
        print("Initializing Jovian UV/Optical Emission Raytracer...")
        
        # Initialize state flags
        self.single_maxwellian_loaded = False
        self.double_maxwellian_loaded = False
        
        # Species names mapping for display
        self.ion_names = {
            'SP': 'S II',    # S+ (singly ionized sulfur)
            'S2P': 'S III',  # S++ (doubly ionized sulfur) - dominant S ion
            'S3P': 'S IV',   # S+++ (triply ionized sulfur)
            'S4P': 'S V',    # S++++ (quadruply ionized sulfur)
            'OP': 'O II',    # O+ (singly ionized oxygen) - dominant ion overall
            'O2P': 'O III',  # O++ (doubly ionized oxygen)
        }
        
        # Default species list for calculations (excludes trace species)
        self.default_species = ['SP', 'S2P', 'S3P', 'OP', 'O2P']
        
        # Set up default paths if not provided
        if plasma_file is None:
            plasma_file = self._get_default_plasma_path()
        if emission_file_single is None:
            emission_file_single = self._get_default_emission_path_single()
        if emission_file_double is None:
            emission_file_double = self._get_default_emission_path_double()
        
        # Load all data files
        self.load_plasma_model(plasma_file)
        self.load_emission_tables_single(emission_file_single)
        self.load_emission_tables_double(emission_file_double)
        
        print("Initialization complete.")
    
    def _get_default_plasma_path(self):
        """Get the default path to the plasma model file."""
        module_dir = Path(__file__).parent
        plasma_dir = module_dir.parent / "3D_Torus_Model"
        plasma_file = plasma_dir / "jovian_plasma_interpolated_381x381x231.h5"
        
        if not plasma_file.exists():
            raise FileNotFoundError(
                f"Plasma model not found at expected location: {plasma_file}\n"
                f"Please ensure the file is in the 3D_Torus_Model directory."
            )
        
        return plasma_file
    
    def _get_default_emission_path_single(self):
        """Get the default path to the single Maxwellian emission tables file."""
        module_dir = Path(__file__).parent
        emiss_dir = module_dir.parent / "Emiss_tables"
        emiss_file = emiss_dir / "CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5"
        
        if not emiss_file.exists():
            raise FileNotFoundError(
                f"Single Maxwellian emission tables not found at: {emiss_file}\n"
                f"Please ensure the file is in the Emiss_tables directory."
            )
        
        return emiss_file
    
    def _get_default_emission_path_double(self):
        """Get the default path to the double Maxwellian emission tables file."""
        module_dir = Path(__file__).parent
        emiss_dir = module_dir.parent / "Emiss_tables"
        emiss_file = emiss_dir / "CHIANTI_11.0.2_emiss_tables_double_maxwellian_24x10x24x12.h5"
        
        if not emiss_file.exists():
            warnings.warn(
                f"Double Maxwellian emission tables not found at: {emiss_file}\n"
                f"Only single Maxwellian calculations will be available."
            )
            return None
        
        return emiss_file
    
    # ==========================================================================
    # DATA LOADING METHODS
    # ==========================================================================
    
    def load_plasma_model(self, filename):
        """
        Load the 3D plasma model from HDF5 file.
        
        The model contains a non-uniform rectilinear Cartesian grid (381x381x231)
        with density fields for various ion species [cm^-3] and temperature fields [eV]
        in Jupiter System III coordinates (co-rotating with Jupiter's magnetic field).
        
        Parameters
        ----------
        filename : str or Path
            Path to the plasma model HDF5 file
            
        Notes
        -----
        The plasma model is based on Bagenal & Delamere (2011) and related works,
        with ion densities derived from UV spectroscopic observations (Cassini UVIS,
        HST/STIS) and Voyager/Galileo in-situ measurements.
        """
        print(f"Loading plasma model from {filename}...")
        
        with h5py.File(filename, 'r') as f:
            # Load coordinate axes (in Jupiter radii)
            self.x_axis = f['coordinates/x'][:].astype(np.float64)
            self.y_axis = f['coordinates/y'][:].astype(np.float64)
            self.z_axis = f['coordinates/z'][:].astype(np.float64)
            
            # Load electron density fields [cm^-3]
            self.nec = f['data/ne_c'][:].astype(np.float64)  # Cold/core electron density
            self.neh = f['data/ne_h'][:].astype(np.float64)  # Hot electron density
            
            # Load ion density fields [cm^-3]
            self.nsp = f['data/nsp'][:].astype(np.float64)    # S+ density
            self.ns2p = f['data/ns2p'][:].astype(np.float64)  # S++ density
            self.ns3p = f['data/ns3p'][:].astype(np.float64)  # S+++ density
            self.nop = f['data/nop'][:].astype(np.float64)    # O+ density
            self.no2p = f['data/no2p'][:].astype(np.float64)  # O++ density
            
            # Load electron temperature fields [eV]
            self.Tec = f['data/Te_c'][:].astype(np.float64)  # Cold electron temperature
            self.Teh = f['data/Te_h'][:].astype(np.float64)  # Hot electron temperature
        
        # Calculate derived quantities
        self.ne_total = self.nec + self.neh
        
        # Calculate hot electron fraction with safe division
        with np.errstate(divide='ignore', invalid='ignore'):
            self.feh = np.where(self.ne_total > 0, self.neh / self.ne_total, 0.0)
        
        # Clean up invalid values (NaN, negative, non-finite -> 0)
        for field in [self.nec, self.neh, self.ne_total, self.feh, 
                      self.nsp, self.ns2p, self.ns3p, self.nop, self.no2p, 
                      self.Tec, self.Teh]:
            mask = ~np.isfinite(field) | (field < 0)
            field[mask] = 0.0
        
        # Create 3D interpolators for all plasma fields
        print("Creating 3D interpolators for plasma parameters...")
        self._create_plasma_interpolators()
        
        print(f"Plasma model loaded: grid shape {self.nec.shape}")
        print(f"  X range: {self.x_axis.min():.1f} to {self.x_axis.max():.1f} R_J")
        print(f"  Y range: {self.y_axis.min():.1f} to {self.y_axis.max():.1f} R_J")
        print(f"  Z range: {self.z_axis.min():.1f} to {self.z_axis.max():.1f} R_J")
    
    def _create_plasma_interpolators(self):
        """
        Create 3D interpolators for all plasma parameters.
        
        Uses scipy's RegularGridInterpolator for trilinear interpolation
        on non-uniform rectilinear grids. Points outside the grid return 0.
        """
        # Common interpolator arguments
        interp_kwargs = dict(method='linear', bounds_error=False, fill_value=0.0)
        grid_points = (self.x_axis, self.y_axis, self.z_axis)
        
        # Electron density interpolators
        self.interp_nec = RegularGridInterpolator(grid_points, self.nec, **interp_kwargs)
        self.interp_neh = RegularGridInterpolator(grid_points, self.neh, **interp_kwargs)
        self.interp_ne_total = RegularGridInterpolator(grid_points, self.ne_total, **interp_kwargs)
        self.interp_feh = RegularGridInterpolator(grid_points, self.feh, **interp_kwargs)
        
        # Ion density interpolators
        self.interp_nsp = RegularGridInterpolator(grid_points, self.nsp, **interp_kwargs)
        self.interp_ns2p = RegularGridInterpolator(grid_points, self.ns2p, **interp_kwargs)
        self.interp_ns3p = RegularGridInterpolator(grid_points, self.ns3p, **interp_kwargs)
        self.interp_nop = RegularGridInterpolator(grid_points, self.nop, **interp_kwargs)
        self.interp_no2p = RegularGridInterpolator(grid_points, self.no2p, **interp_kwargs)
        
        # Temperature interpolators
        self.interp_Tec = RegularGridInterpolator(grid_points, self.Tec, **interp_kwargs)
        self.interp_Teh = RegularGridInterpolator(grid_points, self.Teh, **interp_kwargs)
        
        # Store ion interpolator mapping for convenient access
        self._ion_interp_map = {
            'SP': self.interp_nsp,
            'S2P': self.interp_ns2p,
            'S3P': self.interp_ns3p,
            'OP': self.interp_nop,
            'O2P': self.interp_no2p,
        }
    
    def load_emission_tables_single(self, filename):
        """
        Load single Maxwellian emission tables from HDF5 file.
        
        Parameters
        ----------
        filename : str or Path
            Path to HDF5 file containing single Maxwellian emission tables
            
        Table Structure
        ---------------
        The HDF5 file contains:
        - T : temperature grid [eV], shape (n_T,)
        - n : density grid [cm^-3], shape (n_n,)
        - emiss : per-ion photon emission rates [photons s^-1 ion^-1],
                  shape (n_T, n_n, n_lines)
        - wavelength : emission line wavelengths [Angstroms], shape (n_lines,)
        - species : species labels, shape (n_lines,)
        
        Notes
        -----
        Emission rates are stored as per-ion values (not volume emission rates).
        To get volume emission rate: multiply by ion density [cm^-3].
        To get brightness: multiply by ion column density [cm^-2] and 1e-6.
        """
        print(f"Loading single Maxwellian emission tables from {filename}...")
        
        with h5py.File(filename, 'r') as f:
            # Load parameter grids
            self.temp_arr = f['T'][:].astype(np.float64)
            self.dens_arr = f['n'][:].astype(np.float64)
            
            # Load combined emission array: (n_T, n_n, n_lines)
            emissivity_all = f['emiss'][:]
            
            # Load wavelengths and species labels
            wavelength_all = f['wavelength'][:]
            
            # Handle species array (may be bytes or strings)
            species_all = f['species'][:]
            if species_all.dtype.kind == 'S' or species_all.dtype.kind == 'O':
                species_all = np.array([s.decode('utf-8') if isinstance(s, bytes) else str(s) 
                                       for s in species_all])
        
        # Create log-space grids for interpolation (power-law behavior)
        self.log_temp = np.log10(self.temp_arr)
        self.log_dens = np.log10(self.dens_arr)
        
        # Organize data by species for efficient access
        self.wavelengths_single = {}
        self.emissivities_single = {}
        
        # Get unique species preserving order
        unique_species = []
        for s in species_all:
            if s not in unique_species:
                unique_species.append(s)
        
        # Split combined arrays into per-species dictionaries
        for species_key in unique_species:
            indices = np.where(species_all == species_key)[0]
            self.wavelengths_single[species_key] = wavelength_all[indices].astype(np.float64)
            
            # Reorganize to (n_lines_species, n_T, n_n) for vectorized interpolation
            emiss_species = emissivity_all[:, :, indices]
            self.emissivities_single[species_key] = np.ascontiguousarray(
                np.transpose(emiss_species, (2, 0, 1)).astype(np.float64)
            )
        
        self.single_maxwellian_loaded = True
        
        print("Single Maxwellian tables loaded successfully:")
        print(f"  Temperature range: {self.temp_arr.min():.3f} - {self.temp_arr.max():.1f} eV")
        print(f"  Density range: {self.dens_arr.min():.1f} - {self.dens_arr.max():.0f} cm^-3")
        print(f"  Grid size: {len(self.temp_arr)} x {len(self.dens_arr)}")
        print(f"  Species: {', '.join(sorted(self.wavelengths_single.keys()))}")
    
    def load_emission_tables_double(self, filename):
        """
        Load double Maxwellian emission tables from HDF5 file.
        
        Parameters
        ----------
        filename : str or Path or None
            Path to HDF5 file containing double Maxwellian emission tables.
            If None, double Maxwellian calculations will be unavailable.
            
        Table Structure
        ---------------
        The HDF5 file contains:
        - T_cold : core temperature grid [eV], shape (n_Tc,)
        - T_hot : hot temperature grid [eV], shape (n_Th,)
        - n : total density grid [cm^-3], shape (n_n,)
        - feh : hot electron fraction grid, shape (n_feh,)
        - emiss : per-ion photon emission rates [photons s^-1 ion^-1],
                  shape (n_Tc, n_Th, n_n, n_feh, n_lines)
        - wavelength : emission line wavelengths [Angstroms], shape (n_lines,)
        - species : species labels, shape (n_lines,)
        
        Physical Model
        --------------
        The double Maxwellian electron distribution function is:
        
        f(v) = (1 - feh) * f_Maxwell(v, Tec) + feh * f_Maxwell(v, Teh)
        
        High-excitation lines are strongly enhanced by hot electrons due to their
        larger collision cross sections at high energies.
        """
        if filename is None:
            self.double_maxwellian_loaded = False
            return
        
        print(f"Loading double Maxwellian emission tables from {filename}...")
        
        with h5py.File(filename, 'r') as f:
            # Load parameter grids
            self.tec_arr = f['T_cold'][:].astype(np.float64)
            self.teh_arr = f['T_hot'][:].astype(np.float64)
            self.ne_arr = f['n'][:].astype(np.float64)
            self.feh_arr = f['feh'][:].astype(np.float64)
            
            # Load combined emission array: (n_Tc, n_Th, n_n, n_feh, n_lines)
            emissivity_all = f['emiss'][:]
            
            # Load wavelengths and species labels
            wavelength_all = f['wavelength'][:]
            
            # Handle species array
            species_all = f['species'][:]
            if species_all.dtype.kind == 'S' or species_all.dtype.kind == 'O':
                species_all = np.array([s.decode('utf-8') if isinstance(s, bytes) else str(s) 
                                       for s in species_all])
        
        # Create log-space grids for interpolation
        self.log_tec = np.log10(self.tec_arr)
        self.log_teh = np.log10(self.teh_arr)
        self.log_ne = np.log10(self.ne_arr)
        self.log_feh = np.log10(self.feh_arr)
        
        # Organize data by species
        self.wavelengths_double = {}
        self.emissivities_double = {}
        
        # Get unique species preserving order
        unique_species = []
        for s in species_all:
            if s not in unique_species:
                unique_species.append(s)
        
        # Split combined arrays into per-species dictionaries
        # Pre-transpose to optimal order: (n_lines, n_Tec, n_Teh, n_ne, n_feh)
        for species_key in unique_species:
            indices = np.where(species_all == species_key)[0]
            self.wavelengths_double[species_key] = wavelength_all[indices].astype(np.float64)
            
            # Original: (n_Tc, n_Th, n_n, n_feh, n_lines_species)
            # Target:   (n_lines_species, n_Tec, n_Teh, n_ne, n_feh)
            emiss_species = emissivity_all[:, :, :, :, indices]
            self.emissivities_double[species_key] = np.ascontiguousarray(
                np.transpose(emiss_species, (4, 0, 1, 2, 3)).astype(np.float64)
            )
        
        self.double_maxwellian_loaded = True
        
        print("Double Maxwellian tables loaded successfully:")
        print(f"  Core temperature range: {self.tec_arr.min():.3f} - {self.tec_arr.max():.1f} eV")
        print(f"  Hot temperature range: {self.teh_arr.min():.1f} - {self.teh_arr.max():.0f} eV")
        print(f"  Density range: {self.ne_arr.min():.1f} - {self.ne_arr.max():.0f} cm^-3")
        print(f"  Hot fraction range: {self.feh_arr.min():.6f} - {self.feh_arr.max():.4f}")
        print(f"  Grid size: {len(self.tec_arr)} x {len(self.teh_arr)} x {len(self.ne_arr)} x {len(self.feh_arr)}")
        print(f"  Species: {', '.join(sorted(self.wavelengths_double.keys()))}")
    
    # ==========================================================================
    # RAY TRACING METHODS
    # ==========================================================================
    
    def trace_ray(self, start_pos, direction, ds=0.1, max_distance=40.0):
        """
        Trace a ray through the plasma model.
        
        Parameters
        ----------
        start_pos : array-like, shape (3,)
            Starting position [x, y, z] in Jupiter radii
        direction : array-like, shape (3,)
            Direction vector (will be normalized)
        ds : float, optional
            Step size along ray in R_J (default: 0.1)
        max_distance : float, optional
            Maximum distance to trace in R_J (default: 40.0)
        
        Returns
        -------
        s_values : ndarray, shape (n_points,)
            Distances along ray in R_J
        positions : ndarray, shape (n_points, 3)
            Cartesian positions [x, y, z] along ray in R_J
        """
        # Normalize direction vector
        direction = np.asarray(direction, dtype=np.float64)
        direction = direction / np.linalg.norm(direction)
        
        # Create ray points with uniform spacing
        n_points = int(max_distance / ds) + 1
        s_values = np.linspace(0.0, max_distance, n_points)
        
        # Calculate positions along ray using broadcasting
        start_pos = np.asarray(start_pos, dtype=np.float64)
        positions = start_pos + s_values[:, np.newaxis] * direction
        
        return s_values, positions
    
    # ==========================================================================
    # VECTORIZED EMISSION INTERPOLATION METHODS
    # ==========================================================================
    
    def interpolate_emission_rates_single_vectorized(self, Te_arr, ne_arr, species_key,
                                                      min_wav=550.0, max_wav=2100.0):
        """
        Interpolate per-ion photon emission rates for all lines of a species
        at multiple LOS points simultaneously using vectorized bilinear interpolation.
        
        Parameters
        ----------
        Te_arr : ndarray, shape (n_los,)
            Electron temperatures along line of sight [eV]
        ne_arr : ndarray, shape (n_los,)
            Electron densities along line of sight [cm^-3]
        species_key : str
            Species identifier ('SP', 'S2P', 'S3P', 'OP', 'O2P', etc.)
        min_wav : float, optional
            Minimum wavelength to include [Angstroms]
        max_wav : float, optional
            Maximum wavelength to include [Angstroms]
        
        Returns
        -------
        wavelengths : ndarray, shape (n_lines,)
            Emission line wavelengths within wavelength range [Angstroms]
        emission_rates : ndarray, shape (n_los, n_lines)
            Per-ion photon emission rates [photons s^-1 ion^-1]
            for each LOS point and each emission line
        
        Notes
        -----
        Uses bilinear interpolation in log10(T)-log10(n) space.
        Points outside table bounds return 0 (no extrapolation).
        """
        if species_key not in self.wavelengths_single:
            return np.array([]), np.zeros((len(Te_arr), 0))
        
        # Get wavelengths and filter to range
        wav_all = self.wavelengths_single[species_key]
        wav_mask = (wav_all >= min_wav) & (wav_all <= max_wav)
        wavelengths = wav_all[wav_mask]
        
        if len(wavelengths) == 0:
            return np.array([]), np.zeros((len(Te_arr), 0))
        
        # Get emission table for this species: shape (n_lines, n_T, n_n)
        emiss_table_all = self.emissivities_single[species_key]
        emiss_table = emiss_table_all[wav_mask]  # (n_lines_filtered, n_T, n_n)
        
        n_los = len(Te_arr)
        n_lines = len(wavelengths)
        
        # Initialize output array
        emission_rates = np.zeros((n_los, n_lines), dtype=np.float64)
        
        # Check for valid points (positive, finite, within bounds)
        valid = (Te_arr > 0) & (ne_arr > 0) & np.isfinite(Te_arr) & np.isfinite(ne_arr)
        in_bounds = (valid & 
                    (Te_arr >= self.temp_arr[0]) & (Te_arr <= self.temp_arr[-1]) &
                    (ne_arr >= self.dens_arr[0]) & (ne_arr <= self.dens_arr[-1]))
        
        if not np.any(in_bounds):
            return wavelengths, emission_rates
        
        # Convert valid points to log space
        log_T = np.log10(Te_arr[in_bounds])
        log_n = np.log10(ne_arr[in_bounds])
        n_valid = np.sum(in_bounds)
        
        # Find bracketing indices in log space using searchsorted
        i_T = np.searchsorted(self.log_temp, log_T).clip(1, len(self.log_temp) - 1)
        i_n = np.searchsorted(self.log_dens, log_n).clip(1, len(self.log_dens) - 1)
        
        i_T0, i_T1 = i_T - 1, i_T
        i_n0, i_n1 = i_n - 1, i_n
        
        # Compute interpolation weights
        log_T0 = self.log_temp[i_T0]
        log_T1 = self.log_temp[i_T1]
        log_n0 = self.log_dens[i_n0]
        log_n1 = self.log_dens[i_n1]
        
        # Safe division for weights
        dT = log_T1 - log_T0
        dn = log_n1 - log_n0
        w_T = np.where(dT != 0, (log_T - log_T0) / dT, 0.0)
        w_n = np.where(dn != 0, (log_n - log_n0) / dn, 0.0)
        
        # Bilinear interpolation weights: (n_valid,)
        w00 = (1.0 - w_T) * (1.0 - w_n)
        w01 = (1.0 - w_T) * w_n
        w10 = w_T * (1.0 - w_n)
        w11 = w_T * w_n
        
        # Extract corner values for all lines and all valid points
        # emiss_table shape: (n_lines, n_T, n_n)
        # Result for each corner: (n_lines, n_valid)
        E_00 = emiss_table[:, i_T0, i_n0]  # (n_lines, n_valid)
        E_01 = emiss_table[:, i_T0, i_n1]
        E_10 = emiss_table[:, i_T1, i_n0]
        E_11 = emiss_table[:, i_T1, i_n1]
        
        # Apply bilinear interpolation: (n_lines, n_valid)
        emiss_valid = (w00 * E_00 + w01 * E_01 + w10 * E_10 + w11 * E_11)
        
        # Transpose to (n_valid, n_lines) and place into output
        emission_rates[in_bounds, :] = emiss_valid.T
        
        return wavelengths, emission_rates
    
    def interpolate_emission_rates_double_vectorized(self, Tec_arr, Teh_arr, ne_arr, feh_arr,
                                                      species_key, min_wav=550.0, max_wav=2100.0):
        """
        Interpolate per-ion photon emission rates for all lines of a species
        at multiple LOS points using vectorized quadrilinear interpolation.
        
        Parameters
        ----------
        Tec_arr : ndarray, shape (n_los,)
            Core electron temperatures along line of sight [eV]
        Teh_arr : ndarray, shape (n_los,)
            Hot electron temperatures along line of sight [eV]
        ne_arr : ndarray, shape (n_los,)
            Total electron densities along line of sight [cm^-3]
        feh_arr : ndarray, shape (n_los,)
            Hot electron fractions along line of sight
        species_key : str
            Species identifier
        min_wav : float, optional
            Minimum wavelength to include [Angstroms]
        max_wav : float, optional
            Maximum wavelength to include [Angstroms]
        
        Returns
        -------
        wavelengths : ndarray, shape (n_lines,)
            Emission line wavelengths within wavelength range [Angstroms]
        emission_rates : ndarray, shape (n_los, n_lines)
            Per-ion photon emission rates [photons s^-1 ion^-1]
            
        Notes
        -----
        For points where feh < feh_min or Teh < Teh_min, falls back to
        single Maxwellian interpolation for numerical stability.
        """
        if not self.double_maxwellian_loaded:
            return self.interpolate_emission_rates_single_vectorized(
                Tec_arr, ne_arr, species_key, min_wav, max_wav)
        
        if species_key not in self.wavelengths_double:
            return np.array([]), np.zeros((len(Tec_arr), 0))
        
        # Get wavelengths and filter to range
        wav_all = self.wavelengths_double[species_key]
        wav_mask = (wav_all >= min_wav) & (wav_all <= max_wav)
        wavelengths = wav_all[wav_mask]
        
        if len(wavelengths) == 0:
            return np.array([]), np.zeros((len(Tec_arr), 0))
        
        # Get emission table: shape (n_lines, n_Tec, n_Teh, n_ne, n_feh)
        emiss_table_all = self.emissivities_double[species_key]
        emiss_table = emiss_table_all[wav_mask]
        
        n_los = len(Tec_arr)
        n_lines = len(wavelengths)
        
        # Initialize output array
        emission_rates = np.zeros((n_los, n_lines), dtype=np.float64)
        
        # Determine which points should use single vs double Maxwellian
        # Use single Maxwellian when feh is too small or Teh is too small
        use_single = (feh_arr < self.feh_arr[0]) | (Teh_arr < self.teh_arr[0])
        use_double = ~use_single
        
        # Handle single Maxwellian fallback points
        if np.any(use_single) and self.single_maxwellian_loaded:
            wav_s, emiss_s = self.interpolate_emission_rates_single_vectorized(
                Tec_arr[use_single], ne_arr[use_single], species_key, min_wav, max_wav)
            if len(wav_s) > 0:
                emission_rates[use_single, :] = emiss_s
        
        # Handle double Maxwellian points
        if not np.any(use_double):
            return wavelengths, emission_rates
        
        # Check for valid double Maxwellian points
        valid = (use_double & 
                (Tec_arr > 0) & (Teh_arr > 0) & (ne_arr > 0) & 
                (feh_arr >= 0) & (feh_arr <= 1) &
                np.isfinite(Tec_arr) & np.isfinite(Teh_arr) & 
                np.isfinite(ne_arr) & np.isfinite(feh_arr))
        
        in_bounds = (valid &
                    (Tec_arr >= self.tec_arr[0]) & (Tec_arr <= self.tec_arr[-1]) &
                    (Teh_arr >= self.teh_arr[0]) & (Teh_arr <= self.teh_arr[-1]) &
                    (ne_arr >= self.ne_arr[0]) & (ne_arr <= self.ne_arr[-1]) &
                    (feh_arr >= self.feh_arr[0]) & (feh_arr <= self.feh_arr[-1]))
        
        if not np.any(in_bounds):
            return wavelengths, emission_rates
        
        # Convert valid points to log space
        log_Tec = np.log10(Tec_arr[in_bounds])
        log_Teh = np.log10(Teh_arr[in_bounds])
        log_ne = np.log10(ne_arr[in_bounds])
        log_feh = np.log10(feh_arr[in_bounds])
        n_valid = np.sum(in_bounds)
        
        # Find bracketing indices
        i_Tec = np.searchsorted(self.log_tec, log_Tec).clip(1, len(self.log_tec) - 1)
        i_Teh = np.searchsorted(self.log_teh, log_Teh).clip(1, len(self.log_teh) - 1)
        i_ne = np.searchsorted(self.log_ne, log_ne).clip(1, len(self.log_ne) - 1)
        i_feh = np.searchsorted(self.log_feh, log_feh).clip(1, len(self.log_feh) - 1)
        
        i0_Tec, i1_Tec = i_Tec - 1, i_Tec
        i0_Teh, i1_Teh = i_Teh - 1, i_Teh
        i0_ne, i1_ne = i_ne - 1, i_ne
        i0_feh, i1_feh = i_feh - 1, i_feh
        
        # Compute interpolation weights
        dTec = self.log_tec[i1_Tec] - self.log_tec[i0_Tec]
        dTeh = self.log_teh[i1_Teh] - self.log_teh[i0_Teh]
        dne = self.log_ne[i1_ne] - self.log_ne[i0_ne]
        dfeh = self.log_feh[i1_feh] - self.log_feh[i0_feh]
        
        w_Tec = np.where(dTec != 0, (log_Tec - self.log_tec[i0_Tec]) / dTec, 0.0)
        w_Teh = np.where(dTeh != 0, (log_Teh - self.log_teh[i0_Teh]) / dTeh, 0.0)
        w_ne = np.where(dne != 0, (log_ne - self.log_ne[i0_ne]) / dne, 0.0)
        w_feh = np.where(dfeh != 0, (log_feh - self.log_feh[i0_feh]) / dfeh, 0.0)
        
        # Quadrilinear interpolation weights for 16 corners
        w0_Tec, w1_Tec = 1.0 - w_Tec, w_Tec
        w0_Teh, w1_Teh = 1.0 - w_Teh, w_Teh
        w0_ne, w1_ne = 1.0 - w_ne, w_ne
        w0_feh, w1_feh = 1.0 - w_feh, w_feh
        
        # Perform quadrilinear interpolation for all lines
        # emiss_table shape: (n_lines, n_Tec, n_Teh, n_ne, n_feh)
        emiss_valid = np.zeros((n_lines, n_valid), dtype=np.float64)
        
        # Loop over 16 corners of 4D hypercube
        for b_Tec in [0, 1]:
            for b_Teh in [0, 1]:
                for b_ne in [0, 1]:
                    for b_feh in [0, 1]:
                        # Get indices for this corner
                        idx_Tec = i0_Tec if b_Tec == 0 else i1_Tec
                        idx_Teh = i0_Teh if b_Teh == 0 else i1_Teh
                        idx_ne = i0_ne if b_ne == 0 else i1_ne
                        idx_feh = i0_feh if b_feh == 0 else i1_feh
                        
                        # Get weights for this corner
                        wt_Tec = w0_Tec if b_Tec == 0 else w1_Tec
                        wt_Teh = w0_Teh if b_Teh == 0 else w1_Teh
                        wt_ne = w0_ne if b_ne == 0 else w1_ne
                        wt_feh = w0_feh if b_feh == 0 else w1_feh
                        
                        # Combined weight
                        weight = wt_Tec * wt_Teh * wt_ne * wt_feh
                        
                        # Extract corner values: (n_lines, n_valid)
                        corner_vals = emiss_table[:, idx_Tec, idx_Teh, idx_ne, idx_feh]
                        
                        # Add weighted contribution
                        emiss_valid += weight * corner_vals
        
        # Transpose and place into output
        emission_rates[in_bounds, :] = emiss_valid.T
        
        return wavelengths, emission_rates
    
    # ==========================================================================
    # LINE-OF-SIGHT INTEGRATION METHODS
    # ==========================================================================
    
    def integrate_species_emission_single(self, s_values, positions, species_key,
                                          Te_los, ne_los, n_ion_los,
                                          min_wav=550.0, max_wav=2100.0):
        """
        Integrate emission for a single species along line of sight (single Maxwellian).
        
        Parameters
        ----------
        s_values : ndarray, shape (n_los,)
            Distances along ray [R_J]
        positions : ndarray, shape (n_los, 3)
            Positions along ray [R_J]
        species_key : str
            Species identifier
        Te_los : ndarray, shape (n_los,)
            Electron temperature along LOS [eV]
        ne_los : ndarray, shape (n_los,)
            Electron density along LOS [cm^-3]
        n_ion_los : ndarray, shape (n_los,)
            Ion density along LOS [cm^-3]
        min_wav, max_wav : float
            Wavelength range [Angstroms]
        
        Returns
        -------
        wavelengths : ndarray
            Emission line wavelengths [Angstroms]
        brightnesses : ndarray
            Integrated brightness for each line [Rayleighs]
        """
        # Get emission rates at all LOS points: (n_los, n_lines)
        wavelengths, emission_rates = self.interpolate_emission_rates_single_vectorized(
            Te_los, ne_los, species_key, min_wav, max_wav)
        
        if len(wavelengths) == 0:
            return np.array([]), np.array([])
        
        # Calculate integrand: emission_rate * n_ion for each point and line
        # emission_rates: (n_los, n_lines), n_ion_los: (n_los,)
        # Result: (n_los, n_lines)
        integrand = emission_rates * n_ion_los[:, np.newaxis]
        
        # Integrate each line along LOS using Simpson's rule
        # Convert from R_J to cm, multiply by 1e-6 for Rayleighs
        n_lines = len(wavelengths)
        brightnesses = np.zeros(n_lines, dtype=np.float64)
        
        if len(s_values) >= 3:
            # Simpson's rule for each line
            for i in range(n_lines):
                brightnesses[i] = RAYLEIGH_FACTOR * simpson(integrand[:, i], x=s_values)
        elif len(s_values) == 2:
            # Trapezoidal rule
            ds = s_values[1] - s_values[0]
            brightnesses = RAYLEIGH_FACTOR * 0.5 * (integrand[0, :] + integrand[1, :]) * ds
        elif len(s_values) == 1:
            ds = 0.1  # Default step size
            brightnesses = RAYLEIGH_FACTOR * integrand[0, :] * ds
        
        return wavelengths, brightnesses
    
    def integrate_species_emission_double(self, s_values, positions, species_key,
                                          Tec_los, Teh_los, ne_los, feh_los, n_ion_los,
                                          min_wav=550.0, max_wav=2100.0):
        """
        Integrate emission for a single species along line of sight (double Maxwellian).
        
        Parameters
        ----------
        s_values : ndarray, shape (n_los,)
            Distances along ray [R_J]
        positions : ndarray, shape (n_los, 3)
            Positions along ray [R_J]
        species_key : str
            Species identifier
        Tec_los : ndarray, shape (n_los,)
            Core electron temperature along LOS [eV]
        Teh_los : ndarray, shape (n_los,)
            Hot electron temperature along LOS [eV]
        ne_los : ndarray, shape (n_los,)
            Total electron density along LOS [cm^-3]
        feh_los : ndarray, shape (n_los,)
            Hot electron fraction along LOS
        n_ion_los : ndarray, shape (n_los,)
            Ion density along LOS [cm^-3]
        min_wav, max_wav : float
            Wavelength range [Angstroms]
        
        Returns
        -------
        wavelengths : ndarray
            Emission line wavelengths [Angstroms]
        brightnesses : ndarray
            Integrated brightness for each line [Rayleighs]
        """
        # Get emission rates at all LOS points: (n_los, n_lines)
        wavelengths, emission_rates = self.interpolate_emission_rates_double_vectorized(
            Tec_los, Teh_los, ne_los, feh_los, species_key, min_wav, max_wav)
        
        if len(wavelengths) == 0:
            return np.array([]), np.array([])
        
        # Calculate integrand
        integrand = emission_rates * n_ion_los[:, np.newaxis]
        
        # Integrate each line along LOS
        n_lines = len(wavelengths)
        brightnesses = np.zeros(n_lines, dtype=np.float64)
        
        if len(s_values) >= 3:
            for i in range(n_lines):
                brightnesses[i] = RAYLEIGH_FACTOR * simpson(integrand[:, i], x=s_values)
        elif len(s_values) == 2:
            ds = s_values[1] - s_values[0]
            brightnesses = RAYLEIGH_FACTOR * 0.5 * (integrand[0, :] + integrand[1, :]) * ds
        elif len(s_values) == 1:
            ds = 0.1
            brightnesses = RAYLEIGH_FACTOR * integrand[0, :] * ds
        
        return wavelengths, brightnesses
    
    # ==========================================================================
    # INSTRUMENT RESPONSE CONVOLUTION
    # ==========================================================================
    
    def convolve_spectrum_erf(self, wavelength_grid, bin_width, 
                               line_wavelengths, line_brightnesses, fwhm=6.0):
        """
        Convolve discrete emission lines with instrument response function using
        analytic ERF-based Gaussian integration.
        
        This method exactly integrates the Gaussian line profile over finite
        wavelength bins, which is more accurate than simple sampling at bin centers.
        
        Parameters
        ----------
        wavelength_grid : ndarray
            Output wavelength bin centers [Angstroms]
        bin_width : float
            Width of each wavelength bin [Angstroms]
        line_wavelengths : ndarray
            Discrete emission line wavelengths [Angstroms]
        line_brightnesses : ndarray
            Discrete emission line brightnesses [Rayleighs]
        fwhm : float, optional
            Full width at half maximum of instrument response [Angstroms]
            Default: 6.0 (appropriate for Europa-UVS, JUICE-UVS)
        
        Returns
        -------
        spectrum : ndarray
            Convolved spectrum [Rayleighs/Angstrom]
            
        Notes
        -----
        The ERF formulation computes the exact integral of a Gaussian profile
        over each bin using:
        
        I_bin = I_line * 0.5 * [erf((λ_hi - λ_0)/(σ√2)) - erf((λ_lo - λ_0)/(σ√2))]
        
        where σ = FWHM / (2√(2 ln 2)) ≈ FWHM / 2.355
        """
        if len(line_wavelengths) == 0 or len(line_brightnesses) == 0:
            return np.zeros_like(wavelength_grid)
        
        # Convert FWHM to Gaussian sigma
        sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        sigma_sqrt2 = sigma * np.sqrt(2.0)
        
        # Calculate bin edges
        bin_lo = wavelength_grid - bin_width / 2.0
        bin_hi = wavelength_grid + bin_width / 2.0
        
        # Reshape for broadcasting: (n_lines, 1) vs (n_bins,)
        line_wav = line_wavelengths[:, np.newaxis]
        line_bright = line_brightnesses[:, np.newaxis]
        
        # Calculate ERF arguments
        erf_arg_lo = (bin_lo - line_wav) / sigma_sqrt2
        erf_arg_hi = (bin_hi - line_wav) / sigma_sqrt2
        
        # Calculate bin contributions from each line
        bin_contrib = 0.5 * (erf(erf_arg_hi) - erf(erf_arg_lo))
        weighted_contrib = line_bright * bin_contrib
        
        # Sum over all lines and convert to R/Å
        spectrum = np.sum(weighted_contrib, axis=0) / bin_width
        
        return spectrum
    
    # ==========================================================================
    # MAIN SPECTRUM CALCULATION METHODS
    # ==========================================================================
    
    def calculate_spectrum_single(self, slit_pos_vec, norm_vec,
                                  wavelength_range=(550, 2100),
                                  bin_width=1.0, fwhm=6.0, ds=0.1):
        """
        Calculate complete LOS-integrated spectrum for single Maxwellian distribution.
        
        This is the main user-facing method for computing UV/optical emission spectra
        from the Io Plasma Torus using single Maxwellian electron distributions.
        
        Parameters
        ----------
        slit_pos_vec : array-like, shape (3,)
            Starting position [x, y, z] in Jupiter radii
        norm_vec : array-like, shape (3,)
            Line-of-sight direction vector (will be normalized)
        wavelength_range : tuple, optional
            (min, max) wavelength in Angstroms (default: 550-2100)
        bin_width : float, optional
            Spectral bin width in Angstroms (default: 1.0)
        fwhm : float, optional
            Instrumental FWHM in Angstroms (default: 6.0)
        ds : float, optional
            Integration step size in R_J (default: 0.1)
        
        Returns
        -------
        wave_bins : ndarray
            Wavelength bin centers [Angstroms]
        spectrum : ndarray
            Convolved spectrum [Rayleighs/Angstrom]
        line_list : list of tuples
            List of (wavelength, brightness, species) for each emission line
            
        Notes
        -----
        The calculation workflow is:
        1. Trace ray through 3D plasma model
        2. Interpolate plasma parameters (Te, ne, n_ion) at each LOS point
        3. For each species: interpolate emission rates and integrate along LOS
        4. Convolve discrete line brightnesses with instrument response
        """
        print(f"Calculating single Maxwellian spectrum:")
        print(f"  LOS start: [{slit_pos_vec[0]:.1f}, {slit_pos_vec[1]:.1f}, {slit_pos_vec[2]:.1f}] R_J")
        print(f"  Direction: [{norm_vec[0]:.2f}, {norm_vec[1]:.2f}, {norm_vec[2]:.2f}]")
        
        min_wav, max_wav = wavelength_range
        
        # Create output wavelength grid
        n_wav = int((max_wav - min_wav) / bin_width) + 1
        wave_bins = np.linspace(min_wav, max_wav, n_wav)
        
        # Trace ray through plasma model
        s_values, positions = self.trace_ray(slit_pos_vec, norm_vec, ds=ds)
        
        # Interpolate plasma parameters along LOS
        ne_los = self.interp_nec(positions)
        Te_los = self.interp_Tec(positions)
        
        # Clean up invalid values
        ne_los = np.where(np.isfinite(ne_los) & (ne_los > 0), ne_los, 0.0)
        Te_los = np.where(np.isfinite(Te_los) & (Te_los > 0), Te_los, 0.0)
        
        # Collect all emission lines from all species
        all_wavelengths = []
        all_brightnesses = []
        all_species = []
        
        # Process each species
        for species_key in self.default_species:
            if species_key not in self._ion_interp_map:
                continue
            
            # Get ion density along LOS
            n_ion_los = self._ion_interp_map[species_key](positions)
            n_ion_los = np.where(np.isfinite(n_ion_los) & (n_ion_los > 0), n_ion_los, 0.0)
            
            # Skip if no valid plasma
            if not np.any((ne_los > 0) & (Te_los > 0) & (n_ion_los > 0)):
                continue
            
            # Integrate emission for this species
            wavelengths, brightnesses = self.integrate_species_emission_single(
                s_values, positions, species_key, Te_los, ne_los, n_ion_los,
                min_wav, max_wav)
            
            if len(wavelengths) > 0:
                all_wavelengths.extend(wavelengths)
                all_brightnesses.extend(brightnesses)
                all_species.extend([species_key] * len(wavelengths))
        
        # Convert to arrays
        all_wavelengths = np.array(all_wavelengths)
        all_brightnesses = np.array(all_brightnesses)
        
        # Create line list for output
        line_list = [(w, b, s) for w, b, s in zip(all_wavelengths, all_brightnesses, all_species)
                     if b > 1e-10]
        
        print(f"  Processed {len(all_wavelengths)} lines, {len(line_list)} with non-zero brightness")
        
        # Convolve with instrument response
        spectrum = self.convolve_spectrum_erf(wave_bins, bin_width, 
                                              all_wavelengths, all_brightnesses, fwhm)
        
        return wave_bins, spectrum, line_list
    
    def calculate_spectrum_double(self, slit_pos_vec, norm_vec,
                                  wavelength_range=(550, 2100),
                                  bin_width=1.0, fwhm=6.0, ds=0.1):
        """
        Calculate complete LOS-integrated spectrum for double Maxwellian distribution.
        
        This method accounts for suprathermal (hot) electron populations from wave
        heating, which significantly enhance high-excitation transitions.
        
        Parameters
        ----------
        slit_pos_vec : array-like, shape (3,)
            Starting position [x, y, z] in Jupiter radii
        norm_vec : array-like, shape (3,)
            Line-of-sight direction vector
        wavelength_range : tuple, optional
            (min, max) wavelength in Angstroms (default: 550-2100)
        bin_width : float, optional
            Spectral bin width in Angstroms (default: 1.0)
        fwhm : float, optional
            Instrumental FWHM in Angstroms (default: 6.0)
        ds : float, optional
            Integration step size in R_J (default: 0.1)
        
        Returns
        -------
        wave_bins : ndarray
            Wavelength bin centers [Angstroms]
        spectrum : ndarray
            Convolved spectrum [Rayleighs/Angstrom]
        line_list : list of tuples
            List of (wavelength, brightness, species) for each emission line
            
        Raises
        ------
        RuntimeError
            If double Maxwellian tables are not loaded
        """
        if not self.double_maxwellian_loaded:
            raise RuntimeError("Double Maxwellian tables not loaded. "
                             "Call load_emission_tables_double() first.")
        
        print(f"Calculating double Maxwellian spectrum:")
        print(f"  LOS start: [{slit_pos_vec[0]:.1f}, {slit_pos_vec[1]:.1f}, {slit_pos_vec[2]:.1f}] R_J")
        print(f"  Direction: [{norm_vec[0]:.2f}, {norm_vec[1]:.2f}, {norm_vec[2]:.2f}]")
        
        min_wav, max_wav = wavelength_range
        
        # Create output wavelength grid
        n_wav = int((max_wav - min_wav) / bin_width) + 1
        wave_bins = np.linspace(min_wav, max_wav, n_wav)
        
        # Trace ray through plasma model
        s_values, positions = self.trace_ray(slit_pos_vec, norm_vec, ds=ds)
        
        # Interpolate plasma parameters along LOS
        ne_los = self.interp_ne_total(positions)
        feh_los = self.interp_feh(positions)
        Tec_los = self.interp_Tec(positions)
        Teh_los = self.interp_Teh(positions)
        
        # Clean up invalid values
        ne_los = np.where(np.isfinite(ne_los) & (ne_los > 0), ne_los, 0.0)
        feh_los = np.where(np.isfinite(feh_los) & (feh_los >= 0), feh_los, 0.0)
        Tec_los = np.where(np.isfinite(Tec_los) & (Tec_los > 0), Tec_los, 0.0)
        Teh_los = np.where(np.isfinite(Teh_los) & (Teh_los > 0), Teh_los, 0.0)
        
        # Collect all emission lines from all species
        all_wavelengths = []
        all_brightnesses = []
        all_species = []
        
        # Process each species
        for species_key in self.default_species:
            if species_key not in self._ion_interp_map:
                continue
            
            # Get ion density along LOS
            n_ion_los = self._ion_interp_map[species_key](positions)
            n_ion_los = np.where(np.isfinite(n_ion_los) & (n_ion_los > 0), n_ion_los, 0.0)
            
            # Skip if no valid plasma
            if not np.any((ne_los > 0) & (Tec_los > 0) & (n_ion_los > 0)):
                continue
            
            # Integrate emission for this species
            wavelengths, brightnesses = self.integrate_species_emission_double(
                s_values, positions, species_key, 
                Tec_los, Teh_los, ne_los, feh_los, n_ion_los,
                min_wav, max_wav)
            
            if len(wavelengths) > 0:
                all_wavelengths.extend(wavelengths)
                all_brightnesses.extend(brightnesses)
                all_species.extend([species_key] * len(wavelengths))
        
        # Convert to arrays
        all_wavelengths = np.array(all_wavelengths)
        all_brightnesses = np.array(all_brightnesses)
        
        # Create line list for output
        line_list = [(w, b, s) for w, b, s in zip(all_wavelengths, all_brightnesses, all_species)
                     if b > 1e-10]
        
        print(f"  Processed {len(all_wavelengths)} lines, {len(line_list)} with non-zero brightness")
        
        # Convolve with instrument response
        spectrum = self.convolve_spectrum_erf(wave_bins, bin_width,
                                              all_wavelengths, all_brightnesses, fwhm)
        
        return wave_bins, spectrum, line_list
    
    # ==========================================================================
    # PLOTTING METHODS
    # ==========================================================================
    
    def plot_spectrum(self, wave_bins, spectrum, line_list=None,
                      title="Jovian Plasma Torus Emission Spectrum", color='C0'):
        """
        Plot the calculated emission spectrum with optional line markers.
        
        Parameters
        ----------
        wave_bins : ndarray
            Wavelength bin centers [Angstroms]
        spectrum : ndarray
            Spectrum [Rayleighs/Angstrom]
        line_list : list, optional
            List of (wavelength, brightness, species) tuples
        title : str, optional
            Plot title
        color : str, optional
            Line color (default: 'C0')
        
        Returns
        -------
        fig : matplotlib.figure.Figure
            The figure object
        """
        fig, ax = plt.subplots(figsize=(12, 6))
        
        ax.plot(wave_bins, spectrum, color=color, linewidth=1.0, label='Spectrum')
        
        # Mark major emission lines if provided
        if line_list and len(line_list) > 0:
            # Sort by brightness and take top lines
            line_list_sorted = sorted(line_list, key=lambda x: x[1], reverse=True)
            
            for i, (wave, brightness, ion) in enumerate(line_list_sorted[:20]):
                ax.axvline(wave, color='gray', alpha=0.2, linestyle='--', linewidth=0.5)
                if i < 5:
                    ion_label = self.ion_names.get(ion, ion)
                    ax.text(wave, ax.get_ylim()[1]*0.95 if ax.get_ylim()[1] > 0 else 1.0,
                           f'{ion_label}\n{wave:.1f}Å',
                           rotation=90, va='top', ha='right', fontsize=7)
        
        ax.set_xlabel('Wavelength [Å]', fontsize=12)
        ax.set_ylabel('Brightness [Rayleighs/Å]', fontsize=12)
        ax.set_title(title, fontsize=14)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(wave_bins[0], wave_bins[-1])
        ax.set_ylim(bottom=0)
        
        # Add total brightness annotation
        total = simpson(spectrum, x=wave_bins)
        ax.text(0.02, 0.98, f'Total: {total:.1f} R',
                transform=ax.transAxes, va='top', ha='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        return fig


# ==============================================================================
# STANDALONE HELPER FUNCTIONS
# ==============================================================================

def simulate_ipt_spectrum_rayleighs_erf_form(wavelength_grid, bin_width,
                                             line_wavelengths, line_brightnesses,
                                             fwhm=6.0):
    """
    Convolve discrete emission lines with instrument response function.
    
    Uses Error Function (ERF) formulation for exact integration of Gaussian
    line profiles over finite wavelength bins. This standalone function matches
    the cm^-3 version for consistency.
    
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
        Full width at half maximum of instrument response [Angstroms]
        Default: 6.0 (appropriate for Europa-UVS, JUICE-UVS)
    
    Returns
    -------
    spectrum : ndarray
        Convolved spectrum [Rayleighs/Angstrom]
    """
    if len(line_wavelengths) == 0 or len(line_brightnesses) == 0:
        return np.zeros_like(wavelength_grid)
    
    # Convert FWHM to Gaussian sigma
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    sigma_sqrt2 = sigma * np.sqrt(2.0)
    
    # Calculate bin edges
    bin_lo = wavelength_grid - bin_width / 2.0
    bin_hi = wavelength_grid + bin_width / 2.0
    
    # Reshape for broadcasting
    line_wav = np.asarray(line_wavelengths)[:, np.newaxis]
    line_bright = np.asarray(line_brightnesses)[:, np.newaxis]
    
    # Calculate ERF arguments
    erf_arg_lo = (bin_lo - line_wav) / sigma_sqrt2
    erf_arg_hi = (bin_hi - line_wav) / sigma_sqrt2
    
    # Calculate bin contributions
    bin_contrib = 0.5 * (erf(erf_arg_hi) - erf(erf_arg_lo))
    weighted_contrib = line_bright * bin_contrib
    
    # Sum and convert to R/Å
    spectrum = np.sum(weighted_contrib, axis=0) / bin_width
    
    return spectrum


# ==============================================================================
# MODULE TEST
# ==============================================================================

if __name__ == "__main__":
    print("="*70)
    print("IPT Line-of-Sight Emission Model - Module Test")
    print("="*70)
    print()
    print("This module provides the JovianUVEmissionRaytracer class for")
    print("calculating UV/optical emission spectra through the Io Plasma Torus.")
    print()
    print("Usage:")
    print("  from IPT_emiss_MOP_community_code import JovianUVEmissionRaytracer")
    print("  raytracer = JovianUVEmissionRaytracer()")
    print("  wave, spectrum, lines = raytracer.calculate_spectrum_single(...)")
    print()
    print("See basic_example_uv_integration_emission_model_use_tables.py for")
    print("complete usage examples.")