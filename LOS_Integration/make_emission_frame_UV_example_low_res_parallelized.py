#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_emission_frame_UV_example_low_res_parallelized.py
======================================================
Generate UV Emission Contour Maps - Memory-Efficient Parallelized Version

This script creates contour plots of UV emission from the Jovian plasma torus
as a function of viewing position, using parallel processing with SHARED MEMORY
to avoid loading duplicate copies of large data files in each worker process.

MEMORY MANAGEMENT:
The plasma model and emission tables total ~7+ GB of data. Without shared memory,
N worker processes would require N × 7 GB = potentially hundreds of GB of RAM.
This implementation uses multiprocessing.shared_memory to load data ONCE and
share it across all worker processes via OS-level memory mapping.

PHYSICS AND METHODOLOGY:
This module performs line-of-sight (LOS) integration through a 3D plasma model
to calculate UV emission brightness. The calculation workflow for each pixel:

1. RAY TRACING: Trace a ray from the observer position through the plasma torus
   with uniform step size ds (in Jupiter radii). The ray path is defined by:
   position(s) = start_position + s * direction_vector

2. 3D INTERPOLATION: At each point along the ray, plasma parameters are
   interpolated from the 3D model using trilinear interpolation:
   - Electron density ne [cm^-3]
   - Electron temperature Te [eV]
   - Ion densities for S+, S++, S+++, O+, O++ [cm^-3]

3. EMISSION RATE LOOKUP: For each LOS position, per-ion photon emission rates
   [photons s^-1 ion^-1] are obtained from pre-computed CHIANTI tables using
   bilinear interpolation in log10(Te)-log10(ne) space.

4. LINE-OF-SIGHT INTEGRATION: Volume emission rates are integrated along the
   ray using Simpson's rule to obtain column brightness:
   
   Brightness [R] = 1e-6 * R_J_cm * integral[ emission_rate * n_ion * ds ]
   
   where the 1e-6 factor converts from [photons s^-1 cm^-2 sr^-1] to Rayleighs.

5. SPECTRAL CONVOLUTION: Discrete emission lines are convolved with a Gaussian
   instrument response function using analytic ERF-based integration over
   wavelength bins.

OUTPUT PRODUCTS:
Three contour plots showing emission as viewed from the -Y direction:
1. Total UV emission (550-2100 Å) - all species integrated
2. S++ emission around 680 Å (S III line from 2s²2p² ³P - 2s2p³ ³S° transition)
3. O+ emission around 833 Å (O II line from 2s²2p³ ²D° - 2s2p⁴ ²P transition)

GRID CONFIGURATION:
- X position: -10 to +10 R_J in 0.1 R_J steps (201 points)
- Z position: -2.5 to +2.5 R_J in 0.1 R_J steps (51 points)
- Total grid: 201 × 51 = 10,251 independent calculations

PARALLELIZATION STRATEGY:
- Main process loads all data into shared memory ONCE
- Worker processes attach to shared memory (no data copying)
- Each worker creates lightweight interpolators from shared arrays
- Grid points processed in optimal chunks for load balancing
- Memory footprint: ~7 GB total (not per worker)

PHYSICAL ASSUMPTIONS:
- Optically thin plasma (no absorption or self-absorption effects)
- Single Maxwellian electron distribution for emission rate calculations
- Local thermodynamic equilibrium for atomic excitation
- Jupiter System III coordinate system (co-rotating with magnetic field)
- Gaussian instrumental response function with specified FWHM

UNITS:
- Positions: Jupiter radii [R_J] where 1 R_J = 71,492 km
- Temperatures: electron volts [eV]
- Densities: particles per cubic centimeter [cm^-3]
- Wavelengths: Angstroms [Å]
- Emission rates: photons per second per ion [photons s^-1 ion^-1]
- Brightnesses: Rayleighs [R], where 1 R = 10^6 photons s^-1 cm^-2 (4π sr)^-1

DIRECTORY STRUCTURE EXPECTED:
- IPT_Emission_MOP_Community_Code/
  - LOS_Integration/           (this script and raytracer module)
  - Emiss_tables/              (CHIANTI emission tables HDF5 files)
  - 3D_Torus_Model/            (3D plasma model HDF5 file)

REFERENCES:
- CHIANTI database: Dere et al. 1997; Del Zanna et al. 2020
- IPT observations: Steffl et al. 2004; Thomas et al. 2004
- Emission modeling: Nerney et al. 2017, 2020, 2022, 2025

AUTHOR: Edward (Eddie) G. Nerney
INSTITUTION: Laboratory for Atmospheric and Space Physics,
             University of Colorado Boulder
DATE: November 2025
LICENSE: MIT - Open source for academic and research use
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import simpson
from scipy.special import erf
import h5py
import sys
from pathlib import Path
import time
import multiprocessing as mp
from multiprocessing import shared_memory
import os
from typing import Tuple, List, Dict, Optional
import atexit

# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================

R_J = 71492.0          # Jupiter radius in km
R_J_CM = 7.1492e9      # Jupiter radius in cm
RAYLEIGH_FACTOR = 1e-6 * R_J_CM  # Conversion factor for Rayleigh integration

# =============================================================================
# SHARED MEMORY MANAGEMENT
# =============================================================================

# Global registry of shared memory blocks for cleanup
_shared_memory_blocks = []


def cleanup_shared_memory():
    """
    Clean up all shared memory blocks on program exit.
    
    This function is registered with atexit to ensure shared memory
    is properly released even if the program exits unexpectedly.
    """
    for shm in _shared_memory_blocks:
        try:
            shm.close()
            shm.unlink()
        except Exception:
            pass


# Register cleanup function
atexit.register(cleanup_shared_memory)


def create_shared_array(arr: np.ndarray, name: str) -> Tuple[shared_memory.SharedMemory, tuple, np.dtype]:
    """
    Create a shared memory block and copy array data into it.
    
    Parameters
    ----------
    arr : np.ndarray
        Source array to share
    name : str
        Unique name for the shared memory block
        
    Returns
    -------
    shm : SharedMemory
        The shared memory block
    shape : tuple
        Shape of the array
    dtype : np.dtype
        Data type of the array
    """
    shm = shared_memory.SharedMemory(create=True, size=arr.nbytes, name=name)
    shared_arr = np.ndarray(arr.shape, dtype=arr.dtype, buffer=shm.buf)
    shared_arr[:] = arr[:]
    _shared_memory_blocks.append(shm)
    return shm, arr.shape, arr.dtype


def get_shared_array(name: str, shape: tuple, dtype: np.dtype) -> np.ndarray:
    """
    Attach to an existing shared memory block and return array view.
    
    Parameters
    ----------
    name : str
        Name of the shared memory block
    shape : tuple
        Shape of the array
    dtype : np.dtype
        Data type of the array
        
    Returns
    -------
    arr : np.ndarray
        Array view into shared memory (read-only recommended)
    """
    shm = shared_memory.SharedMemory(name=name)
    return np.ndarray(shape, dtype=dtype, buffer=shm.buf)


# =============================================================================
# SHARED DATA SPECIFICATIONS (populated by main process)
# =============================================================================

# These will be set by the main process before spawning workers
_shared_data_specs = {}


# =============================================================================
# LIGHTWEIGHT RAYTRACER CLASS (uses shared memory)
# =============================================================================

class SharedMemoryRaytracer:
    """
    Lightweight raytracer that uses shared memory arrays.
    
    This class creates interpolators from pre-loaded shared memory arrays,
    avoiding the need to load large HDF5 files in each worker process.
    Memory footprint per worker is minimal (~100 MB for interpolator metadata).
    
    Attributes
    ----------
    x_axis, y_axis, z_axis : ndarray
        Coordinate axes of the 3D plasma model [R_J]
    interp_* : RegularGridInterpolator
        Interpolators for plasma parameters
    temp_arr, dens_arr : ndarray
        Temperature and density grids for emission tables
    wavelengths_single, emissivities_single : dict
        Emission data organized by species
    """
    
    def __init__(self, shared_specs: Dict):
        """
        Initialize raytracer from shared memory specifications.
        
        Parameters
        ----------
        shared_specs : dict
            Dictionary containing shared memory names, shapes, and dtypes
            for all required arrays.
        """
        # Load coordinate axes from shared memory
        self.x_axis = get_shared_array(
            shared_specs['x_axis']['name'],
            shared_specs['x_axis']['shape'],
            shared_specs['x_axis']['dtype']
        )
        self.y_axis = get_shared_array(
            shared_specs['y_axis']['name'],
            shared_specs['y_axis']['shape'],
            shared_specs['y_axis']['dtype']
        )
        self.z_axis = get_shared_array(
            shared_specs['z_axis']['name'],
            shared_specs['z_axis']['shape'],
            shared_specs['z_axis']['dtype']
        )
        
        # Load plasma parameter arrays from shared memory
        nec = get_shared_array(
            shared_specs['nec']['name'],
            shared_specs['nec']['shape'],
            shared_specs['nec']['dtype']
        )
        nsp = get_shared_array(
            shared_specs['nsp']['name'],
            shared_specs['nsp']['shape'],
            shared_specs['nsp']['dtype']
        )
        ns2p = get_shared_array(
            shared_specs['ns2p']['name'],
            shared_specs['ns2p']['shape'],
            shared_specs['ns2p']['dtype']
        )
        ns3p = get_shared_array(
            shared_specs['ns3p']['name'],
            shared_specs['ns3p']['shape'],
            shared_specs['ns3p']['dtype']
        )
        nop = get_shared_array(
            shared_specs['nop']['name'],
            shared_specs['nop']['shape'],
            shared_specs['nop']['dtype']
        )
        no2p = get_shared_array(
            shared_specs['no2p']['name'],
            shared_specs['no2p']['shape'],
            shared_specs['no2p']['dtype']
        )
        Tec = get_shared_array(
            shared_specs['Tec']['name'],
            shared_specs['Tec']['shape'],
            shared_specs['Tec']['dtype']
        )
        
        # Create 3D interpolators for plasma parameters
        interp_kwargs = dict(method='linear', bounds_error=False, fill_value=0.0)
        grid_points = (self.x_axis, self.y_axis, self.z_axis)
        
        self.interp_nec = RegularGridInterpolator(grid_points, nec, **interp_kwargs)
        self.interp_Tec = RegularGridInterpolator(grid_points, Tec, **interp_kwargs)
        self.interp_nsp = RegularGridInterpolator(grid_points, nsp, **interp_kwargs)
        self.interp_ns2p = RegularGridInterpolator(grid_points, ns2p, **interp_kwargs)
        self.interp_ns3p = RegularGridInterpolator(grid_points, ns3p, **interp_kwargs)
        self.interp_nop = RegularGridInterpolator(grid_points, nop, **interp_kwargs)
        self.interp_no2p = RegularGridInterpolator(grid_points, no2p, **interp_kwargs)
        
        # Ion interpolator mapping
        self._ion_interp_map = {
            'SP': self.interp_nsp,
            'S2P': self.interp_ns2p,
            'S3P': self.interp_ns3p,
            'OP': self.interp_nop,
            'O2P': self.interp_no2p,
        }
        
        # Load emission table parameters from shared memory
        self.temp_arr = get_shared_array(
            shared_specs['temp_arr']['name'],
            shared_specs['temp_arr']['shape'],
            shared_specs['temp_arr']['dtype']
        )
        self.dens_arr = get_shared_array(
            shared_specs['dens_arr']['name'],
            shared_specs['dens_arr']['shape'],
            shared_specs['dens_arr']['dtype']
        )
        
        # Create log-space grids for interpolation
        self.log_temp = np.log10(self.temp_arr)
        self.log_dens = np.log10(self.dens_arr)
        
        # Load emission data by species
        self.wavelengths_single = {}
        self.emissivities_single = {}
        self.default_species = ['SP', 'S2P', 'S3P', 'OP', 'O2P']
        
        for species in self.default_species:
            wav_key = f'wavelengths_{species}'
            emiss_key = f'emissivities_{species}'
            
            if wav_key in shared_specs and emiss_key in shared_specs:
                self.wavelengths_single[species] = get_shared_array(
                    shared_specs[wav_key]['name'],
                    shared_specs[wav_key]['shape'],
                    shared_specs[wav_key]['dtype']
                )
                self.emissivities_single[species] = get_shared_array(
                    shared_specs[emiss_key]['name'],
                    shared_specs[emiss_key]['shape'],
                    shared_specs[emiss_key]['dtype']
                )
    
    def trace_ray(self, start_pos: np.ndarray, direction: np.ndarray,
                  ds: float = 0.1, max_distance: float = 40.0) -> Tuple[np.ndarray, np.ndarray]:
        """
        Trace a ray through the plasma model.
        
        Parameters
        ----------
        start_pos : ndarray, shape (3,)
            Starting position [x, y, z] in Jupiter radii
        direction : ndarray, shape (3,)
            Direction vector (will be normalized)
        ds : float
            Step size along ray [R_J]
        max_distance : float
            Maximum distance to trace [R_J]
            
        Returns
        -------
        s_values : ndarray
            Distances along ray [R_J]
        positions : ndarray, shape (n_points, 3)
            Positions along ray [R_J]
        """
        direction = np.asarray(direction, dtype=np.float64)
        direction = direction / np.linalg.norm(direction)
        
        n_points = int(max_distance / ds) + 1
        s_values = np.linspace(0.0, max_distance, n_points)
        
        start_pos = np.asarray(start_pos, dtype=np.float64)
        positions = start_pos + s_values[:, np.newaxis] * direction
        
        return s_values, positions
    
    def interpolate_emission_rates_vectorized(self, Te_arr: np.ndarray, ne_arr: np.ndarray,
                                               species_key: str, min_wav: float = 550.0,
                                               max_wav: float = 2100.0) -> Tuple[np.ndarray, np.ndarray]:
        """
        Interpolate per-ion emission rates for all lines at multiple LOS points.
        
        Uses vectorized bilinear interpolation in log10(T)-log10(n) space.
        
        Parameters
        ----------
        Te_arr : ndarray, shape (n_los,)
            Electron temperatures along LOS [eV]
        ne_arr : ndarray, shape (n_los,)
            Electron densities along LOS [cm^-3]
        species_key : str
            Species identifier ('SP', 'S2P', etc.)
        min_wav, max_wav : float
            Wavelength range to include [Angstroms]
            
        Returns
        -------
        wavelengths : ndarray
            Emission line wavelengths [Angstroms]
        emission_rates : ndarray, shape (n_los, n_lines)
            Per-ion emission rates [photons s^-1 ion^-1]
        """
        if species_key not in self.wavelengths_single:
            return np.array([]), np.zeros((len(Te_arr), 0))
        
        # Get wavelengths and filter to range
        wav_all = self.wavelengths_single[species_key]
        wav_mask = (wav_all >= min_wav) & (wav_all <= max_wav)
        wavelengths = wav_all[wav_mask]
        
        if len(wavelengths) == 0:
            return np.array([]), np.zeros((len(Te_arr), 0))
        
        # Get emission table: shape (n_lines, n_T, n_n)
        emiss_table = self.emissivities_single[species_key][wav_mask]
        
        n_los = len(Te_arr)
        n_lines = len(wavelengths)
        emission_rates = np.zeros((n_los, n_lines), dtype=np.float64)
        
        # Find valid points within table bounds
        valid = ((Te_arr > 0) & (ne_arr > 0) & 
                 np.isfinite(Te_arr) & np.isfinite(ne_arr))
        in_bounds = (valid & 
                    (Te_arr >= self.temp_arr[0]) & (Te_arr <= self.temp_arr[-1]) &
                    (ne_arr >= self.dens_arr[0]) & (ne_arr <= self.dens_arr[-1]))
        
        if not np.any(in_bounds):
            return wavelengths, emission_rates
        
        # Convert to log space for interpolation
        log_T = np.log10(Te_arr[in_bounds])
        log_n = np.log10(ne_arr[in_bounds])
        
        # Find bracketing indices
        i_T = np.searchsorted(self.log_temp, log_T).clip(1, len(self.log_temp) - 1)
        i_n = np.searchsorted(self.log_dens, log_n).clip(1, len(self.log_dens) - 1)
        
        i_T0, i_T1 = i_T - 1, i_T
        i_n0, i_n1 = i_n - 1, i_n
        
        # Compute interpolation weights
        log_T0 = self.log_temp[i_T0]
        log_T1 = self.log_temp[i_T1]
        log_n0 = self.log_dens[i_n0]
        log_n1 = self.log_dens[i_n1]
        
        dT = log_T1 - log_T0
        dn = log_n1 - log_n0
        w_T = np.where(dT != 0, (log_T - log_T0) / dT, 0.0)
        w_n = np.where(dn != 0, (log_n - log_n0) / dn, 0.0)
        
        # Bilinear interpolation weights
        w00 = (1.0 - w_T) * (1.0 - w_n)
        w01 = (1.0 - w_T) * w_n
        w10 = w_T * (1.0 - w_n)
        w11 = w_T * w_n
        
        # Extract corner values and interpolate
        E_00 = emiss_table[:, i_T0, i_n0]
        E_01 = emiss_table[:, i_T0, i_n1]
        E_10 = emiss_table[:, i_T1, i_n0]
        E_11 = emiss_table[:, i_T1, i_n1]
        
        emiss_valid = w00 * E_00 + w01 * E_01 + w10 * E_10 + w11 * E_11
        emission_rates[in_bounds, :] = emiss_valid.T
        
        return wavelengths, emission_rates
    
    def integrate_species_emission(self, s_values: np.ndarray, species_key: str,
                                   Te_los: np.ndarray, ne_los: np.ndarray,
                                   n_ion_los: np.ndarray, min_wav: float = 550.0,
                                   max_wav: float = 2100.0) -> Tuple[np.ndarray, np.ndarray]:
        """
        Integrate emission for a species along line of sight.
        
        Parameters
        ----------
        s_values : ndarray
            Distances along ray [R_J]
        species_key : str
            Species identifier
        Te_los : ndarray
            Electron temperature along LOS [eV]
        ne_los : ndarray
            Electron density along LOS [cm^-3]
        n_ion_los : ndarray
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
        wavelengths, emission_rates = self.interpolate_emission_rates_vectorized(
            Te_los, ne_los, species_key, min_wav, max_wav)
        
        if len(wavelengths) == 0:
            return np.array([]), np.array([])
        
        # Calculate integrand: emission_rate × ion_density
        integrand = emission_rates * n_ion_los[:, np.newaxis]
        
        # Integrate along LOS using Simpson's rule
        n_lines = len(wavelengths)
        brightnesses = np.zeros(n_lines, dtype=np.float64)
        
        if len(s_values) >= 3:
            for i in range(n_lines):
                brightnesses[i] = RAYLEIGH_FACTOR * simpson(integrand[:, i], x=s_values)
        elif len(s_values) == 2:
            ds = s_values[1] - s_values[0]
            brightnesses = RAYLEIGH_FACTOR * 0.5 * (integrand[0, :] + integrand[1, :]) * ds
        
        return wavelengths, brightnesses
    
    def convolve_spectrum_erf(self, wavelength_grid: np.ndarray, bin_width: float,
                              line_wavelengths: np.ndarray, line_brightnesses: np.ndarray,
                              fwhm: float = 6.0) -> np.ndarray:
        """
        Convolve emission lines with Gaussian instrument response using ERF method.
        
        Parameters
        ----------
        wavelength_grid : ndarray
            Output wavelength bin centers [Angstroms]
        bin_width : float
            Wavelength bin width [Angstroms]
        line_wavelengths : ndarray
            Discrete line wavelengths [Angstroms]
        line_brightnesses : ndarray
            Discrete line brightnesses [Rayleighs]
        fwhm : float
            Instrumental FWHM [Angstroms]
            
        Returns
        -------
        spectrum : ndarray
            Convolved spectrum [Rayleighs/Angstrom]
        """
        if len(line_wavelengths) == 0 or len(line_brightnesses) == 0:
            return np.zeros_like(wavelength_grid)
        
        sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        sigma_sqrt2 = sigma * np.sqrt(2.0)
        
        bin_lo = wavelength_grid - bin_width / 2.0
        bin_hi = wavelength_grid + bin_width / 2.0
        
        line_wav = line_wavelengths[:, np.newaxis]
        line_bright = line_brightnesses[:, np.newaxis]
        
        erf_arg_lo = (bin_lo - line_wav) / sigma_sqrt2
        erf_arg_hi = (bin_hi - line_wav) / sigma_sqrt2
        
        bin_contrib = 0.5 * (erf(erf_arg_hi) - erf(erf_arg_lo))
        weighted_contrib = line_bright * bin_contrib
        
        spectrum = np.sum(weighted_contrib, axis=0) / bin_width
        
        return spectrum
    
    def calculate_spectrum(self, slit_pos_vec: np.ndarray, norm_vec: np.ndarray,
                           wavelength_range: Tuple[float, float] = (550, 2100),
                           bin_width: float = 1.0, fwhm: float = 6.0,
                           ds: float = 0.1) -> Tuple[np.ndarray, np.ndarray, List]:
        """
        Calculate complete LOS-integrated spectrum.
        
        Parameters
        ----------
        slit_pos_vec : ndarray, shape (3,)
            Starting position [x, y, z] in R_J
        norm_vec : ndarray, shape (3,)
            LOS direction vector
        wavelength_range : tuple
            (min, max) wavelength [Angstroms]
        bin_width : float
            Spectral bin width [Angstroms]
        fwhm : float
            Instrumental FWHM [Angstroms]
        ds : float
            Integration step size [R_J]
            
        Returns
        -------
        wave_bins : ndarray
            Wavelength bin centers [Angstroms]
        spectrum : ndarray
            Convolved spectrum [Rayleighs/Angstrom]
        line_list : list
            List of (wavelength, brightness, species) tuples
        """
        min_wav, max_wav = wavelength_range
        
        # Create output wavelength grid
        n_wav = int((max_wav - min_wav) / bin_width) + 1
        wave_bins = np.linspace(min_wav, max_wav, n_wav)
        
        # Trace ray through plasma
        s_values, positions = self.trace_ray(slit_pos_vec, norm_vec, ds=ds)
        
        # Interpolate plasma parameters along LOS
        ne_los = self.interp_nec(positions)
        Te_los = self.interp_Tec(positions)
        
        # Clean invalid values
        ne_los = np.where(np.isfinite(ne_los) & (ne_los > 0), ne_los, 0.0)
        Te_los = np.where(np.isfinite(Te_los) & (Te_los > 0), Te_los, 0.0)
        
        # Collect emissions from all species
        all_wavelengths = []
        all_brightnesses = []
        all_species = []
        
        for species_key in self.default_species:
            if species_key not in self._ion_interp_map:
                continue
            
            # Get ion density along LOS
            n_ion_los = self._ion_interp_map[species_key](positions)
            n_ion_los = np.where(np.isfinite(n_ion_los) & (n_ion_los > 0), n_ion_los, 0.0)
            
            if not np.any((ne_los > 0) & (Te_los > 0) & (n_ion_los > 0)):
                continue
            
            wavelengths, brightnesses = self.integrate_species_emission(
                s_values, species_key, Te_los, ne_los, n_ion_los, min_wav, max_wav)
            
            if len(wavelengths) > 0:
                all_wavelengths.extend(wavelengths)
                all_brightnesses.extend(brightnesses)
                all_species.extend([species_key] * len(wavelengths))
        
        all_wavelengths = np.array(all_wavelengths)
        all_brightnesses = np.array(all_brightnesses)
        
        line_list = [(w, b, s) for w, b, s in zip(all_wavelengths, all_brightnesses, all_species)
                     if b > 1e-10]
        
        spectrum = self.convolve_spectrum_erf(wave_bins, bin_width,
                                              all_wavelengths, all_brightnesses, fwhm)
        
        return wave_bins, spectrum, line_list


# =============================================================================
# DATA LOADING FUNCTIONS
# =============================================================================

def get_default_paths() -> Tuple[Path, Path]:
    """
    Get default paths for plasma model and emission tables.
    
    Returns
    -------
    plasma_file : Path
        Path to plasma model HDF5 file
    emission_file : Path
        Path to single Maxwellian emission tables
    """
    module_dir = Path(__file__).parent
    plasma_dir = module_dir.parent / "3D_Torus_Model"
    emiss_dir = module_dir.parent / "Emiss_tables"
    
    plasma_file = plasma_dir / "jovian_plasma_interpolated_381x381x231.h5"
    emission_file = emiss_dir / "CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5"
    
    return plasma_file, emission_file


def load_and_share_data() -> Dict:
    """
    Load all required data and create shared memory blocks.
    
    This function is called ONCE in the main process before spawning workers.
    It loads the plasma model and emission tables, then creates shared memory
    blocks so all worker processes can access the same data without copying.
    
    Returns
    -------
    shared_specs : dict
        Dictionary containing shared memory names, shapes, and dtypes
        for all arrays that workers need to access.
    """
    plasma_file, emission_file = get_default_paths()
    
    if not plasma_file.exists():
        raise FileNotFoundError(f"Plasma model not found: {plasma_file}")
    if not emission_file.exists():
        raise FileNotFoundError(f"Emission tables not found: {emission_file}")
    
    shared_specs = {}
    
    # =========================================================================
    # LOAD PLASMA MODEL
    # =========================================================================
    print("Loading plasma model into shared memory...")
    
    with h5py.File(plasma_file, 'r') as f:
        x_axis = f['coordinates/x'][:].astype(np.float64)
        y_axis = f['coordinates/y'][:].astype(np.float64)
        z_axis = f['coordinates/z'][:].astype(np.float64)
        
        nec = f['data/ne_c'][:].astype(np.float64)
        nsp = f['data/nsp'][:].astype(np.float64)
        ns2p = f['data/ns2p'][:].astype(np.float64)
        ns3p = f['data/ns3p'][:].astype(np.float64)
        nop = f['data/nop'][:].astype(np.float64)
        no2p = f['data/no2p'][:].astype(np.float64)
        Tec = f['data/Te_c'][:].astype(np.float64)
    
    # Clean invalid values
    for arr in [nec, nsp, ns2p, ns3p, nop, no2p, Tec]:
        mask = ~np.isfinite(arr) | (arr < 0)
        arr[mask] = 0.0
    
    # Create shared memory for plasma model arrays
    shm, shape, dtype = create_shared_array(x_axis, 'x_axis')
    shared_specs['x_axis'] = {'name': 'x_axis', 'shape': shape, 'dtype': dtype}
    
    shm, shape, dtype = create_shared_array(y_axis, 'y_axis')
    shared_specs['y_axis'] = {'name': 'y_axis', 'shape': shape, 'dtype': dtype}
    
    shm, shape, dtype = create_shared_array(z_axis, 'z_axis')
    shared_specs['z_axis'] = {'name': 'z_axis', 'shape': shape, 'dtype': dtype}
    
    shm, shape, dtype = create_shared_array(nec, 'nec')
    shared_specs['nec'] = {'name': 'nec', 'shape': shape, 'dtype': dtype}
    
    shm, shape, dtype = create_shared_array(nsp, 'nsp')
    shared_specs['nsp'] = {'name': 'nsp', 'shape': shape, 'dtype': dtype}
    
    shm, shape, dtype = create_shared_array(ns2p, 'ns2p')
    shared_specs['ns2p'] = {'name': 'ns2p', 'shape': shape, 'dtype': dtype}
    
    shm, shape, dtype = create_shared_array(ns3p, 'ns3p')
    shared_specs['ns3p'] = {'name': 'ns3p', 'shape': shape, 'dtype': dtype}
    
    shm, shape, dtype = create_shared_array(nop, 'nop')
    shared_specs['nop'] = {'name': 'nop', 'shape': shape, 'dtype': dtype}
    
    shm, shape, dtype = create_shared_array(no2p, 'no2p')
    shared_specs['no2p'] = {'name': 'no2p', 'shape': shape, 'dtype': dtype}
    
    shm, shape, dtype = create_shared_array(Tec, 'Tec')
    shared_specs['Tec'] = {'name': 'Tec', 'shape': shape, 'dtype': dtype}
    
    print(f"  Plasma model loaded: grid shape {nec.shape}")
    print(f"  X range: {x_axis.min():.1f} to {x_axis.max():.1f} R_J")
    print(f"  Y range: {y_axis.min():.1f} to {y_axis.max():.1f} R_J")
    print(f"  Z range: {z_axis.min():.1f} to {z_axis.max():.1f} R_J")
    
    # =========================================================================
    # LOAD EMISSION TABLES (Single Maxwellian only - saves ~4 GB)
    # =========================================================================
    print("\nLoading emission tables into shared memory...")
    
    with h5py.File(emission_file, 'r') as f:
        temp_arr = f['T'][:].astype(np.float64)
        dens_arr = f['n'][:].astype(np.float64)
        emissivity_all = f['emiss'][:]
        wavelength_all = f['wavelength'][:]
        species_all = f['species'][:]
        
        if species_all.dtype.kind == 'S' or species_all.dtype.kind == 'O':
            species_all = np.array([s.decode('utf-8') if isinstance(s, bytes) else str(s) 
                                   for s in species_all])
    
    shm, shape, dtype = create_shared_array(temp_arr, 'temp_arr')
    shared_specs['temp_arr'] = {'name': 'temp_arr', 'shape': shape, 'dtype': dtype}
    
    shm, shape, dtype = create_shared_array(dens_arr, 'dens_arr')
    shared_specs['dens_arr'] = {'name': 'dens_arr', 'shape': shape, 'dtype': dtype}
    
    print(f"  Temperature range: {temp_arr.min():.3f} - {temp_arr.max():.1f} eV")
    print(f"  Density range: {dens_arr.min():.1f} - {dens_arr.max():.0f} cm^-3")
    
    # Organize emission data by species
    unique_species = []
    for s in species_all:
        if s not in unique_species:
            unique_species.append(s)
    
    for species_key in unique_species:
        indices = np.where(species_all == species_key)[0]
        wavelengths = wavelength_all[indices].astype(np.float64)
        
        # Transpose to (n_lines, n_T, n_n) for vectorized interpolation
        emiss_species = emissivity_all[:, :, indices]
        emissivities = np.ascontiguousarray(
            np.transpose(emiss_species, (2, 0, 1)).astype(np.float64)
        )
        
        wav_name = f'wavelengths_{species_key}'
        emiss_name = f'emissivities_{species_key}'
        
        shm, shape, dtype = create_shared_array(wavelengths, wav_name)
        shared_specs[wav_name] = {'name': wav_name, 'shape': shape, 'dtype': dtype}
        
        shm, shape, dtype = create_shared_array(emissivities, emiss_name)
        shared_specs[emiss_name] = {'name': emiss_name, 'shape': shape, 'dtype': dtype}
    
    print(f"  Species loaded: {', '.join(unique_species)}")
    
    # Calculate total shared memory size
    total_bytes = sum(shm.size for shm in _shared_memory_blocks)
    print(f"\nTotal shared memory allocated: {total_bytes / 1e9:.2f} GB")
    
    return shared_specs


# =============================================================================
# WORKER FUNCTIONS
# =============================================================================

# Global raytracer for worker processes
_raytracer = None


def init_worker(shared_specs: Dict):
    """
    Initialize worker process with shared memory raytracer.
    
    Parameters
    ----------
    shared_specs : dict
        Shared memory specifications from main process
    """
    global _raytracer
    _raytracer = SharedMemoryRaytracer(shared_specs)


def integrate_spectral_line(slit_pos_vec: np.ndarray, norm_vec: np.ndarray,
                            center_wavelength: float, fwhm: float = 6.0,
                            ds: float = 0.1) -> float:
    """
    Integrate emission around a specific spectral line.
    
    Parameters
    ----------
    slit_pos_vec : ndarray
        LOS starting position [R_J]
    norm_vec : ndarray
        LOS direction vector
    center_wavelength : float
        Line center wavelength [Angstroms]
    fwhm : float
        Instrumental FWHM [Angstroms]
    ds : float
        Ray tracing step size [R_J]
        
    Returns
    -------
    integrated_emission : float
        Integrated brightness [Rayleighs]
    """
    global _raytracer
    
    sigma = fwhm / 2.35482
    wave_min = center_wavelength - 3.0 * sigma
    wave_max = center_wavelength + 3.0 * sigma
    
    wave_bins, spectrum, _ = _raytracer.calculate_spectrum(
        slit_pos_vec, norm_vec,
        wavelength_range=(wave_min - 1.0, wave_max + 1.0),
        bin_width=0.1,
        fwhm=fwhm,
        ds=ds
    )
    
    mask = (wave_bins >= wave_min) & (wave_bins <= wave_max)
    
    if np.any(mask):
        return simpson(spectrum[mask], x=wave_bins[mask])
    return 0.0


def calculate_emissions_for_point(args: tuple) -> Tuple[float, float, float]:
    """
    Calculate all emission types for a single grid point.
    
    Parameters
    ----------
    args : tuple
        (x0, z0, y_start, norm_vec, s3_wavelength, o2_wavelength, fwhm, ds)
        
    Returns
    -------
    tuple : (total_emission, s3_emission, o2_emission)
        Emissions in Rayleighs
    """
    global _raytracer
    
    x0, z0, y_start, norm_vec, s3_wavelength, o2_wavelength, fwhm, ds = args
    slit_pos_vec = np.array([x0, y_start, z0], dtype=np.float64)
    
    # Total UV emission (550-2100 Å)
    wave_bins, spectrum, _ = _raytracer.calculate_spectrum(
        slit_pos_vec, norm_vec,
        wavelength_range=(550, 2100),
        bin_width=1.0,
        fwhm=fwhm,
        ds=ds
    )
    total_emission = simpson(spectrum, x=wave_bins)
    
    # S++ emission (680 Å)
    s3_emission = integrate_spectral_line(
        slit_pos_vec, norm_vec, s3_wavelength, fwhm=fwhm, ds=ds)
    
    # O+ emission (833 Å)
    o2_emission = integrate_spectral_line(
        slit_pos_vec, norm_vec, o2_wavelength, fwhm=fwhm, ds=ds)
    
    return total_emission, s3_emission, o2_emission


def process_chunk(chunk_data: List[tuple]) -> List[tuple]:
    """
    Process a chunk of grid points.
    
    Parameters
    ----------
    chunk_data : list
        List of (i, j, args) tuples
        
    Returns
    -------
    results : list
        List of (i, j, total, s3, o2) tuples
    """
    results = []
    for i, j, args in chunk_data:
        total, s3, o2 = calculate_emissions_for_point(args)
        results.append((i, j, total, s3, o2))
    return results


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    """
    Main function to generate UV emission contour maps.
    
    Uses shared memory to efficiently parallelize calculations across
    all CPU cores without duplicating large data structures.
    """
    print("=" * 70)
    print("UV Emission Contour Maps - Memory-Efficient Parallelized Version")
    print("=" * 70)
    print()
    
    # =========================================================================
    # SYSTEM CONFIGURATION
    # =========================================================================
    
    num_cores = mp.cpu_count()
    num_workers = max(1, num_cores - 1) if os.name == 'nt' else num_cores
    
    print("System Configuration:")
    print(f"  Total CPU cores detected: {num_cores}")
    print(f"  Worker processes to use: {num_workers}")
    print(f"  Platform: {os.name} ({'Windows' if os.name == 'nt' else 'Unix-like'})")
    print()
    
    # =========================================================================
    # LOAD DATA INTO SHARED MEMORY
    # =========================================================================
    
    print("Loading data into shared memory (one-time operation)...")
    try:
        shared_specs = load_and_share_data()
    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("\nPlease ensure data files are in the correct directories.")
        return
    
    # =========================================================================
    # GRID CONFIGURATION
    # =========================================================================
    
    print("\nSetting up calculation grid...")
    
    x_grid = np.arange(-10.0, 10.1, 0.1)
    z_grid = np.arange(-2.5, 2.6, 0.1)
    
    nx = len(x_grid)
    nz = len(z_grid)
    total_calculations = nx * nz
    
    print(f"Grid dimensions: {nx} × {nz} = {total_calculations} total calculations")
    print(f"X range: {x_grid[0]:.1f} to {x_grid[-1]:.1f} R_J ({nx} points)")
    print(f"Z range: {z_grid[0]:.1f} to {z_grid[-1]:.1f} R_J ({nz} points)")
    
    # Output arrays
    total_emission = np.zeros((nz, nx), dtype=np.float64)
    s3_emission = np.zeros((nz, nx), dtype=np.float64)
    o2_emission = np.zeros((nz, nx), dtype=np.float64)
    
    # =========================================================================
    # OBSERVATION PARAMETERS
    # =========================================================================
    
    y_start = -20.0
    norm_vec = np.array([0.0, 1.0, 0.0], dtype=np.float64)
    s3_wavelength = 680.0
    o2_wavelength = 833.0
    fwhm = 6.0
    ds = 0.1
    
    print(f"\nObservation Parameters:")
    print(f"  Starting Y position: {y_start:.1f} R_J")
    print(f"  Viewing direction: +Y")
    print(f"  S++ line center: {s3_wavelength:.1f} Å")
    print(f"  O+ line center: {o2_wavelength:.1f} Å")
    print(f"  Instrumental FWHM: {fwhm:.2f} Å")
    print(f"  Ray tracing step: {ds:.2f} R_J")
    
    # =========================================================================
    # PREPARE CALCULATION ARGUMENTS
    # =========================================================================
    
    print("\nPreparing calculation tasks...")
    
    all_args = []
    for i, x0 in enumerate(x_grid):
        for j, z0 in enumerate(z_grid):
            args = (x0, z0, y_start, norm_vec, s3_wavelength, 
                    o2_wavelength, fwhm, ds)
            all_args.append((i, j, args))
    
    chunk_size = max(1, total_calculations // (num_workers * 20))
    chunk_size = min(chunk_size, 50)
    
    print(f"\nParallelization Strategy:")
    print(f"  Chunk size: {chunk_size} calculations per task")
    print(f"  Total chunks: {(total_calculations + chunk_size - 1) // chunk_size}")
    print(f"  Memory mode: SHARED (data loaded once, ~7 GB total)")
    
    chunks = []
    for i in range(0, len(all_args), chunk_size):
        chunks.append(all_args[i:i + chunk_size])
    
    # =========================================================================
    # PARALLEL COMPUTATION WITH SHARED MEMORY
    # =========================================================================
    
    print("\nCalculating emissions in parallel...")
    print("Progress: ", end="", flush=True)
    
    start_time = time.time()
    calculations_completed = 0
    last_percent = 0
    
    # Use spawn context for Windows compatibility
    ctx = mp.get_context('spawn')
    
    with ctx.Pool(processes=num_workers, 
                  initializer=init_worker, 
                  initargs=(shared_specs,)) as pool:
        
        for chunk_results in pool.imap_unordered(process_chunk, chunks):
            for i, j, total, s3, o2 in chunk_results:
                total_emission[j, i] = total
                s3_emission[j, i] = s3
                o2_emission[j, i] = o2
            
            calculations_completed += len(chunk_results)
            percent_complete = int((calculations_completed / total_calculations) * 100)
            
            if percent_complete >= last_percent + 5:
                print(f"{percent_complete}%...", end="", flush=True)
                last_percent = percent_complete
    
    elapsed_time = time.time() - start_time
    print(" Done!")
    
    # =========================================================================
    # PERFORMANCE STATISTICS
    # =========================================================================
    
    print(f"\nPerformance Statistics:")
    print(f"  Total calculation time: {elapsed_time:.1f} seconds")
    print(f"  Average time per point: {elapsed_time/total_calculations:.3f} seconds")
    print(f"  Throughput: {total_calculations/elapsed_time:.1f} calculations/second")
    
    # =========================================================================
    # EMISSION STATISTICS
    # =========================================================================
    
    print("\n" + "=" * 70)
    print("Emission Statistics")
    print("-" * 70)
    
    print(f"Total UV emission (550-2100 Å):")
    print(f"  Min: {np.min(total_emission):.1f} R")
    print(f"  Max: {np.max(total_emission):.1f} R")
    print(f"  Mean: {np.mean(total_emission):.1f} R")
    print(f"  Median: {np.median(total_emission):.1f} R")
    
    print(f"\nS++ emission (680 Å ± 3σ):")
    print(f"  Min: {np.min(s3_emission):.2f} R")
    print(f"  Max: {np.max(s3_emission):.2f} R")
    print(f"  Mean: {np.mean(s3_emission):.2f} R")
    print(f"  Median: {np.median(s3_emission):.2f} R")
    
    print(f"\nO+ emission (833 Å ± 3σ):")
    print(f"  Min: {np.min(o2_emission):.2f} R")
    print(f"  Max: {np.max(o2_emission):.2f} R")
    print(f"  Mean: {np.mean(o2_emission):.2f} R")
    print(f"  Median: {np.median(o2_emission):.2f} R")
    
    # =========================================================================
    # CREATE PLOTS
    # =========================================================================
    
    print("\n" + "=" * 70)
    print("Creating Contour Plots")
    print("-" * 70)
    
    X, Z = np.meshgrid(x_grid, z_grid)
    
    # Combined figure
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    ax1 = axes[0]
    cs1 = ax1.contourf(X, Z, total_emission, levels=100, cmap='viridis')
    ax1.contour(X, Z, total_emission, levels=10, colors='black', alpha=0.3, linewidths=0.5)
    ax1.set_xlabel('X Position (R_J)', fontsize=12)
    ax1.set_ylabel('Z Position (R_J)', fontsize=12)
    ax1.set_title('Total UV Emission (550-2100 Å)', fontsize=14)
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)
    plt.colorbar(cs1, ax=ax1, label='Brightness (R)')
    
    ax2 = axes[1]
    cs2 = ax2.contourf(X, Z, s3_emission, levels=100, cmap='viridis')
    ax2.contour(X, Z, s3_emission, levels=10, colors='black', alpha=0.3, linewidths=0.5)
    ax2.set_xlabel('X Position (R_J)', fontsize=12)
    ax2.set_ylabel('Z Position (R_J)', fontsize=12)
    ax2.set_title('S$^{++}$ Emission (680 Å ± 3σ)', fontsize=14)
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)
    plt.colorbar(cs2, ax=ax2, label='Brightness (R)')
    
    ax3 = axes[2]
    cs3 = ax3.contourf(X, Z, o2_emission, levels=100, cmap='viridis')
    ax3.contour(X, Z, o2_emission, levels=10, colors='black', alpha=0.3, linewidths=0.5)
    ax3.set_xlabel('X Position (R_J)', fontsize=12)
    ax3.set_ylabel('Z Position (R_J)', fontsize=12)
    ax3.set_title('O$^{+}$ Emission (833 Å ± 3σ)', fontsize=14)
    ax3.set_aspect('equal')
    ax3.grid(True, alpha=0.3)
    plt.colorbar(cs3, ax=ax3, label='Brightness (R)')
    
    plt.suptitle('UV Emission Maps: View Along +Y Direction (Parallelized Low Resolution)', 
                 fontsize=16, y=1.02)
    plt.tight_layout()
    plt.savefig('emission_contours_low_res_parallel.png', dpi=150, bbox_inches='tight')
    print(f"\nSaved: emission_contours_low_res_parallel.png")
    
    # Individual plots
    print("\nCreating individual contour plots...")
    
    fig1, ax = plt.subplots(figsize=(8, 6))
    cs = ax.contourf(X, Z, total_emission, levels=100, cmap='viridis')
    ax.contour(X, Z, total_emission, levels=10, colors='black', alpha=0.3, linewidths=0.5)
    ax.set_xlabel('X Position (R_J)', fontsize=12)
    ax.set_ylabel('Z Position (R_J)', fontsize=12)
    ax.set_title('Total UV Emission (550-2100 Å) - Parallelized Low Resolution', fontsize=14)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    plt.colorbar(cs, ax=ax, label='Brightness (R)')
    plt.tight_layout()
    plt.savefig('total_emission_low_res_parallel.png', dpi=150)
    print("  Saved: total_emission_low_res_parallel.png")
    
    fig2, ax = plt.subplots(figsize=(8, 6))
    cs = ax.contourf(X, Z, s3_emission, levels=100, cmap='viridis')
    ax.contour(X, Z, s3_emission, levels=10, colors='black', alpha=0.3, linewidths=0.5)
    ax.set_xlabel('X Position (R_J)', fontsize=12)
    ax.set_ylabel('Z Position (R_J)', fontsize=12)
    ax.set_title('S$^{++}$ Emission (680 Å ± 3σ) - Parallelized Low Resolution', fontsize=14)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    plt.colorbar(cs, ax=ax, label='Brightness (R)')
    plt.tight_layout()
    plt.savefig('s3_emission_low_res_parallel.png', dpi=150)
    print("  Saved: s3_emission_low_res_parallel.png")
    
    fig3, ax = plt.subplots(figsize=(8, 6))
    cs = ax.contourf(X, Z, o2_emission, levels=100, cmap='viridis')
    ax.contour(X, Z, o2_emission, levels=10, colors='black', alpha=0.3, linewidths=0.5)
    ax.set_xlabel('X Position (R_J)', fontsize=12)
    ax.set_ylabel('Z Position (R_J)', fontsize=12)
    ax.set_title('O$^{+}$ Emission (833 Å ± 3σ) - Parallelized Low Resolution', fontsize=14)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    plt.colorbar(cs, ax=ax, label='Brightness (R)')
    plt.tight_layout()
    plt.savefig('o2_emission_low_res_parallel.png', dpi=150)
    print("  Saved: o2_emission_low_res_parallel.png")
    
    plt.show()
    
    # =========================================================================
    # SUMMARY
    # =========================================================================
    
    print("\n" + "=" * 70)
    print("Parallelized Emission Frame Generation Complete")
    print("=" * 70)
    
    print("\nPhysical Interpretation:")
    print("- Peak emission occurs near 6 R_J (Io's orbital radius)")
    print("- Vertical extent shows plasma scale height of ~0.5-1 R_J")
    print("- S++ and O+ show different spatial distributions")
    print("- Emission drops rapidly inside 5 R_J and outside 8 R_J")
    
    print("\nNumerical Methods:")
    print("- Ray tracing: Uniform step through 3D domain")
    print("- Plasma interpolation: Trilinear in 3D Cartesian")
    print("- Emission tables: Bilinear in log10(ne)-log10(Te)")
    print("- LOS integration: Simpson's rule")
    print("- Spectral convolution: ERF-based Gaussian")
    print("- Parallelization: Shared memory across workers")


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    if os.name == 'nt':
        mp.freeze_support()
    
    main()