#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
basic_example_optical_cm3_emission_model_use_tables.py

Example Script for Optical Emission Modeling of Io Plasma Torus
================================================================

This script demonstrates the use of the IPT_cm3_emission_model library for calculating
optical emission spectra from the Io Plasma Torus. It shows typical usage for both single
and double Maxwellian electron distributions with parameters appropriate for ground-based
optical spectroscopy.

The example uses emission tables pre-calculated from CHIANTI 11.0.2 atomic database
and demonstrates:
1. Loading emission tables for single and double Maxwellian distributions
2. Calculating discrete emission line brightnesses in optical wavelengths
3. Convolving with ground-based telescope instrument response functions
4. Analyzing key diagnostic optical emission lines
5. Creating publication-quality plots for optical spectroscopy

REQUIRED DATA FILES:
- Single Maxwellian: CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5
- Double Maxwellian: CHIANTI_11.0.2_emiss_tables_double_maxwellian_24x10x24x12.h5
- Located in ../Emiss_tables/ directory relative to this script

OUTPUT:
- Optical emission spectra for single and double Maxwellian cases
- Detailed analysis of key optical emission lines ([S II], [O I], [S III])
- Publication-quality plots saved as PNG files

WAVELENGTH RANGE:
- Optical: 3000-10000 Å (suitable for ground-based telescopes)

AUTHOR: Edward (Eddie) G. Nerney
INSTITUTION: Laboratory for Atmospheric and Space Physics, University of Colorado Boulder
VERSION: 1.0
DATE: November 2025
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import time
from scipy.integrate import simpson

# Import the IPT emission model library
from IPT_cm3_emission_model import (
    EmissionTables,
    calculate_ipt_emiss_tables_single,
    calculate_ipt_emiss_tables_double,
    simulate_ipt_spectrum_rayleighs_erf_form,
    analyze_key_optical_emission_lines
)


def main():
    """
    Example demonstrating optical emission modeling for Io Plasma Torus.
    
    This function shows typical usage of the emission table infrastructure for
    calculating IPT optical spectra observable from ground-based telescopes.
    It includes both single and double Maxwellian cases with appropriate
    parameters for optical spectroscopy.
    
    Returns
    -------
    results : dict
        Dictionary containing:
        - 'xwav': wavelength grid
        - 'ypts_single_maxwellian': single Maxwellian spectrum
        - 'ypts_double_maxwellian': double Maxwellian spectrum
        - 'xwavi_single': discrete line wavelengths (single)
        - 'yptsi_single': discrete line brightnesses (single)
        - 'xwavi_double': discrete line wavelengths (double)
        - 'yptsi_double': discrete line brightnesses (double)
        - 'tables': tables object for further use
        - 'line_analysis': results from optical emission line analysis
    
    Example Usage
    -------------
    >>> results = main()
    >>> plt.plot(results['xwav'], results['ypts_single_maxwellian'])
    >>> plt.xlabel('Wavelength [Å]')
    >>> plt.ylabel('Brightness [R/Å]')
    >>> plt.show()
    
    Notes
    -----
    The example uses typical IPT parameters from observations:
    - Core temperature: 5 eV (from UV/optical line ratios)
    - Density: 2200 cm^-3 (peak torus density)
    - Column densities: ~6 R_J path length through torus
    - Hot population: 270 eV, 0.25% fraction (from in situ data)
    
    Optical wavelength range and resolution parameters are set for
    typical ground-based spectroscopy with R~3000-5000.
    """
    # ============================================================================
    # WAVELENGTH GRID SETUP - OPTICAL RANGE
    # ============================================================================
    # Define output wavelength grid for convolved spectrum
    # Ground-based optical telescope wavelength range and bin width
    min_xwav = 3000.0    # Minimum wavelength [Angstroms] - blue optical
    max_xwav = 10000.0   # Maximum wavelength [Angstroms] - near-infrared limit
    xwav_bin_width = 2.0  # Wavelength bin width [Angstroms] - typical for R~3000
    
    # Create wavelength grid
    n_wav = int((max_xwav - min_xwav) / xwav_bin_width) + 1
    xwav = np.linspace(min_xwav, max_xwav, n_wav)
    
    # Instrument resolution for ground-based optical spectroscopy
    fwhm = 3.0  # FWHM [Angstroms] - typical for R~3000 ground-based spectrograph
    
    # ============================================================================
    # PLASMA PARAMETERS
    # ============================================================================
    # Core/cold electron population parameters
    Tec = 5.0          # Core temperature [eV] - typical IPT value
    nec = 2200.0       # Core density [cm^-3] - peak torus density
    
    # Hot electron population parameters (for double Maxwellian)
    Teh = 270.0        # Hot temperature [eV] - from Voyager measurements
    feh = 0.0025       # Hot electron fraction (0.25% typical)
    fec = 1.0 - feh    # Core electron fraction
    neh = nec * (1.0/fec - 1.0)  # Hot density [cm^-3]
    ne_total = nec + neh          # Total electron density [cm^-3]
    
    # Column densities [cm^-2]
    # Based on Nerney et al. (2017), Steffl et al. (2004b), Thomas et al. (2004)
    # These values correspond to ~10 R_J path length through the Io torus
    column_densities = {
        'S+': 1.2e13,    # S II - singly ionized sulfur
        'S++': 4.2e13,   # S III - dominant sulfur ion (~90% of S)
        'S+++': 5.92e12, # S IV - minor component
        'S++++': 6.0e11, # S V - highly ionized, trace amounts
        'O+': 5.2e13,    # O II - dominant ion overall (~80% of O)
        'O++': 5.92e12   # O III - secondary oxygen ion
    }
    
    # ============================================================================
    # SINGLE MAXWELLIAN CALCULATION
    # ============================================================================
    print("="*64)
    print("SINGLE MAXWELLIAN CALCULATION - OPTICAL WAVELENGTHS")
    print("="*64)
    
    # Initialize emission tables object
    tables = EmissionTables()
    
    # Load single Maxwellian tables
    # Tables should be in ../Emiss_tables/ directory relative to this script
    emiss_dir = Path(__file__).parent.parent / 'Emiss_tables'
    tables.load_single_maxwellian_tables(
        str(emiss_dir / 'CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5')
    )
    
    # Calculate discrete emission lines in optical range
    xwavi_single, yptsi_single = calculate_ipt_emiss_tables_single(
        tables, Tec, nec, column_densities, min_wav=min_xwav, max_wav=max_xwav
    )
    
    # ============================================================================
    # DOUBLE MAXWELLIAN CALCULATION
    # ============================================================================
    print()
    print("="*64)
    print("DOUBLE MAXWELLIAN CALCULATION - OPTICAL WAVELENGTHS")
    print("="*64)
    
    # Load double Maxwellian tables
    tables.load_double_maxwellian_tables(
        str(emiss_dir / 'CHIANTI_11.0.2_emiss_tables_double_maxwellian_24x10x24x12.h5')
    )
    
    # Calculate discrete emission lines in optical range
    xwavi_double, yptsi_double = calculate_ipt_emiss_tables_double(
        tables, Tec, Teh, ne_total, feh, column_densities, 
        min_wav=min_xwav, max_wav=max_xwav
    )
    
    # ============================================================================
    # CONVOLVE WITH INSTRUMENT RESPONSE
    # ============================================================================
    # Create realistic spectra by convolving discrete lines with instrument PSF
    print()
    print("Convolving with ground-based telescope instrument response...")
    
    ypts_single_maxwellian = simulate_ipt_spectrum_rayleighs_erf_form(
        xwav, xwav_bin_width, xwavi_single, yptsi_single, fwhm=fwhm
    )
    
    ypts_double_maxwellian = simulate_ipt_spectrum_rayleighs_erf_form(
        xwav, xwav_bin_width, xwavi_double, yptsi_double, fwhm=fwhm
    )
    
    # ============================================================================
    # SUMMARY OUTPUT
    # ============================================================================
    print()
    print("="*64)
    print("OPTICAL SIMULATION COMPLETE")
    print("="*64)
    print("Plasma parameters:")
    print(f"  Single Maxwellian: Te = {Tec:.1f} eV, ne = {nec:.0f} cm^-3")
    print(f"  Double Maxwellian: Tec = {Tec:.1f} eV, Teh = {Teh:.0f} eV")
    print(f"                     feh = {feh:.4f}, ne_total = {ne_total:.0f} cm^-3")
    print()
    print("Column densities [cm^-2]:")
    for ion, value in column_densities.items():
        print(f"  {ion:<5}: {value:.2e}")
    print()
    print("Optical wavelength range:")
    print(f"  {min_xwav:.0f} - {max_xwav:.0f} Å")
    print(f"  Spectral resolution (FWHM): {fwhm:.1f} Å")
    print(f"  Bin width: {xwav_bin_width:.1f} Å")
    print()
    print("Number of emission lines in optical wavelength range:")
    print(f"  Single Maxwellian: {len(xwavi_single)}")
    print(f"  Double Maxwellian: {len(xwavi_double)}")
    print()
    print("Output optical spectrum statistics:")
    print(f"  Single Maxwellian peak: {ypts_single_maxwellian.max():.2f} R/Å")
    print(f"  Double Maxwellian peak: {ypts_double_maxwellian.max():.2f} R/Å")
    print(f"  Integrated brightness (single): {simpson(ypts_single_maxwellian, x=xwav):.0f} R")
    print(f"  Integrated brightness (double): {simpson(ypts_double_maxwellian, x=xwav):.0f} R")
    
    # ============================================================================
    # PLOTTING - OPTICAL SPECTRA
    # ============================================================================
    print()
    print("Creating optical spectrum plots...")
    
    # Create figure with two subplots
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))
    
    # Single Maxwellian plot - using matplotlib default color C0
    axes[0].plot(xwav, ypts_single_maxwellian, 'C0-', linewidth=1.5)
    axes[0].set_xlabel('Wavelength [Å]', fontsize=12)
    axes[0].set_ylabel('Brightness [R/Å]', fontsize=12)
    axes[0].set_title(
        f'IPT Optical Emission: Single Maxwellian (Te={Tec:.1f} eV, ne={nec:.0f} cm⁻³)', 
        fontsize=14
    )
    axes[0].grid(True, alpha=0.3)
    axes[0].set_xlim(min_xwav, max_xwav)
    axes[0].set_ylim(bottom=0)
    
    # Double Maxwellian plot - using matplotlib default color C1
    axes[1].plot(xwav, ypts_double_maxwellian, 'C1-', linewidth=1.5)
    axes[1].set_xlabel('Wavelength [Å]', fontsize=12)
    axes[1].set_ylabel('Brightness [R/Å]', fontsize=12)
    axes[1].set_title(
        f'IPT Optical Emission: Double Maxwellian (Tec={Tec:.1f} eV, Teh={Teh:.0f} eV, feh={feh:.4f})', 
        fontsize=14
    )
    axes[1].grid(True, alpha=0.3)
    axes[1].set_xlim(min_xwav, max_xwav)
    axes[1].set_ylim(bottom=0)
    
    plt.tight_layout()
    plt.savefig('ipt_optical_emission_single_vs_double.png', dpi=300)
    print("  Saved: ipt_optical_emission_single_vs_double.png")
    
    # Comparison plot
    fig2, ax = plt.subplots(figsize=(12, 6))
    
    ax.plot(xwav, ypts_single_maxwellian, 'C0-', linewidth=1.5, 
            label='Single Maxwellian', alpha=0.8)
    ax.plot(xwav, ypts_double_maxwellian, 'C1-', linewidth=1.5, 
            label='Double Maxwellian', alpha=0.8)
    ax.set_xlabel('Wavelength [Å]', fontsize=12)
    ax.set_ylabel('Brightness [R/Å]', fontsize=12)
    ax.set_title('IPT Optical Emission Comparison: Single vs Double Maxwellian', fontsize=14)
    ax.legend(loc='upper right', fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(min_xwav, max_xwav)
    ax.set_ylim(bottom=0)
    
    plt.tight_layout()
    plt.savefig('ipt_optical_emission_comparison.png', dpi=300)
    print("  Saved: ipt_optical_emission_comparison.png")
    
    # Zoom plots for key spectral regions
    # [S II] doublet region (6700-6750 Å)
    fig3, ax = plt.subplots(figsize=(10, 6))
    mask_sii = (xwav >= 6700) & (xwav <= 6750)
    ax.plot(xwav[mask_sii], ypts_single_maxwellian[mask_sii], 'C0-', 
            linewidth=1.5, label='Single Maxwellian', alpha=0.8)
    ax.plot(xwav[mask_sii], ypts_double_maxwellian[mask_sii], 'C1-', 
            linewidth=1.5, label='Double Maxwellian', alpha=0.8)
    ax.set_xlabel('Wavelength [Å]', fontsize=12)
    ax.set_ylabel('Brightness [R/Å]', fontsize=12)
    ax.set_title('[S II] 6716/6731 Å Doublet Region', fontsize=14)
    ax.legend(loc='upper right', fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(6700, 6750)
    
    plt.tight_layout()
    plt.savefig('ipt_optical_sii_doublet_region.png', dpi=300)
    print("  Saved: ipt_optical_sii_doublet_region.png")
    
    # Red line region (6280-6330 Å) - [O I] 6300 and [S III] 6312
    fig4, ax = plt.subplots(figsize=(10, 6))
    mask_red = (xwav >= 6280) & (xwav <= 6330)
    ax.plot(xwav[mask_red], ypts_single_maxwellian[mask_red], 'C0-', 
            linewidth=1.5, label='Single Maxwellian', alpha=0.8)
    ax.plot(xwav[mask_red], ypts_double_maxwellian[mask_red], 'C1-', 
            linewidth=1.5, label='Double Maxwellian', alpha=0.8)
    ax.set_xlabel('Wavelength [Å]', fontsize=12)
    ax.set_ylabel('Brightness [R/Å]', fontsize=12)
    ax.set_title('[O I] 6300 Å and [S III] 6312 Å Region', fontsize=14)
    ax.legend(loc='upper right', fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(6280, 6330)
    
    plt.tight_layout()
    plt.savefig('ipt_optical_red_line_region.png', dpi=300)
    print("  Saved: ipt_optical_red_line_region.png")
    
    plt.show()
    
    # ================================================================
    # ANALYZE SPECIFIC OPTICAL EMISSION LINES
    # ================================================================
    # Perform detailed analysis of key diagnostic optical emission lines
    # and generate high-resolution plots for scientific interpretation
    line_analysis_results = analyze_key_optical_emission_lines(
        xwav,
        ypts_single_maxwellian,
        ypts_double_maxwellian,
        fwhm_instrument=fwhm  
    )
    
    # Return results for further analysis
    return {
        'xwav': xwav,
        'ypts_single_maxwellian': ypts_single_maxwellian,
        'ypts_double_maxwellian': ypts_double_maxwellian,
        'xwavi_single': xwavi_single,
        'yptsi_single': yptsi_single,
        'xwavi_double': xwavi_double,
        'yptsi_double': yptsi_double,
        'tables': tables,  # Return tables object for further use
        'line_analysis': line_analysis_results
    }


if __name__ == "__main__":
    # Run the optical emission calculation
    print("="*64)
    print("OPTICAL EMISSION MODEL FOR IO PLASMA TORUS")
    print("Ground-Based Telescope Spectroscopy")
    print("Table-Based Implementation using CHIANTI 11.0.2")
    print("="*64)
    print()
    
    start_time = time.perf_counter()
    
    results = main()
    
    end_time = time.perf_counter()
    
    elapsed_time = end_time - start_time
    print(f"Execution time: {elapsed_time:.4f} seconds")
    
    print()
    print("="*64)
    print("Results stored in 'results' dictionary")
    print("Available for further analysis and visualization")
    print("="*64)