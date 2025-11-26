#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
basic_example_optical_integration_emission_model_use_tables.py

Example Optical Emission Calculation Using Line-of-Sight Integration
=====================================================================

This script demonstrates how to use the Jovian UV/Optical emission raytracer to
calculate synthetic optical spectra through the Io Plasma Torus using proper
line-of-sight integration with species-by-species interpolation from CHIANTI
emission tables.

The example calculates optical emission spectra for test cases with both single
and double Maxwellian electron distributions:
1. Single Maxwellian: equatorial and off-equator lines of sight
2. Double Maxwellian: same geometries with hot electron population

All lines of sight start at x=6 R_J and look through the plasma torus in the
+y direction, sampling the peak emission region near Io's orbit.

PHYSICAL MODEL:
- Ray tracing through 3D plasma model with trilinear interpolation
- Vectorized bilinear interpolation of single Maxwellian emission rates (2D tables)
- Vectorized quadrilinear interpolation of double Maxwellian emission rates (4D tables)
- Per-ion photon emission rates multiplied by ion density and 1e-6 Rayleigh factor
- Simpson's rule integration over line of sight
- Analytic ERF-based Gaussian convolution for instrument response

WAVELENGTH RANGE:
- Optical: 3000-10000 Å (suitable for ground-based telescopes)
- Key diagnostic lines: [S II] 6716/6731, [O III] 5007, [S III] 6312

DIRECTORY STRUCTURE EXPECTED:
- IPT_Emission_MOP_Community_Code/
  - LOS_Integration/           (this script and module)
  - Emiss_tables/              (CHIANTI emission tables)
  - 3D_Torus_Model/            (3D plasma model in HDF5 format)

OUTPUT FILES:
- optical_spectrum_single_equatorial.png: Single Maxwellian equatorial spectrum
- optical_spectrum_single_off_equator.png: Single Maxwellian off-equator spectrum
- optical_spectrum_single_comparison.png: Comparison of equatorial vs off-equator
- optical_spectrum_double_equatorial.png: Double Maxwellian equatorial spectrum
- optical_spectrum_double_off_equator.png: Double Maxwellian off-equator spectrum
- optical_spectrum_single_vs_double.png: Comparison of single vs double Maxwellian
- optical_sii_doublet_region.png: Zoom on [S II] 6716/6731 Å doublet
- optical_red_line_region.png: Zoom on [S III] 6312 Å and vicinity

AUTHOR: Edward (Eddie) G. Nerney
INSTITUTION: Laboratory for Atmospheric and Space Physics, University of Colorado Boulder
DATE: November 2025
LICENSE: Open source for academic and research use
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import sys
import time
from pathlib import Path

# Add the current directory to path to import the raytracer module
sys.path.insert(0, str(Path(__file__).parent))

# Import the raytracer module
from IPT_emiss_MOP_community_code import JovianUVEmissionRaytracer


def main():
    """
    Main function demonstrating raytracer usage for optical wavelengths.
    
    This function demonstrates the complete workflow for calculating optical
    emission spectra through the Io Plasma Torus using ground-based telescope
    parameters:
    
    1. Initialize the raytracer with plasma model and emission tables
    2. Calculate optical spectra for single Maxwellian distribution
    3. Calculate optical spectra for double Maxwellian distribution  
    4. Create comparison plots showing hot electron enhancement
    5. Generate zoom plots for key diagnostic line regions
    6. Print diagnostic information about plasma conditions
    
    The plasma model and emission tables are loaded from default locations
    based on the expected directory structure.
    
    Returns
    -------
    results : dict
        Dictionary containing all calculated spectra and line lists
        
    Notes
    -----
    Optical wavelength range and resolution parameters are set for typical
    ground-based spectroscopy with R~3000-5000. Key diagnostic lines include:
    
    - [S II] 6716/6731 Å: Density diagnostic doublet
    - [S III] 6312 Å: Temperature diagnostic  
    - [O III] 4959/5007 Å: Nebular doublet
    - [O II] 3726/3729 Å: Density diagnostic doublet
    """
    print("="*70)
    print("Jovian Optical Emission Raytracer - Line-of-Sight Integration")
    print("="*70)
    print()
    print("This example demonstrates optical emission calculations through")
    print("the Io Plasma Torus using proper line-of-sight integration")
    print("with species-by-species interpolation from CHIANTI tables.")
    print()
    print("Wavelength range: 3000-10000 Å (ground-based optical)")
    print()
    
    start_time = time.perf_counter()
    
    # ========================================================================
    # INITIALIZE RAYTRACER
    # ========================================================================
    
    print("Loading plasma model and emission tables...")
    try:
        raytracer = JovianUVEmissionRaytracer()
    except FileNotFoundError as e:
        print(f"\nError: {e}")
        print("\nPlease ensure the following directory structure exists:")
        print("  IPT_Emission_MOP_Community_Code/")
        print("    - LOS_Integration/ (containing this script)")
        print("    - Emiss_tables/CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5")
        print("    - Emiss_tables/CHIANTI_11.0.2_emiss_tables_double_maxwellian_24x10x24x12.h5")
        print("    - 3D_Torus_Model/jovian_plasma_interpolated_381x381x231.h5")
        return None
    
    # ========================================================================
    # DEFINE OBSERVATION GEOMETRY AND OPTICAL PARAMETERS
    # ========================================================================
    
    # Optical wavelength parameters (ground-based telescope)
    wavelength_range = (3000, 10000)  # Optical range [Angstroms]
    bin_width = 2.0  # Spectral bin width [Angstroms] - typical for R~3000
    fwhm = 3.0  # Instrumental FWHM [Angstroms] - ground-based spectrograph
    ds = 0.01  # Integration step size [R_J]
    
    # Test case 1: Equatorial line of sight
    slit_pos_vec1 = np.array([6.0, -20.0, 0.0])  # Start at x=6 R_J, equator
    norm_vec1 = np.array([0.0, 1.0, 0.0])        # Look in +y direction
    
    # Test case 2: Off-equator line of sight
    slit_pos_vec2 = np.array([6.0, -20.0, 0.5])  # 0.5 R_J above equator
    norm_vec2 = np.array([0.0, 1.0, 0.0])        # Look in +y direction
    
    # ========================================================================
    # SINGLE MAXWELLIAN CALCULATIONS
    # ========================================================================
    
    print("\n" + "="*70)
    print("Test Case 1: Equatorial Line of Sight (Single Maxwellian)")
    print("-"*70)
    print(f"Starting position: ({slit_pos_vec1[0]:.1f}, {slit_pos_vec1[1]:.1f}, {slit_pos_vec1[2]:.1f}) R_J")
    print(f"Direction vector: ({norm_vec1[0]:.1f}, {norm_vec1[1]:.1f}, {norm_vec1[2]:.1f})")
    print(f"Wavelength range: {wavelength_range[0]:.0f} - {wavelength_range[1]:.0f} Å")
    print()
    
    wave_bins1_single, spectrum1_single, lines1_single = raytracer.calculate_spectrum_single(
        slit_pos_vec1, norm_vec1,
        wavelength_range=wavelength_range,
        bin_width=bin_width,
        fwhm=fwhm,
        ds=ds
    )
    
    # Calculate summary statistics
    total_brightness1_single = simpson(spectrum1_single, x=wave_bins1_single)
    peak_idx1_single = np.argmax(spectrum1_single)
    
    print(f"\nSingle Maxwellian Results:")
    print(f"  Total brightness: {total_brightness1_single:.1f} Rayleighs")
    print(f"  Peak brightness: {spectrum1_single.max():.2f} R/Å at {wave_bins1_single[peak_idx1_single]:.1f} Å")
    print(f"  Number of emission lines: {len(lines1_single)}")
    
    # Test case 2: Off-equator (single Maxwellian)
    print("\n" + "="*70)
    print("Test Case 2: Off-Equator Line of Sight (Single Maxwellian)")
    print("-"*70)
    print(f"Starting position: ({slit_pos_vec2[0]:.1f}, {slit_pos_vec2[1]:.1f}, {slit_pos_vec2[2]:.1f}) R_J")
    print(f"Direction vector: ({norm_vec2[0]:.1f}, {norm_vec2[1]:.1f}, {norm_vec2[2]:.1f})")
    print()
    
    wave_bins2_single, spectrum2_single, lines2_single = raytracer.calculate_spectrum_single(
        slit_pos_vec2, norm_vec2,
        wavelength_range=wavelength_range,
        bin_width=bin_width,
        fwhm=fwhm,
        ds=ds
    )
    
    total_brightness2_single = simpson(spectrum2_single, x=wave_bins2_single)
    peak_idx2_single = np.argmax(spectrum2_single)
    
    print(f"\nSingle Maxwellian Results:")
    print(f"  Total brightness: {total_brightness2_single:.1f} Rayleighs")
    print(f"  Peak brightness: {spectrum2_single.max():.2f} R/Å at {wave_bins2_single[peak_idx2_single]:.1f} Å")
    
    # ========================================================================
    # DOUBLE MAXWELLIAN CALCULATIONS
    # ========================================================================
    
    # Initialize variables for double Maxwellian results
    wave_bins1_double = None
    spectrum1_double = None
    lines1_double = None
    wave_bins2_double = None
    spectrum2_double = None
    lines2_double = None
    total_brightness1_double = 0.0
    total_brightness2_double = 0.0
    enhancement1 = 1.0
    enhancement2 = 1.0
    
    if raytracer.double_maxwellian_loaded:
        # Test case 3: Equatorial (double Maxwellian)
        print("\n" + "="*70)
        print("Test Case 3: Equatorial LOS (Double Maxwellian)")
        print("-"*70)
        
        wave_bins1_double, spectrum1_double, lines1_double = raytracer.calculate_spectrum_double(
            slit_pos_vec1, norm_vec1,
            wavelength_range=wavelength_range,
            bin_width=bin_width,
            fwhm=fwhm,
            ds=ds
        )
        
        total_brightness1_double = simpson(spectrum1_double, x=wave_bins1_double)
        peak_idx1_double = np.argmax(spectrum1_double)
        enhancement1 = total_brightness1_double / total_brightness1_single if total_brightness1_single > 0 else 1.0
        
        print(f"\nDouble Maxwellian Results:")
        print(f"  Total brightness: {total_brightness1_double:.1f} Rayleighs")
        print(f"  Peak brightness: {spectrum1_double.max():.2f} R/Å at {wave_bins1_double[peak_idx1_double]:.1f} Å")
        print(f"  Enhancement factor: {enhancement1:.3f}")
        
        # Test case 4: Off-equator (double Maxwellian)
        print("\n" + "="*70)
        print("Test Case 4: Off-Equator LOS (Double Maxwellian)")
        print("-"*70)
        
        wave_bins2_double, spectrum2_double, lines2_double = raytracer.calculate_spectrum_double(
            slit_pos_vec2, norm_vec2,
            wavelength_range=wavelength_range,
            bin_width=bin_width,
            fwhm=fwhm,
            ds=ds
        )
        
        total_brightness2_double = simpson(spectrum2_double, x=wave_bins2_double)
        peak_idx2_double = np.argmax(spectrum2_double)
        enhancement2 = total_brightness2_double / total_brightness2_single if total_brightness2_single > 0 else 1.0
        
        print(f"\nDouble Maxwellian Results:")
        print(f"  Total brightness: {total_brightness2_double:.1f} Rayleighs")
        print(f"  Peak brightness: {spectrum2_double.max():.2f} R/Å at {wave_bins2_double[peak_idx2_double]:.1f} Å")
        print(f"  Enhancement factor: {enhancement2:.3f}")
    
    # ========================================================================
    # GENERATE PLOTS
    # ========================================================================
    
    print("\n" + "="*70)
    print("Generating Optical Spectrum Plots")
    print("-"*70)
    
    # Plot 1: Single Maxwellian equatorial
    fig1 = raytracer.plot_spectrum(
        wave_bins1_single, spectrum1_single, lines1_single,
        title="Single Maxwellian: Equatorial LOS - Optical (x=6 R_J, z=0)",
        color='C0'
    )
    plt.savefig('optical_spectrum_single_equatorial.png', dpi=150)
    print("Saved: optical_spectrum_single_equatorial.png")
    
    # Plot 2: Single Maxwellian off-equator
    fig2 = raytracer.plot_spectrum(
        wave_bins2_single, spectrum2_single, lines2_single,
        title="Single Maxwellian: Off-Equator LOS - Optical (x=6 R_J, z=0.5 R_J)",
        color='C0'
    )
    plt.savefig('optical_spectrum_single_off_equator.png', dpi=150)
    print("Saved: optical_spectrum_single_off_equator.png")
    
    # Plot 3: Comparison - equatorial vs off-equator (single Maxwellian)
    fig3, ax3 = plt.subplots(figsize=(12, 6))
    ax3.plot(wave_bins1_single, spectrum1_single, 'C0', 
             label='Equatorial (z=0)', linewidth=1.0)
    ax3.plot(wave_bins2_single, spectrum2_single, 'C1', 
             label='Off-equator (z=0.5 R_J)', linewidth=1.0, alpha=0.7)
    ax3.set_xlabel('Wavelength [Å]', fontsize=12)
    ax3.set_ylabel('Brightness [Rayleighs/Å]', fontsize=12)
    ax3.set_title('Single Maxwellian: Equatorial vs Off-Equator Optical Spectra', fontsize=14)
    ax3.grid(True, alpha=0.3)
    ax3.legend(loc='upper right')
    ax3.set_xlim(wavelength_range[0], wavelength_range[1])
    ax3.set_ylim(bottom=0)
    
    ax3.text(0.02, 0.98, f'Equatorial: {total_brightness1_single:.1f} R\n'
                        f'Off-equator: {total_brightness2_single:.1f} R',
             transform=ax3.transAxes, va='top', ha='left',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig('optical_spectrum_single_comparison.png', dpi=150)
    print("Saved: optical_spectrum_single_comparison.png")
    
    # Double Maxwellian plots (if available)
    if raytracer.double_maxwellian_loaded and spectrum1_double is not None:
        # Plot 4: Double Maxwellian equatorial
        fig4 = raytracer.plot_spectrum(
            wave_bins1_double, spectrum1_double, lines1_double,
            title="Double Maxwellian: Equatorial LOS - Optical",
            color='C1'
        )
        plt.savefig('optical_spectrum_double_equatorial.png', dpi=150)
        print("Saved: optical_spectrum_double_equatorial.png")
        
        # Plot 5: Double Maxwellian off-equator
        fig5 = raytracer.plot_spectrum(
            wave_bins2_double, spectrum2_double, lines2_double,
            title="Double Maxwellian: Off-Equator LOS - Optical",
            color='C1'
        )
        plt.savefig('optical_spectrum_double_off_equator.png', dpi=150)
        print("Saved: optical_spectrum_double_off_equator.png")
        
        # Plot 6: Comparison - single vs double Maxwellian (equatorial)
        fig6, ax6 = plt.subplots(figsize=(12, 6))
        ax6.plot(wave_bins1_single, spectrum1_single, 'C0', 
                 label='Single Maxwellian', linewidth=1.0)
        ax6.plot(wave_bins1_double, spectrum1_double, 'C1', 
                 label='Double Maxwellian', linewidth=1.0, alpha=0.7)
        ax6.set_xlabel('Wavelength [Å]', fontsize=12)
        ax6.set_ylabel('Brightness [Rayleighs/Å]', fontsize=12)
        ax6.set_title('Equatorial LOS: Single vs Double Maxwellian - Optical', fontsize=14)
        ax6.grid(True, alpha=0.3)
        ax6.legend(loc='upper right')
        ax6.set_xlim(wavelength_range[0], wavelength_range[1])
        ax6.set_ylim(bottom=0)
        
        ax6.text(0.02, 0.98, f'Single: {total_brightness1_single:.1f} R\n'
                            f'Double: {total_brightness1_double:.1f} R\n'
                            f'Enhancement: {enhancement1:.3f}×',
                 transform=ax6.transAxes, va='top', ha='left',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        plt.savefig('optical_spectrum_single_vs_double.png', dpi=150)
        print("Saved: optical_spectrum_single_vs_double.png")
    
    # ========================================================================
    # ZOOM PLOTS FOR KEY DIAGNOSTIC LINE REGIONS
    # ========================================================================
    
    print("\nGenerating zoom plots for key diagnostic regions...")
    
    # [S II] doublet region (6700-6750 Å)
    fig7, ax7 = plt.subplots(figsize=(10, 6))
    mask_sii = (wave_bins1_single >= 6700) & (wave_bins1_single <= 6750)
    ax7.plot(wave_bins1_single[mask_sii], spectrum1_single[mask_sii], 'C0', 
            linewidth=1.5, label='Single Maxwellian', alpha=0.8)
    if spectrum1_double is not None:
        ax7.plot(wave_bins1_double[mask_sii], spectrum1_double[mask_sii], 'C1', 
                linewidth=1.5, label='Double Maxwellian', alpha=0.8)
    ax7.set_xlabel('Wavelength [Å]', fontsize=12)
    ax7.set_ylabel('Brightness [Rayleighs/Å]', fontsize=12)
    ax7.set_title('[S II] 6716/6731 Å Doublet Region - Density Diagnostic', fontsize=14)
    ax7.legend(loc='upper right', fontsize=11)
    ax7.grid(True, alpha=0.3)
    ax7.set_xlim(6700, 6750)
    ax7.set_ylim(bottom=0)
    
    # Add line identification annotations
    ax7.axvline(6716.4, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    ax7.axvline(6730.8, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    ax7.text(6716.4, ax7.get_ylim()[1]*0.9, '[S II]\n6716 Å', 
             ha='center', va='top', fontsize=9)
    ax7.text(6730.8, ax7.get_ylim()[1]*0.9, '[S II]\n6731 Å', 
             ha='center', va='top', fontsize=9)
    
    plt.tight_layout()
    plt.savefig('optical_sii_doublet_region.png', dpi=300)
    print("Saved: optical_sii_doublet_region.png")
    
    # Red line region (6280-6330 Å) - [S III] 6312
    fig8, ax8 = plt.subplots(figsize=(10, 6))
    mask_red = (wave_bins1_single >= 6280) & (wave_bins1_single <= 6330)
    ax8.plot(wave_bins1_single[mask_red], spectrum1_single[mask_red], 'C0', 
            linewidth=1.5, label='Single Maxwellian', alpha=0.8)
    if spectrum1_double is not None:
        ax8.plot(wave_bins1_double[mask_red], spectrum1_double[mask_red], 'C1', 
                linewidth=1.5, label='Double Maxwellian', alpha=0.8)
    ax8.set_xlabel('Wavelength [Å]', fontsize=12)
    ax8.set_ylabel('Brightness [Rayleighs/Å]', fontsize=12)
    ax8.set_title('[S III] 6312 Å Region - Temperature Diagnostic', fontsize=14)
    ax8.legend(loc='upper right', fontsize=11)
    ax8.grid(True, alpha=0.3)
    ax8.set_xlim(6280, 6330)
    ax8.set_ylim(bottom=0)
    
    # Add line identification
    ax8.axvline(6312.1, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    ax8.text(6312.1, ax8.get_ylim()[1]*0.9, '[S III]\n6312 Å', 
             ha='center', va='top', fontsize=9)
    
    plt.tight_layout()
    plt.savefig('optical_red_line_region.png', dpi=300)
    print("Saved: optical_red_line_region.png")
    
    # [O III] nebular doublet region (4940-5020 Å)
    fig9, ax9 = plt.subplots(figsize=(10, 6))
    mask_oiii = (wave_bins1_single >= 4940) & (wave_bins1_single <= 5020)
    ax9.plot(wave_bins1_single[mask_oiii], spectrum1_single[mask_oiii], 'C0', 
            linewidth=1.5, label='Single Maxwellian', alpha=0.8)
    if spectrum1_double is not None:
        ax9.plot(wave_bins1_double[mask_oiii], spectrum1_double[mask_oiii], 'C1', 
                linewidth=1.5, label='Double Maxwellian', alpha=0.8)
    ax9.set_xlabel('Wavelength [Å]', fontsize=12)
    ax9.set_ylabel('Brightness [Rayleighs/Å]', fontsize=12)
    ax9.set_title('[O III] 4959/5007 Å Nebular Doublet Region', fontsize=14)
    ax9.legend(loc='upper right', fontsize=11)
    ax9.grid(True, alpha=0.3)
    ax9.set_xlim(4940, 5020)
    ax9.set_ylim(bottom=0)
    
    # Add line identification
    ax9.axvline(4958.9, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    ax9.axvline(5006.8, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    ax9.text(4958.9, ax9.get_ylim()[1]*0.9, '[O III]\n4959 Å', 
             ha='center', va='top', fontsize=9)
    ax9.text(5006.8, ax9.get_ylim()[1]*0.9, '[O III]\n5007 Å', 
             ha='center', va='top', fontsize=9)
    
    plt.tight_layout()
    plt.savefig('optical_oiii_doublet_region.png', dpi=300)
    print("Saved: optical_oiii_doublet_region.png")
    
    # ========================================================================
    # STRONGEST EMISSION LINES
    # ========================================================================
    
    print("\n" + "="*70)
    print("Strongest Optical Emission Lines (Equatorial, Single Maxwellian)")
    print("-"*70)
    
    if len(lines1_single) > 0:
        lines_sorted = sorted(lines1_single, key=lambda x: x[1], reverse=True)[:20]
        
        ion_display = {
            'SP': 'S II',
            'S2P': 'S III',
            'S3P': 'S IV',
            'S4P': 'S V',
            'OP': 'O II',
            'O2P': 'O III'
        }
        
        print(f"{'Ion':<8} {'Wavelength':<12} {'Brightness':<12}")
        print("-" * 35)
        for wave, brightness, ion in lines_sorted:
            ion_label = ion_display.get(ion, ion)
            print(f"{ion_label:<8} {wave:>8.2f} Å   {brightness:>8.2f} R")
    
    # ========================================================================
    # KEY DIAGNOSTIC LINE RATIOS
    # ========================================================================
    
    print("\n" + "="*70)
    print("Key Diagnostic Line Ratios")
    print("-"*70)
    
    # Find specific diagnostic lines
    def find_line_brightness(line_list, target_wav, tolerance=5.0):
        """Find brightness of line closest to target wavelength."""
        for wave, brightness, ion in line_list:
            if abs(wave - target_wav) < tolerance:
                return brightness
        return 0.0
    
    # [S II] 6716/6731 ratio (density diagnostic)
    sii_6716 = find_line_brightness(lines1_single, 6716.4, tolerance=3.0)
    sii_6731 = find_line_brightness(lines1_single, 6730.8, tolerance=3.0)
    if sii_6731 > 0:
        sii_ratio = sii_6716 / sii_6731
        print(f"[S II] 6716/6731 ratio (density diagnostic):")
        print(f"  Single Maxwellian: {sii_ratio:.3f}")
        print(f"  (Low density limit ~1.5, high density limit ~0.44)")
    
    # [O III] temperature diagnostic (if lines are present)
    oiii_5007 = find_line_brightness(lines1_single, 5006.8, tolerance=3.0)
    oiii_4959 = find_line_brightness(lines1_single, 4958.9, tolerance=3.0)
    if oiii_4959 > 0:
        oiii_ratio = oiii_5007 / oiii_4959
        print(f"\n[O III] 5007/4959 ratio (should be ~3):")
        print(f"  Single Maxwellian: {oiii_ratio:.3f}")
    
    # ========================================================================
    # DIAGNOSTIC INFORMATION
    # ========================================================================
    
    print("\n" + "="*70)
    print("Diagnostic Information: Plasma Parameters Along Ray (Equatorial)")
    print("-"*70)
    
    # Trace ray and sample plasma parameters
    s_test, pos_test = raytracer.trace_ray(slit_pos_vec1, norm_vec1, ds=1.0)
    
    if len(s_test) > 0:
        print("\nSampling plasma parameters at key positions:")
        print(f"{'Position (R_J)':<30} {'ρ (R_J)':<10} {'ne (cm⁻³)':<12} {'Te (eV)':<10}")
        print("-" * 65)
        
        # Sample at key points
        if len(s_test) > 20:
            sample_indices = [10, len(s_test)//4, len(s_test)//2, 3*len(s_test)//4, len(s_test)-10]
        else:
            sample_indices = [0, len(s_test)//2, len(s_test)-1]
        
        for idx in sample_indices:
            if idx < len(pos_test):
                pos = pos_test[idx]
                ne = raytracer.interp_nec(pos.reshape(1, -1))[0]
                Te = raytracer.interp_Tec(pos.reshape(1, -1))[0]
                rho = np.sqrt(pos[0]**2 + pos[1]**2)
                
                pos_str = f"[{pos[0]:6.1f}, {pos[1]:6.1f}, {pos[2]:6.1f}]"
                print(f"{pos_str:<30} {rho:>8.2f}   {ne:>10.2f}   {Te:>8.2f}")
    
    # ========================================================================
    # PHYSICAL INTERPRETATION
    # ========================================================================
    
    print("\n" + "="*70)
    print("Physical Interpretation - Optical Observations")
    print("-"*70)
    print("The plasma torus optical emission shows:")
    print("- Peak emission near 6 R_J (Io's orbital radius)")
    print("- Scale height of ~0.5-1 R_J")
    print("- Dominant optical emission from [S II] and [O II] forbidden lines")
    print("- [S II] 6716/6731 doublet ratio sensitive to electron density")
    print("- [O III] 5007/4959 ratio fixed by atomic physics (~3:1)")
    print("- [S III] 6312 Å line useful for temperature diagnostics")
    print("- Temperature ~5-10 eV in the cold torus")
    print("- Electron density ~100-2000 cm⁻³ at peak")
    if raytracer.double_maxwellian_loaded:
        print(f"- Hot electrons enhance high-excitation optical lines")
        print(f"- Overall brightness enhancement: {enhancement1:.1%}")
    
    # ========================================================================
    # TIMING AND SUMMARY
    # ========================================================================
    
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    
    print("\n" + "="*70)
    print("Optical Emission Calculation Complete")
    print("="*70)
    print(f"\nExecution time: {elapsed_time:.2f} seconds")
    print("\nOutput files generated:")
    print("  - optical_spectrum_single_equatorial.png")
    print("  - optical_spectrum_single_off_equator.png")
    print("  - optical_spectrum_single_comparison.png")
    if raytracer.double_maxwellian_loaded:
        print("  - optical_spectrum_double_equatorial.png")
        print("  - optical_spectrum_double_off_equator.png")
        print("  - optical_spectrum_single_vs_double.png")
    print("  - optical_sii_doublet_region.png")
    print("  - optical_red_line_region.png")
    print("  - optical_oiii_doublet_region.png")
    
    print("\nCalculation methodology:")
    print("  - Ray tracing through 3D plasma model with trilinear interpolation")
    print("  - Vectorized bilinear interpolation of single Maxwellian emission rates (2D)")
    if raytracer.double_maxwellian_loaded:
        print("  - Vectorized quadrilinear interpolation of double Maxwellian emission rates (4D)")
    print("  - Per-ion photon emission rates × ion density × 1e-6 → Rayleighs")
    print("  - Simpson's rule for line-of-sight integration")
    print("  - ERF-based Gaussian convolution for instrument response")
    print(f"  - Integration step size: ds = {ds} R_J")
    print(f"  - Instrumental FWHM: {fwhm} Å (R~3000 ground-based)")
    print(f"  - Wavelength range: {wavelength_range[0]}-{wavelength_range[1]} Å")
    
    # Show all plots
    plt.show()
    
    # Return results dictionary for further analysis
    results = {
        'wave_bins_single': wave_bins1_single,
        'spectrum_single_equatorial': spectrum1_single,
        'spectrum_single_off_equator': spectrum2_single,
        'lines_single_equatorial': lines1_single,
        'lines_single_off_equator': lines2_single,
        'total_brightness_single_equatorial': total_brightness1_single,
        'total_brightness_single_off_equator': total_brightness2_single,
        'raytracer': raytracer,
        'wavelength_range': wavelength_range,
        'bin_width': bin_width,
        'fwhm': fwhm
    }
    
    if raytracer.double_maxwellian_loaded and spectrum1_double is not None:
        results.update({
            'wave_bins_double': wave_bins1_double,
            'spectrum_double_equatorial': spectrum1_double,
            'spectrum_double_off_equator': spectrum2_double,
            'lines_double_equatorial': lines1_double,
            'lines_double_off_equator': lines2_double,
            'total_brightness_double_equatorial': total_brightness1_double,
            'total_brightness_double_off_equator': total_brightness2_double,
            'enhancement_equatorial': enhancement1,
            'enhancement_off_equator': enhancement2
        })
    
    return results


if __name__ == "__main__":
    # Run the optical emission calculation
    print("="*70)
    print("OPTICAL EMISSION MODEL FOR IO PLASMA TORUS")
    print("Ground-Based Telescope Spectroscopy")
    print("Line-of-Sight Integration with CHIANTI 11.0.2")
    print("="*70)
    print()
    
    results = main()