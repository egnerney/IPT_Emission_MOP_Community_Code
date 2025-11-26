#!/usr/bin/env python3
"""
basic_example_uv_integration_emission_model_use_tables_aligned.py
==========================================================
Example UV Emission Calculation Using Line-of-Sight Integration

This script demonstrates how to use the Jovian UV emission raytracer
to calculate synthetic UV spectra through the Io Plasma Torus using
proper line-of-sight integration with species-by-species interpolation
from CHIANTI emission tables.

The example calculates UV emission spectra for test cases with both
single and double Maxwellian electron distributions:
1. Single Maxwellian: equatorial and off-equator lines of sight
2. Double Maxwellian: same geometries with hot electron population

All lines of sight start at x=6 R_J and look through the plasma torus
in the +y direction, sampling the peak emission region near Io's orbit.

Physical Model:
- Trilinear interpolation of 3D plasma model (densities, temperatures)
- Bilinear interpolation of single Maxwellian volume emission rates (2D tables)
- Quadlinear interpolation of double Maxwellian volume emission rates (4D tables)
- Species-by-species calculation with proper handling of grid boundaries
- Simpson's rule integration over line of sight
- Gaussian instrumental response convolution

Directory Structure Expected:
- IPT_Emission_MOP_Community_Code/
  - LOS_Integration/           (this script and module)
  - Emiss_tables/              (CHIANTI emission tables)
  - 3D_Torus_Model/            (3D plasma model in HDF5 format)

Output:
- UV emission spectra plots for single and double Maxwellian
- Comparison plots showing hot electron enhancement
- List of strongest emission lines
- Diagnostic plasma parameters

Author: Edward (Eddie) G. Nerney
Institution: Laboratory for Atmospheric and Space Physics, University of Colorado Boulder
Date: November 2025
License: Open source for academic and research use
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import sys
import os
from pathlib import Path

# Add the parent directory to path to import the raytracer module
sys.path.insert(0, str(Path(__file__).parent))

# Import the raytracer module
from IPT_emiss_MOP_community_code import JovianUVEmissionRaytracer


def main():
    """
    Main function demonstrating raytracer usage with test cases.
    
    This function:
    1. Initializes the raytracer with plasma model and emission tables
    2. Calculates UV spectra for single Maxwellian distribution
    3. Calculates UV spectra for double Maxwellian distribution
    4. Creates comparison plots showing hot electron enhancement
    5. Prints diagnostic information about plasma conditions
    
    The plasma model and emission tables are loaded from default locations
    based on the expected directory structure. Files are located using
    cross-platform path handling (works on Windows, Mac, Linux).
    """
    print("="*70)
    print("Jovian UV Emission Raytracer - Line-of-Sight Integration")
    print("="*70)
    print()
    print("This example demonstrates UV emission calculations through")
    print("the Io Plasma Torus using proper line-of-sight integration")
    print("with species-by-species interpolation from CHIANTI tables.")
    print()
    
    # Initialize raytracer with default file locations
    print("Loading plasma model and emission tables...")
    try:
        raytracer = JovianUVEmissionRaytracer()
    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("\nPlease ensure the following directory structure exists:")
        print("  IPT_Emission_MOP_Community_Code/")
        print("    - LOS_Integration/ (containing this script)")
        print("    - Emiss_tables/CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5")
        print("    - Emiss_tables/CHIANTI_11.0.2_emiss_tables_double_maxwellian_24x10x24x12.h5")
        print("    - 3D_Torus_Model/jovian_plasma_interpolated_381x381x231.h5")
        return
    
    # ========================================================================
    # SINGLE MAXWELLIAN TEST CASES
    # ========================================================================
    
    # Test case 1: Equatorial line of sight (single Maxwellian)
    print("\n" + "="*70)
    print("Test Case 1: Equatorial Line of Sight (Single Maxwellian)")
    print("-"*70)
    
    # Define observation geometry
    #slit_pos_vec1 = np.array([6.0, -20.0, 0.0])  # Start at x=6 R_J, equator
    #norm_vec1 = np.array([0.0, 1.0, 0.0])        # Look in +y direction
    # Above rotated about z axis by ~+72 degrees 
    slit_pos_vec1 = np.array([20.875, -0.474, 0.0]) # pointing at aligned centrifugal for peak emission about at 6 R_J warm torus
    norm_vec1 = np.array([-0.951, 0.309190879555009, 0.0])       # pointing at aligned centrifugal for peak emission about at 6 R_J warm torus
    
    print(f"Starting position: ({slit_pos_vec1[0]:.1f}, {slit_pos_vec1[1]:.1f}, {slit_pos_vec1[2]:.1f}) R_J")
    print(f"Direction vector: ({norm_vec1[0]:.1f}, {norm_vec1[1]:.1f}, {norm_vec1[2]:.1f})")
    print()
    
    # Calculate UV spectrum (single Maxwellian)
    wave_bins1_single, spectrum1_single, lines1_single = raytracer.calculate_spectrum_single(
        slit_pos_vec1, norm_vec1,
        wavelength_range=(550, 2100),
        bin_width=1.0,
        fwhm=6.0,
        ds=0.01
    )
    
    # Calculate total brightness
    total_brightness1_single = simpson(spectrum1_single, x=wave_bins1_single)
    peak_idx1_single = spectrum1_single.argmax()
    
    print(f"\nSingle Maxwellian Results:")
    print(f"  Total brightness: {total_brightness1_single:.1f} Rayleighs")
    print(f"  Peak brightness: {spectrum1_single.max():.2f} R/Å at {wave_bins1_single[peak_idx1_single]:.1f} Å")
    
    # Test case 2: Off-equator line of sight (single Maxwellian)
    print("\n" + "="*70)
    print("Test Case 2: Off-Equator Line of Sight (Single Maxwellian)")
    print("-"*70)
    
    #slit_pos_vec2 = np.array([6.0, -20.0, 0.5])  # 0.5 R_J above equator
    #norm_vec2 = np.array([0.0, 1.0, 0.0])
    # Above rotated about z axis by ~+72 degrees 
    slit_pos_vec2 = np.array([20.875, -0.474, 0.5]) # pointing at + 0.5 (half scale height) above aligned centrifugal for peak emission about at 6 R_J warm torus
    norm_vec2 = np.array([-0.951, 0.309190879555009, 0.0])       # pointing at +0.5 above (half scale height) aligned centrifugal for peak emission about at 6 R_J warm torus
    
    
    print(f"Starting position: ({slit_pos_vec2[0]:.1f}, {slit_pos_vec2[1]:.1f}, {slit_pos_vec2[2]:.1f}) R_J")
    print(f"Direction vector: ({norm_vec2[0]:.1f}, {norm_vec2[1]:.1f}, {norm_vec2[2]:.1f})")
    print()
    
    wave_bins2_single, spectrum2_single, lines2_single = raytracer.calculate_spectrum_single(
        slit_pos_vec2, norm_vec2,
        wavelength_range=(550, 2100),
        bin_width=1.0,
        fwhm=6.0,
        ds=0.01
    )
    
    total_brightness2_single = simpson(spectrum2_single, x=wave_bins2_single)
    peak_idx2_single = spectrum2_single.argmax()
    
    print(f"\nSingle Maxwellian Results:")
    print(f"  Total brightness: {total_brightness2_single:.1f} Rayleighs")
    print(f"  Peak brightness: {spectrum2_single.max():.2f} R/Å at {wave_bins2_single[peak_idx2_single]:.1f} Å")
    
    # ========================================================================
    # DOUBLE MAXWELLIAN TEST CASES
    # ========================================================================
    
    if raytracer.double_maxwellian_loaded:
        
        # Test case 3: Equatorial line of sight (double Maxwellian)
        print("\n" + "="*70)
        print(f"Test Case 3: Equatorial LOS (Double Maxwellian)")
        print("-"*70)
        
        wave_bins1_double, spectrum1_double, lines1_double = raytracer.calculate_spectrum_double(
            slit_pos_vec1, norm_vec1,
            wavelength_range=(550, 2100),
            bin_width=1.0,
            fwhm=6.0,
            ds=0.01
        )
        
        total_brightness1_double = simpson(spectrum1_double, x=wave_bins1_double)
        peak_idx1_double = spectrum1_double.argmax()
        
        print(f"\nDouble Maxwellian Results:")
        print(f"  Total brightness: {total_brightness1_double:.1f} Rayleighs")
        print(f"  Peak brightness: {spectrum1_double.max():.2f} R/Å at {wave_bins1_double[peak_idx1_double]:.1f} Å")
        print(f"  Enhancement factor: {total_brightness1_double/total_brightness1_single:.3f}")
        
        # Test case 4: Off-equator line of sight (double Maxwellian)
        print("\n" + "="*70)
        print(f"Test Case 4: Off-Equator LOS (Double Maxwellian)")
        print("-"*70)
        
        wave_bins2_double, spectrum2_double, lines2_double = raytracer.calculate_spectrum_double(
            slit_pos_vec2, norm_vec2,
            wavelength_range=(550, 2100),
            bin_width=1.0,
            fwhm=6.0,
            ds=0.01
        )
        
        total_brightness2_double = simpson(spectrum2_double, x=wave_bins2_double)
        peak_idx2_double = spectrum2_double.argmax()
        
        print(f"\nDouble Maxwellian Results:")
        print(f"  Total brightness: {total_brightness2_double:.1f} Rayleighs")
        print(f"  Peak brightness: {spectrum2_double.max():.2f} R/Å at {wave_bins2_double[peak_idx2_double]:.1f} Å")
        print(f"  Enhancement factor: {total_brightness2_double/total_brightness2_single:.3f}")
    
    # ========================================================================
    # SINGLE MAXWELLIAN PLOTS
    # ========================================================================
    
    # Plot single Maxwellian: equatorial
    fig1 = raytracer.plot_spectrum(
        wave_bins1_single, spectrum1_single, lines1_single,
        title="Single Maxwellian: Equatorial LOS (x=6 R_J, z=0)"
    )
    plt.savefig('spectrum_single_equatorial.png', dpi=150)
    print("\nSaved: spectrum_single_equatorial.png")
    
    # Plot single Maxwellian: off-equator
    fig2 = raytracer.plot_spectrum(
        wave_bins2_single, spectrum2_single, lines2_single,
        title="Single Maxwellian: Off-Equator LOS (x=6 R_J, z=0.5 R_J)"
    )
    plt.savefig('spectrum_single_off_equator.png', dpi=150)
    print("Saved: spectrum_single_off_equator.png")
    
    # Comparison plot: single Maxwellian equatorial vs off-equator
    fig3, ax3 = plt.subplots(figsize=(12, 6))
    ax3.plot(wave_bins1_single, spectrum1_single, 'b-', label='Equatorial (z=0)', linewidth=1.0)
    ax3.plot(wave_bins2_single, spectrum2_single, 'r-', label='Off-equator (z=0.5 R_J)',
             linewidth=1.0, alpha=0.7)
    ax3.set_xlabel('Wavelength (Å)', fontsize=12)
    ax3.set_ylabel('Brightness (Rayleighs/Å)', fontsize=12)
    ax3.set_title('Single Maxwellian: Equatorial vs Off-Equator UV Spectra', fontsize=14)
    ax3.grid(True, alpha=0.3)
    ax3.legend(loc='upper right')
    ax3.set_xlim(550, 2100)
    ax3.set_ylim(bottom=0)
    
    ax3.text(0.02, 0.98, f'Equatorial: {total_brightness1_single:.1f} R\n'
                        f'Off-equator: {total_brightness2_single:.1f} R',
             transform=ax3.transAxes, va='top', ha='left',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig('spectrum_single_comparison.png', dpi=150)
    print("Saved: spectrum_single_comparison.png")
    
    # ========================================================================
    # DOUBLE MAXWELLIAN PLOTS (if available)
    # ========================================================================
    
    if raytracer.double_maxwellian_loaded:
        # Plot double Maxwellian: equatorial
        fig4 = raytracer.plot_spectrum(
            wave_bins1_double, spectrum1_double, lines1_double,
            title=f"Double Maxwellian: Equatorial LOS"
        )
        plt.savefig('spectrum_double_equatorial.png', dpi=150)
        print("Saved: spectrum_double_equatorial.png")
        
        # Plot double Maxwellian: off-equator
        fig5 = raytracer.plot_spectrum(
            wave_bins2_double, spectrum2_double, lines2_double,
            title=f"Double Maxwellian: Off-Equator LOS"
        )
        plt.savefig('spectrum_double_off_equator.png', dpi=150)
        print("Saved: spectrum_double_off_equator.png")
        
        # Comparison plot: single vs double Maxwellian (equatorial)
        fig6, ax6 = plt.subplots(figsize=(12, 6))
        ax6.plot(wave_bins1_single, spectrum1_single, 'b-', label='Single Maxwellian', linewidth=1.0)
        ax6.plot(wave_bins1_double, spectrum1_double, 'r-', label='Double Maxwellian',
                 linewidth=1.0, alpha=0.7)
        ax6.set_xlabel('Wavelength (Å)', fontsize=12)
        ax6.set_ylabel('Brightness (Rayleighs/Å)', fontsize=12)
        ax6.set_title('Equatorial LOS: Single vs Double Maxwellian', fontsize=14)
        ax6.grid(True, alpha=0.3)
        ax6.legend(loc='upper right')
        ax6.set_xlim(550, 2100)
        ax6.set_ylim(bottom=0)
        
        enhancement = total_brightness1_double / total_brightness1_single
        ax6.text(0.02, 0.98, f'Single: {total_brightness1_single:.1f} R\n'
                            f'Double: {total_brightness1_double:.1f} R\n'
                            f'Enhancement: {enhancement:.3f}×',
                 transform=ax6.transAxes, va='top', ha='left',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        plt.savefig('spectrum_single_vs_double.png', dpi=150)
        print("Saved: spectrum_single_vs_double.png")
    
    # ========================================================================
    # STRONGEST EMISSION LINES
    # ========================================================================
    
    print("\n" + "="*70)
    print("Strongest UV Emission Lines (Equatorial, Single Maxwellian)")
    print("-"*70)
    
    if len(lines1_single) > 0:
        lines_sorted = sorted(lines1_single, key=lambda x: x[1], reverse=True)[:15]
        
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
    print("Physical Interpretation")
    print("-"*70)
    print("The plasma torus shows:")
    print("- Peak emission near 6 R_J (Io's orbital radius)")
    print("- Scale height of ~0.5-1 R_J")
    print("- Dominant emission from S III and O II ions")
    print("- Temperature ~5-10 eV in the cold torus")
    print("- Electron density ~100-2000 cm⁻³ at peak")
    if raytracer.double_maxwellian_loaded:
        print(f"- Hot electrons enhance high-excitation lines")
        print(f"- Overall brightness enhancement: {enhancement:.1%}")
    
    # Show all plots
    plt.show()
    
    # ========================================================================
    # SUMMARY
    # ========================================================================
    
    print("\n" + "="*70)
    print("UV Emission Calculation Complete")
    print("="*70)
    print("\nOutput files generated:")
    print("  - spectrum_single_equatorial.png")
    print("  - spectrum_single_off_equator.png")
    print("  - spectrum_single_comparison.png")
    if raytracer.double_maxwellian_loaded:
        print("  - spectrum_double_equatorial.png")
        print("  - spectrum_double_off_equator.png")
        print("  - spectrum_single_vs_double.png")
    
    print("\nAll calculations performed using:")
    print("  - Trilinear interpolation of 3D plasma model")
    print("  - Bilinear interpolation of single Maxwellian volume emission rates (2D)")
    if raytracer.double_maxwellian_loaded:
        print("  - Quadlinear interpolation of double Maxwellian volume emission rates (4D)")
    print("  - Simpson's rule for line-of-sight integration")
    print("  - Gaussian instrumental response (FWHM=6.0 Å)")
    print("  - ds = 0.01 R_J integration step size")
    print("  - Species-by-species calculation with proper grid handling")
    print("  - No extrapolation outside table bounds (set to 0)")


if __name__ == "__main__":
    main()