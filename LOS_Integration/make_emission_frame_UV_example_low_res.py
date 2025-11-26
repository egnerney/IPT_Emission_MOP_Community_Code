#!/usr/bin/env python3
"""
make_emission_frame_UV_example_low_res.py
==========================================
Generate UV Emission Contour Maps - Low Resolution Version

This script creates contour plots of UV emission from the Jovian plasma torus
as a function of viewing position, using uniform grid spacing for faster 
computation suitable for quick exploration or testing.

The script generates three contour plots:
1. Total UV emission (550-2100 Å)
2. S++ emission around 680 Å (S III line)
3. O+ emission around 833 Å (O II line)

Grid Resolution:
- X position: -10 to +10 R_J in 0.1 R_J steps (201 points)
- Z position: -2.5 to +2.5 R_J in 0.1 R_J steps (51 points)
- Total grid: 201 × 51 = 10,251 calculations

Line Integration:
- Integrates emission within ±3σ of line center
- σ = FWHM / 2.35482 for Gaussian profile
- Uses Simpson's rule for wavelength integration

Output:
- Three contour plots with viridis colormap (100 levels)
- Saved as PNG files

Author: Edward (Eddie) G. Nerney Nov 2025

License: MIT
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import sys
from pathlib import Path
import time

# Add the parent directory to path to import the raytracer module
sys.path.insert(0, str(Path(__file__).parent))

# Import the raytracer module
from IPT_emiss_MOP_community_code import JovianUVEmissionRaytracer


def integrate_line_emission(raytracer, slit_pos_vec, norm_vec, 
                           center_wavelength, fwhm=6.0):
    """
    Integrate emission around a specific spectral line.
    
    Integrates within ±3σ of the line center where σ = FWHM/2.35482.
    
    Parameters:
    -----------
    raytracer : JovianUVEmissionRaytracer
        The raytracer object
    slit_pos_vec : array [x0, y0, z0]
        Starting position in R_J
    norm_vec : array [nx, ny, nz]
        Normalized direction vector
    center_wavelength : float
        Central wavelength of the line in Angstroms
    fwhm : float
        Instrumental FWHM in Angstroms (default 6.0 Å)
    
    Returns:
    --------
    integrated_emission : float
        Integrated emission in Rayleighs
    """
    # Calculate sigma from FWHM
    # For Gaussian: FWHM = 2.35482 * sigma
    sigma = fwhm / 2.35482
    
    # Integration range: ±3σ around line center
    wave_min = center_wavelength - 3 * sigma
    wave_max = center_wavelength + 3 * sigma
    
    # Create wavelength grid for integration (0.1 Å sampling)
    wave_grid = np.arange(wave_min, wave_max + 0.1, 0.1)
    
    # Calculate spectrum over this range
    wave_bins, spectrum, _ = raytracer.calculate_spectrum(
        slit_pos_vec, norm_vec,
        wavelength_range=(wave_min - 1, wave_max + 1),
        bin_width=0.1,
        fwhm=fwhm, ds=0.1
    )
    
    # Find the portion within our integration range
    mask = (wave_bins >= wave_min) & (wave_bins <= wave_max)
    
    if np.any(mask):
        # Integrate using Simpson's rule
        integrated = simpson(spectrum[mask], x=wave_bins[mask])
        return integrated
    else:
        return 0.0


def main():
    """
    Main function to generate UV emission contour maps.
    
    Creates contour plots showing emission as a function of viewing position
    for three cases: total UV, S++ 680 Å, and O+ 833 Å.
    """
    print("="*70)
    print("UV Emission Contour Maps - Low Resolution")
    print("="*70)
    print()
    
    # Initialize raytracer
    print("Loading plasma model and emission tables...")
    try:
        raytracer = JovianUVEmissionRaytracer()
    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("\nPlease ensure the data files are in the correct directories.")
        return
    
    print("\nSetting up calculation grid...")
    
    # Define position grids - uniform spacing for low resolution
    # X grid: -10 to +10 R_J in 0.1 R_J steps
    x_grid = np.arange(-10.0, 10.1, 0.1)
    # Z grid: -2.5 to +2.5 R_J in 0.1 R_J steps  
    z_grid = np.arange(-2.5, 2.6, 0.1)
    
    nx = len(x_grid)
    nz = len(z_grid)
    print(f"Grid dimensions: {nx} x {nz} = {nx*nz} total calculations")
    print(f"X range: {x_grid[0]:.1f} to {x_grid[-1]:.1f} R_J ({nx} points)")
    print(f"Z range: {z_grid[0]:.1f} to {z_grid[-1]:.1f} R_J ({nz} points)")
    
    # Initialize arrays for storing results
    total_emission = np.zeros((nz, nx))
    s3_emission = np.zeros((nz, nx))    # S++ at 680 Å
    o2_emission = np.zeros((nz, nx))    # O+ at 833 Å
    
    # Fixed parameters
    y_start = -20.0  # Starting y position (R_J)
    norm_vec = np.array([0.0, 1.0, 0.0])  # Always looking in +y direction
    
    # Wavelength centers for specific lines
    s3_wavelength = 680.0   # S III (S++) line
    o2_wavelength = 833.0   # O II (O+) line
    
    print(f"\nFixed parameters:")
    print(f"  Starting y position: {y_start:.1f} R_J")
    print(f"  Viewing direction: +y (through torus)")
    print(f"  S++ line center: {s3_wavelength:.1f} Å")
    print(f"  O+ line center: {o2_wavelength:.1f} Å")
    print(f"  Integration range: ±3σ around line centers")
    
    # Calculate emissions for each grid point
    print("\nCalculating emissions...")
    print("Progress: ", end="", flush=True)
    
    start_time = time.time()
    total_calculations = nx * nz
    calculation_count = 0
    
    for i, x0 in enumerate(x_grid):
        for j, z0 in enumerate(z_grid):
            # Define starting position
            slit_pos_vec = np.array([x0, y_start, z0])
            
            # Calculate total UV emission (550-2100 Å)
            wave_bins, spectrum, _ = raytracer.calculate_spectrum(
                slit_pos_vec, norm_vec,
                wavelength_range=(550, 2100),
                bin_width=1.0,
                fwhm=6.0, ds=0.1
            )
            total_emission[j, i] = simpson(spectrum, x=wave_bins)
            
            # Calculate S++ emission around 680 Å
            s3_emission[j, i] = integrate_line_emission(
                raytracer, slit_pos_vec, norm_vec, 
                s3_wavelength, fwhm=6.0
            )
            
            # Calculate O+ emission around 833 Å
            o2_emission[j, i] = integrate_line_emission(
                raytracer, slit_pos_vec, norm_vec,
                o2_wavelength, fwhm=6.0
            )
            
            # Progress indicator
            calculation_count += 1
            if calculation_count % 100 == 0:
                percent = (calculation_count / total_calculations) * 100
                print(f"{percent:.0f}%...", end="", flush=True)
    
    elapsed_time = time.time() - start_time
    print(f" Done!")
    print(f"Calculation time: {elapsed_time:.1f} seconds")
    print(f"Average time per point: {elapsed_time/total_calculations:.3f} seconds")
    
    # Print emission statistics
    print("\n" + "="*70)
    print("Emission Statistics")
    print("-"*70)
    print(f"Total UV emission (550-2100 Å):")
    print(f"  Min: {np.min(total_emission):.1f} R")
    print(f"  Max: {np.max(total_emission):.1f} R")
    print(f"  Mean: {np.mean(total_emission):.1f} R")
    
    print(f"\nS++ emission (680 Å ± 3σ):")
    print(f"  Min: {np.min(s3_emission):.2f} R")
    print(f"  Max: {np.max(s3_emission):.2f} R")
    print(f"  Mean: {np.mean(s3_emission):.2f} R")
    
    print(f"\nO+ emission (833 Å ± 3σ):")
    print(f"  Min: {np.min(o2_emission):.2f} R")
    print(f"  Max: {np.max(o2_emission):.2f} R")
    print(f"  Mean: {np.mean(o2_emission):.2f} R")
    
    # Create contour plots
    print("\n" + "="*70)
    print("Creating Contour Plots")
    print("-"*70)
    
    # Create figure with three subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # Plot 1: Total UV emission
    ax1 = axes[0]
    X, Z = np.meshgrid(x_grid, z_grid)
    cs1 = ax1.contourf(X, Z, total_emission, levels=100, cmap='viridis')
    ax1.contour(X, Z, total_emission, levels=10, colors='black', alpha=0.3, linewidths=0.5)
    ax1.set_xlabel('X Position (R_J)', fontsize=12)
    ax1.set_ylabel('Z Position (R_J)', fontsize=12)
    ax1.set_title('Total UV Emission (550-2100 Å)', fontsize=14)
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)
    cbar1 = plt.colorbar(cs1, ax=ax1, label='Brightness (R)')
    
    
    # Plot 2: S++ emission
    ax2 = axes[1]
    cs2 = ax2.contourf(X, Z, s3_emission, levels=100, cmap='viridis')
    ax2.contour(X, Z, s3_emission, levels=10, colors='black', alpha=0.3, linewidths=0.5)
    ax2.set_xlabel('X Position (R_J)', fontsize=12)
    ax2.set_ylabel('Z Position (R_J)', fontsize=12)
    ax2.set_title('S$^{++}$ Emission (680 Å ± 3σ)', fontsize=14)
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)
    cbar2 = plt.colorbar(cs2, ax=ax2, label='Brightness (R)')
    
    
    # Plot 3: O+ emission
    ax3 = axes[2]
    cs3 = ax3.contourf(X, Z, o2_emission, levels=100, cmap='viridis')
    ax3.contour(X, Z, o2_emission, levels=10, colors='black', alpha=0.3, linewidths=0.5)
    ax3.set_xlabel('X Position (R_J)', fontsize=12)
    ax3.set_ylabel('Z Position (R_J)', fontsize=12)
    ax3.set_title('O$^{+}$ Emission (833 Å ± 3σ)', fontsize=14)
    ax3.set_aspect('equal')
    ax3.grid(True, alpha=0.3)
    cbar3 = plt.colorbar(cs3, ax=ax3, label='Brightness (R)')
    

    
    plt.suptitle('UV Emission Maps: View Along +Y Direction (Low Resolution)', 
                 fontsize=16, y=1.02)
    plt.tight_layout()
    
    # Save the figure
    output_file = 'emission_contours_low_res.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nSaved contour plots to: {output_file}")
    
    # Create individual plots for better visibility
    print("\nCreating individual contour plots...")
    
    # Individual plot for total emission
    fig1, ax = plt.subplots(figsize=(8, 6))
    cs = ax.contourf(X, Z, total_emission, levels=100, cmap='viridis')
    ax.contour(X, Z, total_emission, levels=10, colors='black', alpha=0.3, linewidths=0.5)
    ax.set_xlabel('X Position (R_J)', fontsize=12)
    ax.set_ylabel('Z Position (R_J)', fontsize=12)
    ax.set_title('Total UV Emission (550-2100 Å) - Low Resolution', fontsize=14)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    plt.colorbar(cs, ax=ax, label='Brightness (R)')

    plt.tight_layout()
    plt.savefig('total_emission_low_res.png', dpi=150)
    print("  Saved: total_emission_low_res.png")
    
    # Individual plot for S++ emission
    fig2, ax = plt.subplots(figsize=(8, 6))
    cs = ax.contourf(X, Z, s3_emission, levels=100, cmap='viridis')
    ax.contour(X, Z, s3_emission, levels=10, colors='black', alpha=0.3, linewidths=0.5)
    ax.set_xlabel('X Position (R_J)', fontsize=12)
    ax.set_ylabel('Z Position (R_J)', fontsize=12)
    ax.set_title('S$^{++}$ Emission (680 Å ± 3σ) - Low Resolution', fontsize=14)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    plt.colorbar(cs, ax=ax, label='Brightness (R)')

    plt.tight_layout()
    plt.savefig('s3_emission_low_res.png', dpi=150)
    print("  Saved: s3_emission_low_res.png")
    
    # Individual plot for O+ emission
    fig3, ax = plt.subplots(figsize=(8, 6))
    cs = ax.contourf(X, Z, o2_emission, levels=100, cmap='viridis')
    ax.contour(X, Z, o2_emission, levels=10, colors='black', alpha=0.3, linewidths=0.5)
    ax.set_xlabel('X Position (R_J)', fontsize=12)
    ax.set_ylabel('Z Position (R_J)', fontsize=12)
    ax.set_title('O$^{+}$ Emission (833 Å ± 3σ) - Low Resolution', fontsize=14)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    plt.colorbar(cs, ax=ax, label='Brightness (R)')
   
    plt.tight_layout()
    plt.savefig('o2_emission_low_res.png', dpi=150)
    print("  Saved: o2_emission_low_res.png")
    
    # Show all plots
    plt.show()
    
    print("\n" + "="*70)
    print("Emission Frame Generation Complete")
    print("="*70)
    print("\nPhysical Interpretation:")
    print("- White dashed circle marks Io's orbital radius (~5.9 R_J)")
    print("- Peak emission occurs near 6 R_J in the equatorial plane")
    print("- Vertical extent shows scale height of ~0.5-1 R_J")
    print("- S++ and O+ show different spatial distributions")
    print("- Emission drops rapidly inside 5 R_J and outside 8 R_J")


if __name__ == "__main__":
    main()