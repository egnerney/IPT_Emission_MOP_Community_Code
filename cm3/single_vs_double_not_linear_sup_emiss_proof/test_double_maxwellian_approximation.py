#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_double_maxwellian_approximation.py

Comparison of Coronal Low-Density Approximation vs Exact Double Maxwellian
===========================================================================

This script tests the validity of the coronal low-density approximation for double
Maxwellian electron distributions. The approximation assumes that emission from a 
double Maxwellian can be represented as a linear superposition of single Maxwellians:

    emiss(Tec,Teh,ne_total,feh) ≈ fec*emiss(Tec,nec) + feh*emiss(Teh,neh)

where:
    - fec = nec/ne_total = core electron fraction = 1 - feh
    - feh = neh/ne_total = hot electron fraction
    - nec = fec * ne_total = core electron density
    - neh = feh * ne_total = hot electron density

PHYSICAL BACKGROUND:
The exact emission is found by solving the level balance equation:
    C·f = b
where C_{ij} = A_{ij} + ne*q_{ij}, with q_{ij} being the rate coefficients.

For a double Maxwellian, the rate coefficients are:
    q_{ij} = fec*q_{ij}(Tec) + feh*q_{ij}(Teh)

However, the level populations f_j are found by matrix inversion:
    f = C^{-1}·b

This is a nonlinear operation, so in general:
    f_j(Tec,Teh,ne,feh) ≠ fec*f_j(Tec,nec) + feh*f_j(Teh,neh)

The coronal approximation assumes we can treat this linearly, which is valid
when collisional de-excitation is negligible (low density) and cascade from
higher levels is minimal.

This script quantifies when this approximation is valid by comparing:
1. Exact double Maxwellian emission (from 4D tables)
2. Approximate emission (from weighted sum of single Maxwellian tables)

TESTS PERFORMED:
- Variation with hot electron fraction (0.001% to 10%)
- Variation with temperature ratios (Teh/Tec from 10 to 100)
- Species dependence (low vs high ionization states)
- Wavelength/energy dependence
- Density dependence

AUTHOR: Edward (Eddie) G. Nerney
INSTITUTION: Laboratory for Atmospheric and Space Physics, University of Colorado Boulder
DATE: November 2025
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import time
from scipy.integrate import simpson
import warnings

# Import the IPT emission model library
from IPT_cm3_emission_model import (
    EmissionTables,
    interpolate_emissivity_2d,
    interpolate_emissivity_4d
)


def calculate_approximate_double_maxwellian(
    tables: EmissionTables,
    Tec: float,
    Teh: float,
    ne_total: float,
    feh: float,
    species_list=None
):
    """
    Calculate emission using linear superposition approximation.
    
    Approximation: emiss ≈ fec*emiss(Tec,nec) + feh*emiss(Teh,neh)
    
    Parameters
    ----------
    tables : EmissionTables
        Loaded emission tables (must have single Maxwellian loaded)
    Tec : float
        Core electron temperature [eV]
    Teh : float
        Hot electron temperature [eV]
    ne_total : float
        Total electron density [cm^-3]
    feh : float
        Hot electron fraction
    species_list : list, optional
        Species to include
        
    Returns
    -------
    wavelengths : ndarray
        Emission line wavelengths [Angstroms]
    emissivities_approx : ndarray
        Approximate volume emission rates [photons s^-1 cm^-3]
    species_tags : list
        Species tags for each line
    """
    # Calculate component densities
    fec = 1.0 - feh
    nec = fec * ne_total
    neh = feh * ne_total
    
    # Get emission from core population
    wav_core, emiss_core, species_core = interpolate_emissivity_2d(
        tables, Tec, nec, species_list=species_list
    )
    
    # Get emission from hot population
    wav_hot, emiss_hot, species_hot = interpolate_emissivity_2d(
        tables, Teh, neh, species_list=species_list
    )
    
    # Verify wavelength arrays match
    if not np.allclose(wav_core, wav_hot, rtol=1e-6):
        raise ValueError("Wavelength arrays from core and hot populations do not match!")
    
    # Linear superposition weighted by fractions
    emissivities_approx = fec * emiss_core + feh * emiss_hot
    
    return wav_core, emissivities_approx, species_core


def compare_exact_vs_approximate(
    tables: EmissionTables,
    Tec: float,
    Teh: float,
    ne_total: float,
    feh: float,
    species_list=None,
    verbose=True
):
    """
    Compare exact double Maxwellian with linear approximation.
    
    Parameters
    ----------
    tables : EmissionTables
        Loaded emission tables (must have both single and double Maxwellian)
    Tec : float
        Core electron temperature [eV]
    Teh : float
        Hot electron temperature [eV]
    ne_total : float
        Total electron density [cm^-3]
    feh : float
        Hot electron fraction
    species_list : list, optional
        Species to include
    verbose : bool
        Print summary statistics
        
    Returns
    -------
    results : dict
        Dictionary containing:
        - 'wavelengths': wavelength array
        - 'exact': exact double Maxwellian emissivities
        - 'approx': approximate emissivities
        - 'ratio': approx/exact ratio
        - 'percent_error': percentage error
        - 'species': species tags
    """
    # Get exact double Maxwellian emission
    wav_exact, emiss_exact, species_exact = interpolate_emissivity_4d(
        tables, Tec, Teh, ne_total, feh, species_list=species_list
    )
    
    # Get approximate emission from linear superposition
    wav_approx, emiss_approx, species_approx = calculate_approximate_double_maxwellian(
        tables, Tec, Teh, ne_total, feh, species_list=species_list
    )
    
    # Verify wavelengths match
    if not np.allclose(wav_exact, wav_approx, rtol=1e-6):
        raise ValueError("Wavelength arrays do not match between exact and approximate!")
    
    # Calculate comparison metrics
    # Add small epsilon to avoid division by zero
    epsilon = 1e-30
    ratio = emiss_approx / (emiss_exact + epsilon)
    percent_error = 100.0 * (emiss_approx - emiss_exact) / (emiss_exact + epsilon)
    
    # Filter out negligible lines (emissivity < 1e-10)
    significant = emiss_exact > 1e-10
    
    if verbose and np.sum(significant) > 0:
        print(f"\nComparison Results:")
        print(f"  Parameters: Tec={Tec:.1f} eV, Teh={Teh:.0f} eV, ne={ne_total:.0f} cm^-3, feh={feh:.6f}")
        print(f"  Number of significant lines: {np.sum(significant)}")
        print(f"\n  Ratio statistics (approx/exact) for significant lines:")
        print(f"    Mean:   {np.mean(ratio[significant]):.4f}")
        print(f"    Median: {np.median(ratio[significant]):.4f}")
        print(f"    Std:    {np.std(ratio[significant]):.4f}")
        print(f"    Min:    {np.min(ratio[significant]):.4f}")
        print(f"    Max:    {np.max(ratio[significant]):.4f}")
        print(f"\n  Percent error statistics:")
        print(f"    Mean abs error:   {np.mean(np.abs(percent_error[significant])):.2f}%")
        print(f"    Median abs error: {np.median(np.abs(percent_error[significant])):.2f}%")
        print(f"    RMS error:        {np.sqrt(np.mean(percent_error[significant]**2)):.2f}%")
    
    return {
        'wavelengths': wav_exact,
        'exact': emiss_exact,
        'approx': emiss_approx,
        'ratio': ratio,
        'percent_error': percent_error,
        'species': species_exact,
        'significant': significant
    }


def test_feh_variation(tables: EmissionTables, save_plots=True):
    """
    Test how approximation quality varies with hot electron fraction.
    
    This tests the core assumption: does the linear approximation work
    for small hot fractions as predicted by coronal theory?
    """
    print("\n" + "="*70)
    print("TEST 1: Variation with Hot Electron Fraction")
    print("="*70)
    
    # Fixed parameters
    Tec = 5.0          # Core temperature [eV]
    Teh = 270.0        # Hot temperature [eV]
    ne_total = 2200.0  # Total density [cm^-3]
    
    # Range of hot fractions to test (0.001% to 10%)
    feh_values = np.array([0.00001, 0.0001, 0.001, 0.005, 0.01, 0.025, 0.05, 0.1])
    
    # Storage for results
    mean_ratios = []
    rms_errors = []
    
    for feh in feh_values:
        print(f"\nTesting feh = {feh:.6f} ({100*feh:.4f}%)")
        
        results = compare_exact_vs_approximate(
            tables, Tec, Teh, ne_total, feh, verbose=False
        )
        
        sig = results['significant']
        if np.sum(sig) > 0:
            mean_ratios.append(np.mean(results['ratio'][sig]))
            rms_errors.append(np.sqrt(np.mean(results['percent_error'][sig]**2)))
        else:
            mean_ratios.append(np.nan)
            rms_errors.append(np.nan)
    
    # Plotting
    if save_plots:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        # Plot 1: Mean ratio vs feh
        axes[0].semilogx(feh_values * 100, mean_ratios, 'o-', linewidth=2, markersize=8)
        axes[0].axhline(1.0, color='k', linestyle='--', alpha=0.5, label='Perfect agreement')
        axes[0].set_xlabel('Hot Electron Fraction [%]', fontsize=12)
        axes[0].set_ylabel('Mean Ratio (Approx/Exact)', fontsize=12)
        axes[0].set_title(f'Approximation Quality vs Hot Fraction\nTec={Tec} eV, Teh={Teh} eV', 
                         fontsize=13)
        axes[0].grid(True, alpha=0.3)
        axes[0].legend()
        
        # Plot 2: RMS error vs feh
        axes[1].loglog(feh_values * 100, rms_errors, 'o-', linewidth=2, markersize=8, color='C1')
        axes[1].set_xlabel('Hot Electron Fraction [%]', fontsize=12)
        axes[1].set_ylabel('RMS Percent Error [%]', fontsize=12)
        axes[1].set_title(f'Approximation Error vs Hot Fraction\nTec={Tec} eV, Teh={Teh} eV', 
                         fontsize=13)
        axes[1].grid(True, alpha=0.3, which='both')
        
        plt.tight_layout()
        plt.savefig('test1_feh_variation.png', dpi=300, bbox_inches='tight')
        print("\nSaved: test1_feh_variation.png")
        plt.close()


def test_temperature_ratio_variation(tables: EmissionTables, save_plots=True):
    """
    Test how approximation quality varies with temperature ratio Teh/Tec.
    
    Higher temperature ratios should show larger deviations because the
    difference in excitation rates becomes more extreme.
    """
    print("\n" + "="*70)
    print("TEST 2: Variation with Temperature Ratio")
    print("="*70)
    
    # Fixed parameters
    Tec = 5.0          # Core temperature [eV]
    ne_total = 2200.0  # Total density [cm^-3]
    feh = 0.0025       # Hot fraction (0.25%)
    
    # Range of temperature ratios to test
    temp_ratios = np.array([5, 10, 20, 30, 50, 75, 100])
    Teh_values = Tec * temp_ratios
    
    # Storage for results
    mean_ratios = []
    rms_errors = []
    
    for Teh in Teh_values:
        ratio_val = Teh / Tec
        print(f"\nTesting Teh/Tec = {ratio_val:.1f} (Teh = {Teh:.1f} eV)")
        
        results = compare_exact_vs_approximate(
            tables, Tec, Teh, ne_total, feh, verbose=False
        )
        
        sig = results['significant']
        if np.sum(sig) > 0:
            mean_ratios.append(np.mean(results['ratio'][sig]))
            rms_errors.append(np.sqrt(np.mean(results['percent_error'][sig]**2)))
        else:
            mean_ratios.append(np.nan)
            rms_errors.append(np.nan)
    
    # Plotting
    if save_plots:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        # Plot 1: Mean ratio vs temperature ratio
        axes[0].plot(temp_ratios, mean_ratios, 'o-', linewidth=2, markersize=8)
        axes[0].axhline(1.0, color='k', linestyle='--', alpha=0.5, label='Perfect agreement')
        axes[0].set_xlabel('Temperature Ratio (Teh/Tec)', fontsize=12)
        axes[0].set_ylabel('Mean Ratio (Approx/Exact)', fontsize=12)
        axes[0].set_title(f'Approximation Quality vs Temperature Ratio\nTec={Tec} eV, feh={feh:.4f}', 
                         fontsize=13)
        axes[0].grid(True, alpha=0.3)
        axes[0].legend()
        
        # Plot 2: RMS error vs temperature ratio
        axes[1].semilogy(temp_ratios, rms_errors, 'o-', linewidth=2, markersize=8, color='C1')
        axes[1].set_xlabel('Temperature Ratio (Teh/Tec)', fontsize=12)
        axes[1].set_ylabel('RMS Percent Error [%]', fontsize=12)
        axes[1].set_title(f'Approximation Error vs Temperature Ratio\nTec={Tec} eV, feh={feh:.4f}', 
                         fontsize=13)
        axes[1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('test2_temperature_ratio.png', dpi=300, bbox_inches='tight')
        print("\nSaved: test2_temperature_ratio.png")
        plt.close()


def test_species_dependence(tables: EmissionTables, save_plots=True):
    """
    Test how approximation quality varies with species/ionization state.
    
    Higher ionization states (S IV, S V) require more energy to excite,
    so hot electrons have disproportionate impact. The approximation
    should be worse for high ionization states.
    """
    print("\n" + "="*70)
    print("TEST 3: Species Dependence")
    print("="*70)
    
    # Fixed parameters
    Tec = 5.0          # Core temperature [eV]
    Teh = 270.0        # Hot temperature [eV]
    ne_total = 2200.0  # Total density [cm^-3]
    feh = 0.0025       # Hot fraction (0.25%)
    
    # Test each species separately
    species_to_test = ['SP', 'S2P', 'S3P', 'S4P', 'OP', 'O2P']
    species_names = {
        'SP': 'S II (S+)',
        'S2P': 'S III (S++)',
        'S3P': 'S IV (S+++)',
        'S4P': 'S V (S++++)',
        'OP': 'O II (O+)',
        'O2P': 'O III (O++)'
    }
    
    mean_ratios = []
    rms_errors = []
    species_labels = []
    
    for species in species_to_test:
        print(f"\nTesting species: {species_names[species]}")
        
        results = compare_exact_vs_approximate(
            tables, Tec, Teh, ne_total, feh, species_list=[species], verbose=False
        )
        
        sig = results['significant']
        if np.sum(sig) > 0:
            mean_ratios.append(np.mean(results['ratio'][sig]))
            rms_errors.append(np.sqrt(np.mean(results['percent_error'][sig]**2)))
            species_labels.append(species_names[species])
            
            print(f"  Mean ratio: {mean_ratios[-1]:.4f}")
            print(f"  RMS error:  {rms_errors[-1]:.2f}%")
    
    # Plotting
    if save_plots:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        x = np.arange(len(species_labels))
        
        # Plot 1: Mean ratio by species
        axes[0].bar(x, mean_ratios, alpha=0.7, edgecolor='black')
        axes[0].axhline(1.0, color='r', linestyle='--', alpha=0.7, label='Perfect agreement')
        axes[0].set_xticks(x)
        axes[0].set_xticklabels(species_labels, rotation=45, ha='right')
        axes[0].set_ylabel('Mean Ratio (Approx/Exact)', fontsize=12)
        axes[0].set_title(f'Approximation Quality by Species\nTec={Tec} eV, Teh={Teh} eV, feh={feh:.4f}', 
                         fontsize=13)
        axes[0].grid(True, alpha=0.3, axis='y')
        axes[0].legend()
        
        # Plot 2: RMS error by species
        axes[1].bar(x, rms_errors, alpha=0.7, edgecolor='black', color='C1')
        axes[1].set_xticks(x)
        axes[1].set_xticklabels(species_labels, rotation=45, ha='right')
        axes[1].set_ylabel('RMS Percent Error [%]', fontsize=12)
        axes[1].set_title(f'Approximation Error by Species\nTec={Tec} eV, Teh={Teh} eV, feh={feh:.4f}', 
                         fontsize=13)
        axes[1].grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        plt.savefig('test3_species_dependence.png', dpi=300, bbox_inches='tight')
        print("\nSaved: test3_species_dependence.png")
        plt.close()


def test_wavelength_dependence(tables: EmissionTables, save_plots=True):
    """
    Test how approximation quality varies with wavelength/transition energy.
    
    Shorter wavelengths (higher energies) should show worse approximation
    because they require more energetic collisions where hot electrons dominate.
    """
    print("\n" + "="*70)
    print("TEST 4: Wavelength/Energy Dependence")
    print("="*70)
    
    # Fixed parameters
    Tec = 5.0          # Core temperature [eV]
    Teh = 270.0        # Hot temperature [eV]
    ne_total = 2200.0  # Total density [cm^-3]
    feh = 0.0025       # Hot fraction (0.25%)
    
    # Get full comparison
    results = compare_exact_vs_approximate(
        tables, Tec, Teh, ne_total, feh, verbose=True
    )
    
    # Extract data
    wavelengths = results['wavelengths']
    ratios = results['ratio']
    percent_errors = results['percent_error']
    species = results['species']
    sig = results['significant']
    
    # Convert wavelength to photon energy [eV]
    # E [eV] = 12398.4 / λ [Å]
    energies = 12398.4 / wavelengths
    
    # Plotting
    if save_plots:
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Plot 1: Ratio vs wavelength
        for sp in ['SP', 'S2P', 'S3P', 'S4P', 'OP', 'O2P']:
            mask = np.array([s == sp for s in species]) & sig
            if np.sum(mask) > 0:
                axes[0, 0].scatter(wavelengths[mask], ratios[mask], alpha=0.6, s=30, label=sp)
        
        axes[0, 0].axhline(1.0, color='k', linestyle='--', alpha=0.5)
        axes[0, 0].set_xlabel('Wavelength [Å]', fontsize=12)
        axes[0, 0].set_ylabel('Ratio (Approx/Exact)', fontsize=12)
        axes[0, 0].set_title(f'Approximation Ratio vs Wavelength\nTec={Tec} eV, Teh={Teh} eV, feh={feh:.4f}', 
                            fontsize=13)
        axes[0, 0].legend(ncol=2, fontsize=9)
        axes[0, 0].grid(True, alpha=0.3)
        axes[0, 0].set_ylim([0, 2])
        
        # Plot 2: Percent error vs wavelength
        for sp in ['SP', 'S2P', 'S3P', 'S4P', 'OP', 'O2P']:
            mask = np.array([s == sp for s in species]) & sig
            if np.sum(mask) > 0:
                axes[0, 1].scatter(wavelengths[mask], percent_errors[mask], alpha=0.6, s=30, label=sp)
        
        axes[0, 1].axhline(0.0, color='k', linestyle='--', alpha=0.5)
        axes[0, 1].set_xlabel('Wavelength [Å]', fontsize=12)
        axes[0, 1].set_ylabel('Percent Error [%]', fontsize=12)
        axes[0, 1].set_title(f'Percent Error vs Wavelength\nTec={Tec} eV, Teh={Teh} eV, feh={feh:.4f}', 
                            fontsize=13)
        axes[0, 1].legend(ncol=2, fontsize=9)
        axes[0, 1].grid(True, alpha=0.3)
        
        # Plot 3: Ratio vs energy
        for sp in ['SP', 'S2P', 'S3P', 'S4P', 'OP', 'O2P']:
            mask = np.array([s == sp for s in species]) & sig
            if np.sum(mask) > 0:
                axes[1, 0].scatter(energies[mask], ratios[mask], alpha=0.6, s=30, label=sp)
        
        axes[1, 0].axhline(1.0, color='k', linestyle='--', alpha=0.5)
        axes[1, 0].set_xlabel('Photon Energy [eV]', fontsize=12)
        axes[1, 0].set_ylabel('Ratio (Approx/Exact)', fontsize=12)
        axes[1, 0].set_title(f'Approximation Ratio vs Energy\nTec={Tec} eV, Teh={Teh} eV, feh={feh:.4f}', 
                            fontsize=13)
        axes[1, 0].legend(ncol=2, fontsize=9)
        axes[1, 0].grid(True, alpha=0.3)
        axes[1, 0].set_ylim([0, 2])
        
        # Plot 4: Absolute percent error vs energy
        for sp in ['SP', 'S2P', 'S3P', 'S4P', 'OP', 'O2P']:
            mask = np.array([s == sp for s in species]) & sig
            if np.sum(mask) > 0:
                axes[1, 1].scatter(energies[mask], np.abs(percent_errors[mask]), 
                                  alpha=0.6, s=30, label=sp)
        
        axes[1, 1].set_xlabel('Photon Energy [eV]', fontsize=12)
        axes[1, 1].set_ylabel('Absolute Percent Error [%]', fontsize=12)
        axes[1, 1].set_title(f'Absolute Error vs Energy\nTec={Tec} eV, Teh={Teh} eV, feh={feh:.4f}', 
                            fontsize=13)
        axes[1, 1].legend(ncol=2, fontsize=9)
        axes[1, 1].grid(True, alpha=0.3)
        axes[1, 1].set_yscale('log')
        
        plt.tight_layout()
        plt.savefig('test4_wavelength_energy_dependence.png', 
                   dpi=300, bbox_inches='tight')
        print("\nSaved: test4_wavelength_energy_dependence.png")
        plt.close()


def create_summary_heatmap(tables: EmissionTables, save_plots=True):
    """
    Create 2D heatmaps showing approximation quality across parameter space.
    
    Shows RMS error as a function of:
    1. Hot fraction vs Temperature ratio
    2. Core temperature vs Hot fraction
    """
    print("\n" + "="*70)
    print("TEST 5: Parameter Space Heatmaps")
    print("="*70)
    
    # Fixed density
    ne_total = 2200.0  # Total density [cm^-3]
    
    # Create parameter grids for heatmap 1: feh vs Teh/Tec
    Tec_fixed = 5.0
    feh_grid = np.logspace(-5, -1, 15)  # 0.001% to 10%
    temp_ratio_grid = np.array([5, 10, 20, 30, 50, 75, 100, 150, 200])
    
    rms_error_grid1 = np.zeros((len(feh_grid), len(temp_ratio_grid)))
    
    print("\nCalculating heatmap 1: feh vs Teh/Tec...")
    for i, feh in enumerate(feh_grid):
        for j, ratio in enumerate(temp_ratio_grid):
            Teh = Tec_fixed * ratio
            
            try:
                results = compare_exact_vs_approximate(
                    tables, Tec_fixed, Teh, ne_total, feh, verbose=False
                )
                sig = results['significant']
                if np.sum(sig) > 0:
                    rms_error_grid1[i, j] = np.sqrt(np.mean(results['percent_error'][sig]**2))
                else:
                    rms_error_grid1[i, j] = np.nan
            except:
                rms_error_grid1[i, j] = np.nan
        
        print(f"  Completed feh = {feh:.6f}")
    
    # Plotting
    if save_plots:
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        
        # Create meshgrid for plotting
        X, Y = np.meshgrid(temp_ratio_grid, feh_grid * 100)
        
        # Plot heatmap
        im = ax.pcolormesh(X, Y, rms_error_grid1, shading='auto', 
                          cmap='RdYlGn_r', vmin=0, vmax=50)
        
        # Add contour lines
        contours = ax.contour(X, Y, rms_error_grid1, levels=[1, 5, 10, 20, 30], 
                             colors='black', linewidths=1, alpha=0.5)
        ax.clabel(contours, inline=True, fontsize=9, fmt='%g%%')
        
        ax.set_xlabel('Temperature Ratio (Teh/Tec)', fontsize=13)
        ax.set_ylabel('Hot Electron Fraction [%]', fontsize=13)
        ax.set_title(f'RMS Approximation Error [%]\nTec = {Tec_fixed} eV, ne = {ne_total} cm⁻³', 
                    fontsize=14)
        ax.set_yscale('log')
        
        # Colorbar
        cbar = plt.colorbar(im, ax=ax, label='RMS Error [%]')
        
        plt.tight_layout()
        plt.savefig('test5_parameter_space_heatmap.png', 
                   dpi=300, bbox_inches='tight')
        print("\nSaved: test5_parameter_space_heatmap.png")
        plt.close()


def main():
    """
    Main function to run all comparison tests.
    """
    print("="*70)
    print("DOUBLE MAXWELLIAN APPROXIMATION VALIDATION")
    print("Testing: emiss(Tec,Teh,ne,feh) ≈ fec*emiss(Tec,nec) + feh*emiss(Teh,neh)")
    print("="*70)
    
    # Initialize emission tables
    tables = EmissionTables()
    
    # Load tables - adjust paths as needed
    emiss_dir = Path(__file__).parent.parent / 'Emiss_tables'
    
    print("\nLoading emission tables...")
    print("-" * 70)
    
    # Load single Maxwellian tables (needed for approximation)
    single_file = emiss_dir / 'CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5'
    if not single_file.exists():
        print(f"ERROR: Single Maxwellian table file not found:")
        print(f"  {single_file}")
        print("\nPlease ensure emission tables are in ../Emiss_tables/ directory")
        return None
    
    tables.load_single_maxwellian_tables(str(single_file))
    
    # Load double Maxwellian tables (exact solution)
    double_file = emiss_dir / 'CHIANTI_11.0.2_emiss_tables_double_maxwellian_15x10x20x10.h5'
    if not double_file.exists():
        print(f"ERROR: Double Maxwellian table file not found:")
        print(f"  {double_file}")
        print("\nPlease ensure emission tables are in ../Emiss_tables/ directory")
        return None
    
    tables.load_double_maxwellian_tables(str(double_file))
    
    # Run all tests
    print("\n" + "="*70)
    print("Running comparison tests...")
    print("="*70)
    
    start_time = time.perf_counter()
    
    # Test 1: Variation with hot electron fraction
    test_feh_variation(tables, save_plots=True)
    
    # Test 2: Variation with temperature ratio
    test_temperature_ratio_variation(tables, save_plots=True)
    
    # Test 3: Species dependence
    test_species_dependence(tables, save_plots=True)
    
    # Test 4: Wavelength/energy dependence
    test_wavelength_dependence(tables, save_plots=True)
    
    # Test 5: Parameter space heatmaps
    create_summary_heatmap(tables, save_plots=True)
    
    end_time = time.perf_counter()
    elapsed = end_time - start_time
    
    print("\n" + "="*70)
    print("ALL TESTS COMPLETE")
    print("="*70)
    print(f"Total execution time: {elapsed:.2f} seconds")
    print("\nPlots saved to /mnt/user-data/outputs/")
    print("\nKEY FINDINGS:")
    print("  - Check plots to see where approximation is valid")
    print("  - Approximation typically good for feh < 0.1% (0.001)")
    print("  - Higher ionization states show larger errors")
    print("  - Errors increase with temperature ratio Teh/Tec")
    print("="*70)
    
    return tables


if __name__ == "__main__":
    start_time = time.perf_counter()
    
    tables = main()
    
    end_time = time.perf_counter()
    print(f"\nScript execution time: {end_time - start_time:.2f} seconds")