/**
 * @file basic_example_optical_integration_emission_model_use_tables.cpp
 * @brief Example Optical Emission Calculation Using Line-of-Sight Integration
 *
 * This program demonstrates how to use the Jovian UV/Optical emission raytracer to
 * calculate synthetic optical spectra through the Io Plasma Torus using proper
 * line-of-sight integration with species-by-species interpolation from CHIANTI
 * emission tables.
 *
 * The example calculates optical emission spectra for test cases with both single
 * and double Maxwellian electron distributions:
 * 1. Single Maxwellian: equatorial and off-equator lines of sight
 * 2. Double Maxwellian: same geometries with hot electron population
 *
 * All lines of sight start at x=6 R_J and look through the plasma torus in the
 * +y direction, sampling the peak emission region near Io's orbit.
 *
 * @section physical_model PHYSICAL MODEL
 * - Ray tracing through 3D plasma model with trilinear interpolation
 * - Vectorized bilinear interpolation of single Maxwellian emission rates (2D tables)
 * - Vectorized quadrilinear interpolation of double Maxwellian emission rates (4D tables)
 * - Per-ion photon emission rates multiplied by ion density and 1e-6 Rayleigh factor
 * - Simpson's rule integration over line of sight
 * - Analytic ERF-based Gaussian convolution for instrument response
 *
 * @section wavelength_range WAVELENGTH RANGE
 * - Optical: 3000-10000 Å (suitable for ground-based telescopes)
 * - Key diagnostic lines: [S II] 6716/6731, [O III] 5007, [S III] 6312
 *
 * @section directory_structure DIRECTORY STRUCTURE EXPECTED
 * - IPT_Emission_MOP_Community_Code/
 *   - LOS_Integration/Cpp_Code           (this executable)
 *   - Emiss_tables/                       (CHIANTI emission tables)
 *   - 3D_Torus_Model/                     (3D plasma model in HDF5 format)
 *
 * @section output_files OUTPUT FILES
 * - optical_spectrum_single_equatorial.dat: Single Maxwellian equatorial spectrum data
 * - optical_spectrum_single_off_equator.dat: Single Maxwellian off-equator spectrum data
 * - optical_spectrum_double_equatorial.dat: Double Maxwellian equatorial spectrum data
 * - optical_spectrum_double_off_equator.dat: Double Maxwellian off-equator spectrum data
 * - optical_lines_single_equatorial.dat: Emission line list
 * - plot_optical_spectra.py: Python script to generate publication-quality plots
 *
 * @author Edward (Eddie) G. Nerney
 * @institution Laboratory for Atmospheric and Space Physics, University of Colorado Boulder
 * @date November 2025
 * @license Open source for academic and research use
 */

#include "IPT_emiss_MOP_community_code.hpp"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <cmath>

/**
 * @brief Save spectrum data to ASCII file.
 */
void save_spectrum(const std::string& filename,
                   const std::vector<double>& wavelengths,
                   const std::vector<double>& spectrum,
                   const std::string& header) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << " for writing." << std::endl;
        return;
    }
    
    file << "# " << header << "\n";
    file << "# Wavelength [Angstrom]    Brightness [R/Angstrom]\n";
    file << std::scientific << std::setprecision(6);
    
    for (size_t i = 0; i < wavelengths.size(); ++i) {
        file << wavelengths[i] << "    " << spectrum[i] << "\n";
    }
    
    file.close();
    std::cout << "Saved: " << filename << std::endl;
}

/**
 * @brief Save emission line list to ASCII file.
 */
void save_line_list(const std::string& filename,
                    const std::vector<ipt::EmissionLine>& lines,
                    const ipt::JovianUVEmissionRaytracer& raytracer) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << " for writing." << std::endl;
        return;
    }
    
    file << "# Optical Emission Line List\n";
    file << "# Wavelength [Angstrom]    Brightness [R]    Species\n";
    file << std::fixed << std::setprecision(2);
    
    // Sort lines by brightness
    std::vector<ipt::EmissionLine> sorted_lines = lines;
    std::sort(sorted_lines.begin(), sorted_lines.end(),
              [](const ipt::EmissionLine& a, const ipt::EmissionLine& b) {
                  return a.brightness > b.brightness;
              });
    
    for (const auto& line : sorted_lines) {
        if (line.brightness > 0.001) {
            file << line.wavelength << "    " 
                 << std::scientific << std::setprecision(4) << line.brightness << "    "
                 << raytracer.get_ion_name(line.species) << "\n";
        }
    }
    
    file.close();
    std::cout << "Saved: " << filename << std::endl;
}

/**
 * @brief Find brightness of line closest to target wavelength.
 */
double find_line_brightness(const std::vector<ipt::EmissionLine>& lines,
                            double target_wav, double tolerance = 5.0) {
    for (const auto& line : lines) {
        if (std::fabs(line.wavelength - target_wav) < tolerance) {
            return line.brightness;
        }
    }
    return 0.0;
}

/**
 * @brief Generate Python plotting script for optical spectra.
 */
void generate_optical_plot_script() {
    std::ofstream file("plot_optical_spectra.py");
    if (!file.is_open()) {
        std::cerr << "Error: Could not create plot_optical_spectra.py" << std::endl;
        return;
    }
    
    file << R"(#!/usr/bin/env python3
"""
Plot optical emission spectra from C++ raytracer output.

Author: Edward (Eddie) G. Nerney
Institution: LASP, University of Colorado Boulder
"""

import numpy as np
import matplotlib.pyplot as plt

def load_spectrum(filename):
    """Load spectrum data from ASCII file."""
    try:
        data = np.loadtxt(filename, comments='#')
        if data.size == 0:
            return None, None
        return data[:, 0], data[:, 1]
    except Exception as e:
        print(f"Warning: Could not load {filename}: {e}")
        return None, None

# Load all spectra
wav_s_eq, spec_s_eq = load_spectrum('optical_spectrum_single_equatorial.dat')
wav_s_off, spec_s_off = load_spectrum('optical_spectrum_single_off_equator.dat')
wav_d_eq, spec_d_eq = load_spectrum('optical_spectrum_double_equatorial.dat')
wav_d_off, spec_d_off = load_spectrum('optical_spectrum_double_off_equator.dat')

# Plot 1: Single Maxwellian equatorial
if wav_s_eq is not None:
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(wav_s_eq, spec_s_eq, 'C0-', linewidth=0.8)
    ax.set_xlabel('Wavelength [Å]', fontsize=12)
    ax.set_ylabel('Brightness [R/Å]', fontsize=12)
    ax.set_title('Single Maxwellian - Equatorial LOS - Optical (x=6 R$_J$, z=0)', fontsize=14)
    ax.set_xlim(3000, 10000)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('optical_spectrum_single_equatorial.png', dpi=150)
    plt.close()
    print('Saved: optical_spectrum_single_equatorial.png')

# Plot 2: Single Maxwellian off-equator
if wav_s_off is not None:
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(wav_s_off, spec_s_off, 'C0-', linewidth=0.8)
    ax.set_xlabel('Wavelength [Å]', fontsize=12)
    ax.set_ylabel('Brightness [R/Å]', fontsize=12)
    ax.set_title('Single Maxwellian - Off-Equator LOS - Optical (x=6 R$_J$, z=0.5 R$_J$)', fontsize=14)
    ax.set_xlim(3000, 10000)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('optical_spectrum_single_off_equator.png', dpi=150)
    plt.close()
    print('Saved: optical_spectrum_single_off_equator.png')

# Plot 3: Single Maxwellian comparison
if wav_s_eq is not None and wav_s_off is not None:
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(wav_s_eq, spec_s_eq, 'C0-', linewidth=0.8, label='Equatorial (z=0)')
    ax.plot(wav_s_off, spec_s_off, 'C1-', linewidth=0.8, alpha=0.7, label='Off-equator (z=0.5 R$_J$)')
    ax.set_xlabel('Wavelength [Å]', fontsize=12)
    ax.set_ylabel('Brightness [R/Å]', fontsize=12)
    ax.set_title('Single Maxwellian - Equatorial vs Off-Equator - Optical', fontsize=14)
    ax.set_xlim(3000, 10000)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('optical_spectrum_single_comparison.png', dpi=150)
    plt.close()
    print('Saved: optical_spectrum_single_comparison.png')

# Plot 4: Double Maxwellian equatorial
if wav_d_eq is not None:
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(wav_d_eq, spec_d_eq, 'C1-', linewidth=0.8)
    ax.set_xlabel('Wavelength [Å]', fontsize=12)
    ax.set_ylabel('Brightness [R/Å]', fontsize=12)
    ax.set_title('Double Maxwellian - Equatorial LOS - Optical (x=6 R$_J$, z=0)', fontsize=14)
    ax.set_xlim(3000, 10000)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('optical_spectrum_double_equatorial.png', dpi=150)
    plt.close()
    print('Saved: optical_spectrum_double_equatorial.png')

# Plot 5: Double Maxwellian off-equator
if wav_d_off is not None:
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(wav_d_off, spec_d_off, 'C1-', linewidth=0.8)
    ax.set_xlabel('Wavelength [Å]', fontsize=12)
    ax.set_ylabel('Brightness [R/Å]', fontsize=12)
    ax.set_title('Double Maxwellian - Off-Equator LOS - Optical (x=6 R$_J$, z=0.5 R$_J$)', fontsize=14)
    ax.set_xlim(3000, 10000)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('optical_spectrum_double_off_equator.png', dpi=150)
    plt.close()
    print('Saved: optical_spectrum_double_off_equator.png')

# Plot 6: Single vs Double Maxwellian comparison
if wav_s_eq is not None and wav_d_eq is not None:
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(wav_s_eq, spec_s_eq, 'C0-', linewidth=0.8, label='Single Maxwellian')
    ax.plot(wav_d_eq, spec_d_eq, 'C1-', linewidth=0.8, alpha=0.7, label='Double Maxwellian')
    ax.set_xlabel('Wavelength [Å]', fontsize=12)
    ax.set_ylabel('Brightness [R/Å]', fontsize=12)
    ax.set_title('Equatorial LOS: Single vs Double Maxwellian - Optical', fontsize=14)
    ax.set_xlim(3000, 10000)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('optical_spectrum_single_vs_double.png', dpi=150)
    plt.close()
    print('Saved: optical_spectrum_single_vs_double.png')

# Plot 7: [S II] doublet region (6700-6750 Å) - Density diagnostic
if wav_s_eq is not None:
    mask = (wav_s_eq >= 6700) & (wav_s_eq <= 6750)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(wav_s_eq[mask], spec_s_eq[mask], 'C0-', linewidth=1.5, 
            label='Single Maxwellian', alpha=0.8)
    if wav_d_eq is not None:
        ax.plot(wav_d_eq[mask], spec_d_eq[mask], 'C1-', linewidth=1.5,
                label='Double Maxwellian', alpha=0.8)
    ax.set_xlabel('Wavelength [Å]', fontsize=12)
    ax.set_ylabel('Brightness [R/Å]', fontsize=12)
    ax.set_title('[S II] 6716/6731 Å Doublet Region - Density Diagnostic', fontsize=14)
    ax.legend(loc='upper right', fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(6700, 6750)
    ax.set_ylim(bottom=0)
    # Line identifications
    ax.axvline(6716.4, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.axvline(6730.8, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    ymax = ax.get_ylim()[1]
    ax.text(6716.4, ymax*0.9, '[S II]\n6716 Å', ha='center', va='top', fontsize=9)
    ax.text(6730.8, ymax*0.9, '[S II]\n6731 Å', ha='center', va='top', fontsize=9)
    plt.tight_layout()
    plt.savefig('optical_sii_doublet_region.png', dpi=300)
    plt.close()
    print('Saved: optical_sii_doublet_region.png')

# Plot 8: [S III] 6312 Å region - Temperature diagnostic
if wav_s_eq is not None:
    mask = (wav_s_eq >= 6280) & (wav_s_eq <= 6330)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(wav_s_eq[mask], spec_s_eq[mask], 'C0-', linewidth=1.5,
            label='Single Maxwellian', alpha=0.8)
    if wav_d_eq is not None:
        ax.plot(wav_d_eq[mask], spec_d_eq[mask], 'C1-', linewidth=1.5,
                label='Double Maxwellian', alpha=0.8)
    ax.set_xlabel('Wavelength [Å]', fontsize=12)
    ax.set_ylabel('Brightness [R/Å]', fontsize=12)
    ax.set_title('[S III] 6312 Å Region - Temperature Diagnostic', fontsize=14)
    ax.legend(loc='upper right', fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(6280, 6330)
    ax.set_ylim(bottom=0)
    ax.axvline(6312.1, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    ymax = ax.get_ylim()[1]
    ax.text(6312.1, ymax*0.9, '[S III]\n6312 Å', ha='center', va='top', fontsize=9)
    plt.tight_layout()
    plt.savefig('optical_red_line_region.png', dpi=300)
    plt.close()
    print('Saved: optical_red_line_region.png')

# Plot 9: [O III] nebular doublet region (4940-5020 Å)
if wav_s_eq is not None:
    mask = (wav_s_eq >= 4940) & (wav_s_eq <= 5020)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(wav_s_eq[mask], spec_s_eq[mask], 'C0-', linewidth=1.5,
            label='Single Maxwellian', alpha=0.8)
    if wav_d_eq is not None:
        ax.plot(wav_d_eq[mask], spec_d_eq[mask], 'C1-', linewidth=1.5,
                label='Double Maxwellian', alpha=0.8)
    ax.set_xlabel('Wavelength [Å]', fontsize=12)
    ax.set_ylabel('Brightness [R/Å]', fontsize=12)
    ax.set_title('[O III] 4959/5007 Å Nebular Doublet Region', fontsize=14)
    ax.legend(loc='upper right', fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(4940, 5020)
    ax.set_ylim(bottom=0)
    ax.axvline(4958.9, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.axvline(5006.8, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    ymax = ax.get_ylim()[1]
    ax.text(4958.9, ymax*0.9, '[O III]\n4959 Å', ha='center', va='top', fontsize=9)
    ax.text(5006.8, ymax*0.9, '[O III]\n5007 Å', ha='center', va='top', fontsize=9)
    plt.tight_layout()
    plt.savefig('optical_oiii_doublet_region.png', dpi=300)
    plt.close()
    print('Saved: optical_oiii_doublet_region.png')

print('\nAll optical plots generated successfully.')
)";
    
    file.close();
    std::cout << "Saved: plot_optical_spectra.py" << std::endl;
}

/**
 * @brief Print strongest emission lines in optical range.
 */
void print_strongest_optical_lines(const std::vector<ipt::EmissionLine>& lines,
                                   const ipt::JovianUVEmissionRaytracer& raytracer,
                                   size_t n_lines = 20) {
    if (lines.empty()) {
        std::cout << "  No emission lines found." << std::endl;
        return;
    }
    
    std::vector<ipt::EmissionLine> sorted_lines = lines;
    std::sort(sorted_lines.begin(), sorted_lines.end(),
              [](const ipt::EmissionLine& a, const ipt::EmissionLine& b) {
                  return a.brightness > b.brightness;
              });
    
    std::cout << "\n  Strongest optical emission lines:" << std::endl;
    std::cout << "  " << std::setw(8) << "Ion" << std::setw(14) << "Wavelength" 
              << std::setw(15) << "Brightness" << std::endl;
    std::cout << "  " << std::setw(8) << "" << std::setw(14) << "[Å]" 
              << std::setw(15) << "[R]" << std::endl;
    std::cout << "  " << std::string(35, '-') << std::endl;
    
    size_t count = std::min(n_lines, sorted_lines.size());
    for (size_t i = 0; i < count; ++i) {
        const auto& line = sorted_lines[i];
        std::cout << "  " << std::setw(8) << raytracer.get_ion_name(line.species)
                  << std::fixed << std::setprecision(2) << std::setw(14) << line.wavelength
                  << std::scientific << std::setprecision(2) << std::setw(15) << line.brightness
                  << std::endl;
    }
}

int main() {
    std::cout << "======================================================================\n";
    std::cout << "OPTICAL EMISSION MODEL FOR IO PLASMA TORUS\n";
    std::cout << "Ground-Based Telescope Spectroscopy\n";
    std::cout << "Line-of-Sight Integration with CHIANTI 11.0.2\n";
    std::cout << "======================================================================\n\n";
    
    std::cout << "This example demonstrates optical emission calculations through\n";
    std::cout << "the Io Plasma Torus using proper line-of-sight integration\n";
    std::cout << "with species-by-species interpolation from CHIANTI tables.\n\n";
    
    std::cout << "Wavelength range: 3000-10000 Å (ground-based optical)\n\n";
    
    std::cout << "Loading plasma model and emission tables...\n";
    
    auto program_start = std::chrono::high_resolution_clock::now();
    
    try {
        // Initialize raytracer
        ipt::JovianUVEmissionRaytracer raytracer;
        
        // Define test line-of-sight geometries
        std::array<double, 3> equatorial_pos = {6.0, -20.0, 0.0};
        std::array<double, 3> off_equator_pos = {6.0, -20.0, 0.5};
        std::array<double, 3> direction = {0.0, 1.0, 0.0};
        
        // Optical wavelength parameters (ground-based telescope)
        std::pair<double, double> wavelength_range = {3000.0, 10000.0};
        double bin_width = 0.6;   // Spectral bin width [Angstroms] - typical for ground based
        double fwhm = 1.2;        // Instrumental FWHM [Angstroms] - ground-based spectrograph
        double ds = 0.1;          // Integration step size [R_J]
        
        // =====================================================================
        // Test Case 1: Equatorial Line of Sight (Single Maxwellian)
        // =====================================================================
        std::cout << "\n======================================================================\n";
        std::cout << "Test Case 1: Equatorial Line of Sight (Single Maxwellian)\n";
        std::cout << "----------------------------------------------------------------------\n";
        std::cout << "Starting position: (6.0, -20.0, 0.0) R_J\n";
        std::cout << "Direction vector: (0.0, 1.0, 0.0)\n";
        std::cout << "Wavelength range: 3000 - 10000 Å\n\n";
        
        auto start_time = std::chrono::high_resolution_clock::now();
        auto result_single_eq = raytracer.calculate_spectrum_single(
            equatorial_pos, direction, wavelength_range, bin_width, fwhm, ds);
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        // Find peak
        double peak_brightness_s_eq = 0.0;
        double peak_wavelength_s_eq = 0.0;
        for (size_t i = 0; i < result_single_eq.spectrum.size(); ++i) {
            if (result_single_eq.spectrum[i] > peak_brightness_s_eq) {
                peak_brightness_s_eq = result_single_eq.spectrum[i];
                peak_wavelength_s_eq = result_single_eq.wave_bins[i];
            }
        }
        
        std::cout << "\nSingle Maxwellian Results:\n";
        std::cout << "  Total brightness: " << std::fixed << std::setprecision(1) 
                  << result_single_eq.total_brightness << " Rayleighs\n";
        std::cout << "  Peak brightness: " << std::scientific << std::setprecision(2) 
                  << peak_brightness_s_eq << " R/Å at " << std::fixed << std::setprecision(1) 
                  << peak_wavelength_s_eq << " Å\n";
        std::cout << "  Number of emission lines: " << result_single_eq.lines.size() << "\n";
        std::cout << "  Calculation time: " << duration.count() << " ms\n";
        
        // =====================================================================
        // Test Case 2: Off-Equator Line of Sight (Single Maxwellian)
        // =====================================================================
        std::cout << "\n======================================================================\n";
        std::cout << "Test Case 2: Off-Equator Line of Sight (Single Maxwellian)\n";
        std::cout << "----------------------------------------------------------------------\n";
        std::cout << "Starting position: (6.0, -20.0, 0.5) R_J\n";
        std::cout << "Direction vector: (0.0, 1.0, 0.0)\n\n";
        
        start_time = std::chrono::high_resolution_clock::now();
        auto result_single_off = raytracer.calculate_spectrum_single(
            off_equator_pos, direction, wavelength_range, bin_width, fwhm, ds);
        end_time = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        double peak_brightness_s_off = 0.0;
        double peak_wavelength_s_off = 0.0;
        for (size_t i = 0; i < result_single_off.spectrum.size(); ++i) {
            if (result_single_off.spectrum[i] > peak_brightness_s_off) {
                peak_brightness_s_off = result_single_off.spectrum[i];
                peak_wavelength_s_off = result_single_off.wave_bins[i];
            }
        }
        
        std::cout << "\nSingle Maxwellian Results:\n";
        std::cout << "  Total brightness: " << std::fixed << std::setprecision(1) 
                  << result_single_off.total_brightness << " Rayleighs\n";
        std::cout << "  Peak brightness: " << std::scientific << std::setprecision(2) 
                  << peak_brightness_s_off << " R/Å at " << std::fixed << std::setprecision(1) 
                  << peak_wavelength_s_off << " Å\n";
        std::cout << "  Calculation time: " << duration.count() << " ms\n";
        
        // =====================================================================
        // Test Case 3: Equatorial LOS (Double Maxwellian)
        // =====================================================================
        ipt::SpectrumResult result_double_eq;
        double enhancement_eq = 1.0;
        if (raytracer.is_double_maxwellian_loaded()) {
            std::cout << "\n======================================================================\n";
            std::cout << "Test Case 3: Equatorial LOS (Double Maxwellian)\n";
            std::cout << "----------------------------------------------------------------------\n";
            
            start_time = std::chrono::high_resolution_clock::now();
            result_double_eq = raytracer.calculate_spectrum_double(
                equatorial_pos, direction, wavelength_range, bin_width, fwhm, ds);
            end_time = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
            
            double peak_brightness_d_eq = 0.0;
            double peak_wavelength_d_eq = 0.0;
            for (size_t i = 0; i < result_double_eq.spectrum.size(); ++i) {
                if (result_double_eq.spectrum[i] > peak_brightness_d_eq) {
                    peak_brightness_d_eq = result_double_eq.spectrum[i];
                    peak_wavelength_d_eq = result_double_eq.wave_bins[i];
                }
            }
            
            enhancement_eq = (result_single_eq.total_brightness > 0) ? 
                result_double_eq.total_brightness / result_single_eq.total_brightness : 1.0;
            
            std::cout << "\nDouble Maxwellian Results:\n";
            std::cout << "  Total brightness: " << std::fixed << std::setprecision(1) 
                      << result_double_eq.total_brightness << " Rayleighs\n";
            std::cout << "  Peak brightness: " << std::scientific << std::setprecision(2) 
                      << peak_brightness_d_eq << " R/Å at " << std::fixed << std::setprecision(1) 
                      << peak_wavelength_d_eq << " Å\n";
            std::cout << "  Enhancement factor: " << std::fixed << std::setprecision(3) 
                      << enhancement_eq << "\n";
            std::cout << "  Calculation time: " << duration.count() << " ms\n";
        }
        
        // =====================================================================
        // Test Case 4: Off-Equator LOS (Double Maxwellian)
        // =====================================================================
        ipt::SpectrumResult result_double_off;
        double enhancement_off = 1.0;
        if (raytracer.is_double_maxwellian_loaded()) {
            std::cout << "\n======================================================================\n";
            std::cout << "Test Case 4: Off-Equator LOS (Double Maxwellian)\n";
            std::cout << "----------------------------------------------------------------------\n";
            
            start_time = std::chrono::high_resolution_clock::now();
            result_double_off = raytracer.calculate_spectrum_double(
                off_equator_pos, direction, wavelength_range, bin_width, fwhm, ds);
            end_time = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
            
            double peak_brightness_d_off = 0.0;
            double peak_wavelength_d_off = 0.0;
            for (size_t i = 0; i < result_double_off.spectrum.size(); ++i) {
                if (result_double_off.spectrum[i] > peak_brightness_d_off) {
                    peak_brightness_d_off = result_double_off.spectrum[i];
                    peak_wavelength_d_off = result_double_off.wave_bins[i];
                }
            }
            
            enhancement_off = (result_single_off.total_brightness > 0) ? 
                result_double_off.total_brightness / result_single_off.total_brightness : 1.0;
            
            std::cout << "\nDouble Maxwellian Results:\n";
            std::cout << "  Total brightness: " << std::fixed << std::setprecision(1) 
                      << result_double_off.total_brightness << " Rayleighs\n";
            std::cout << "  Peak brightness: " << std::scientific << std::setprecision(2) 
                      << peak_brightness_d_off << " R/Å at " << std::fixed << std::setprecision(1) 
                      << peak_wavelength_d_off << " Å\n";
            std::cout << "  Enhancement factor: " << std::fixed << std::setprecision(3) 
                      << enhancement_off << "\n";
            std::cout << "  Calculation time: " << duration.count() << " ms\n";
        }
        
        // =====================================================================
        // Save Data Files
        // =====================================================================
        std::cout << "\n======================================================================\n";
        std::cout << "Saving Data Files\n";
        std::cout << "----------------------------------------------------------------------\n";
        
        save_spectrum("optical_spectrum_single_equatorial.dat", 
                      result_single_eq.wave_bins, result_single_eq.spectrum,
                      "Single Maxwellian equatorial optical spectrum (x=6 R_J, z=0)");
        
        save_spectrum("optical_spectrum_single_off_equator.dat",
                      result_single_off.wave_bins, result_single_off.spectrum,
                      "Single Maxwellian off-equator optical spectrum (x=6 R_J, z=0.5 R_J)");
        
        if (raytracer.is_double_maxwellian_loaded()) {
            save_spectrum("optical_spectrum_double_equatorial.dat",
                          result_double_eq.wave_bins, result_double_eq.spectrum,
                          "Double Maxwellian equatorial optical spectrum (x=6 R_J, z=0)");
            
            save_spectrum("optical_spectrum_double_off_equator.dat",
                          result_double_off.wave_bins, result_double_off.spectrum,
                          "Double Maxwellian off-equator optical spectrum (x=6 R_J, z=0.5 R_J)");
        }
        
        save_line_list("optical_lines_single_equatorial.dat", result_single_eq.lines, raytracer);
        
        generate_optical_plot_script();
        
        // =====================================================================
        // Strongest Optical Emission Lines
        // =====================================================================
        std::cout << "\n======================================================================\n";
        std::cout << "Strongest Optical Emission Lines (Equatorial, Single Maxwellian)\n";
        std::cout << "----------------------------------------------------------------------\n";
        
        print_strongest_optical_lines(result_single_eq.lines, raytracer, 20);
        
        // =====================================================================
        // Key Diagnostic Line Ratios
        // =====================================================================
        std::cout << "\n======================================================================\n";
        std::cout << "Key Diagnostic Line Ratios\n";
        std::cout << "----------------------------------------------------------------------\n";
        
        // [S II] 6716/6731 ratio (density diagnostic)
        double sii_6716 = find_line_brightness(result_single_eq.lines, 6716.4, 3.0);
        double sii_6731 = find_line_brightness(result_single_eq.lines, 6730.8, 3.0);
        if (sii_6731 > 0) {
            double sii_ratio = sii_6716 / sii_6731;
            std::cout << "[S II] 6716/6731 ratio (density diagnostic):\n";
            std::cout << "  Single Maxwellian: " << std::fixed << std::setprecision(3) 
                      << sii_ratio << "\n";
            std::cout << "  (Low density limit ~1.5, high density limit ~0.44)\n";
        }
        
        // [O III] temperature diagnostic
        double oiii_5007 = find_line_brightness(result_single_eq.lines, 5006.8, 3.0);
        double oiii_4959 = find_line_brightness(result_single_eq.lines, 4958.9, 3.0);
        if (oiii_4959 > 0) {
            double oiii_ratio = oiii_5007 / oiii_4959;
            std::cout << "\n[O III] 5007/4959 ratio (should be ~3):\n";
            std::cout << "  Single Maxwellian: " << std::fixed << std::setprecision(3) 
                      << oiii_ratio << "\n";
        }
        
        // =====================================================================
        // Physical Interpretation
        // =====================================================================
        std::cout << "\n======================================================================\n";
        std::cout << "Physical Interpretation - Optical Observations\n";
        std::cout << "----------------------------------------------------------------------\n";
        std::cout << "The plasma torus optical emission shows:\n";
        std::cout << "- Peak emission near 6 R_J (Io's orbital radius)\n";
        std::cout << "- Scale height of ~0.5-1 R_J\n";
        std::cout << "- Dominant optical emission from [S II] and [O II] forbidden lines\n";
        std::cout << "- [S II] 6716/6731 doublet ratio sensitive to electron density\n";
        std::cout << "- [O III] 5007/4959 ratio fixed by atomic physics (~3:1)\n";
        std::cout << "- [S III] 6312 Å line useful for temperature diagnostics\n";
        std::cout << "- Temperature ~5-10 eV in the cold torus\n";
        std::cout << "- Electron density ~100-2000 cm^-3 at peak\n";
        if (raytracer.is_double_maxwellian_loaded()) {
            std::cout << "- Hot electrons enhance high-excitation optical lines\n";
            std::cout << "- Overall brightness enhancement: " << std::fixed << std::setprecision(1)
                      << (enhancement_eq - 1.0) * 100.0 << "%\n";
        }
        
        // =====================================================================
        // Timing and Summary
        // =====================================================================
        auto program_end = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(program_end - program_start);
        
        std::cout << "\n======================================================================\n";
        std::cout << "Optical Emission Calculation Complete\n";
        std::cout << "======================================================================\n\n";
        std::cout << "Execution time: " << std::fixed << std::setprecision(2) 
                  << total_duration.count() / 1000.0 << " seconds\n\n";
        std::cout << "Output files generated:\n";
        std::cout << "  - optical_spectrum_single_equatorial.dat\n";
        std::cout << "  - optical_spectrum_single_off_equator.dat\n";
        if (raytracer.is_double_maxwellian_loaded()) {
            std::cout << "  - optical_spectrum_double_equatorial.dat\n";
            std::cout << "  - optical_spectrum_double_off_equator.dat\n";
        }
        std::cout << "  - optical_lines_single_equatorial.dat\n";
        std::cout << "  - plot_optical_spectra.py\n\n";
        
        std::cout << "Calculation methodology:\n";
        std::cout << "  - Ray tracing through 3D plasma model with trilinear interpolation\n";
        std::cout << "  - Vectorized bilinear interpolation of single Maxwellian emission rates (2D)\n";
        if (raytracer.is_double_maxwellian_loaded()) {
            std::cout << "  - Vectorized quadrilinear interpolation of double Maxwellian emission rates (4D)\n";
        }
        std::cout << "  - Per-ion photon emission rates × ion density × 1e-6 → Rayleighs\n";
        std::cout << "  - Simpson's rule for line-of-sight integration\n";
        std::cout << "  - ERF-based Gaussian convolution for instrument response\n";
        std::cout << "  - Integration step size: ds = " << ds << " R_J\n";
        std::cout << "  - Instrumental FWHM: " << fwhm << " Å (R~3000 ground-based)\n";
        std::cout << "  - Wavelength range: " << wavelength_range.first << "-" 
                  << wavelength_range.second << " Å\n\n";
        
        std::cout << "To generate plots, run: python3 plot_optical_spectra.py\n";
        std::cout << "======================================================================\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
