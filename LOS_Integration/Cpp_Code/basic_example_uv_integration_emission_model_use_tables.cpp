/**
 * @file basic_example_uv_integration_emission_model_use_tables.cpp
 * @brief Example UV Emission Calculation Using Line-of-Sight Integration
 *
 * This program demonstrates how to use the Jovian UV emission raytracer to calculate
 * synthetic UV spectra through the Io Plasma Torus using proper line-of-sight
 * integration with species-by-species interpolation from CHIANTI emission tables.
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
    
    file << "# Emission Line List\n";
    file << "# Wavelength [Angstrom]    Brightness [R]    Species\n";
    file << std::fixed << std::setprecision(2);
    
    // Sort lines by brightness
    std::vector<ipt::EmissionLine> sorted_lines = lines;
    std::sort(sorted_lines.begin(), sorted_lines.end(),
              [](const ipt::EmissionLine& a, const ipt::EmissionLine& b) {
                  return a.brightness > b.brightness;
              });
    
    for (const auto& line : sorted_lines) {
        if (line.brightness > 0.01) {
            file << line.wavelength << "    " 
                 << std::scientific << std::setprecision(4) << line.brightness << "    "
                 << raytracer.get_ion_name(line.species) << "\n";
        }
    }
    
    file.close();
    std::cout << "Saved: " << filename << std::endl;
}

/**
 * @brief Generate Python plotting script.
 */
void generate_plot_script() {
    std::ofstream file("plot_spectra.py");
    if (!file.is_open()) {
        std::cerr << "Error: Could not create plot_spectra.py" << std::endl;
        return;
    }
    
    file << R"(#!/usr/bin/env python3
"""
Plot UV emission spectra from C++ raytracer output.

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
wav_s_eq, spec_s_eq = load_spectrum('spectrum_single_equatorial.dat')
wav_s_off, spec_s_off = load_spectrum('spectrum_single_off_equator.dat')
wav_d_eq, spec_d_eq = load_spectrum('spectrum_double_equatorial.dat')
wav_d_off, spec_d_off = load_spectrum('spectrum_double_off_equator.dat')

# Plot 1: Single Maxwellian equatorial
if wav_s_eq is not None:
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(wav_s_eq, spec_s_eq, 'C0-', linewidth=0.8)
    ax.set_xlabel('Wavelength [Å]', fontsize=12)
    ax.set_ylabel('Brightness [R/Å]', fontsize=12)
    ax.set_title('Single Maxwellian - Equatorial LOS (x=6 R$_J$, z=0)', fontsize=14)
    ax.set_xlim(550, 2100)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('spectrum_single_equatorial.png', dpi=150)
    plt.close()
    print('Saved: spectrum_single_equatorial.png')

# Plot 2: Single Maxwellian off-equator
if wav_s_off is not None:
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(wav_s_off, spec_s_off, 'C0-', linewidth=0.8)
    ax.set_xlabel('Wavelength [Å]', fontsize=12)
    ax.set_ylabel('Brightness [R/Å]', fontsize=12)
    ax.set_title('Single Maxwellian - Off-Equator LOS (x=6 R$_J$, z=0.5 R$_J$)', fontsize=14)
    ax.set_xlim(550, 2100)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('spectrum_single_off_equator.png', dpi=150)
    plt.close()
    print('Saved: spectrum_single_off_equator.png')

# Plot 3: Single Maxwellian comparison
if wav_s_eq is not None and wav_s_off is not None:
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(wav_s_eq, spec_s_eq, 'C0-', linewidth=0.8, label='Equatorial (z=0)')
    ax.plot(wav_s_off, spec_s_off, 'C1-', linewidth=0.8, label='Off-equator (z=0.5 R$_J$)')
    ax.set_xlabel('Wavelength [Å]', fontsize=12)
    ax.set_ylabel('Brightness [R/Å]', fontsize=12)
    ax.set_title('Single Maxwellian - Equatorial vs Off-Equator', fontsize=14)
    ax.set_xlim(550, 2100)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('spectrum_single_comparison.png', dpi=150)
    plt.close()
    print('Saved: spectrum_single_comparison.png')

# Plot 4: Double Maxwellian equatorial
if wav_d_eq is not None:
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(wav_d_eq, spec_d_eq, 'C1-', linewidth=0.8)
    ax.set_xlabel('Wavelength [Å]', fontsize=12)
    ax.set_ylabel('Brightness [R/Å]', fontsize=12)
    ax.set_title('Double Maxwellian - Equatorial LOS (x=6 R$_J$, z=0)', fontsize=14)
    ax.set_xlim(550, 2100)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('spectrum_double_equatorial.png', dpi=150)
    plt.close()
    print('Saved: spectrum_double_equatorial.png')

# Plot 5: Double Maxwellian off-equator
if wav_d_off is not None:
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(wav_d_off, spec_d_off, 'C1-', linewidth=0.8)
    ax.set_xlabel('Wavelength [Å]', fontsize=12)
    ax.set_ylabel('Brightness [R/Å]', fontsize=12)
    ax.set_title('Double Maxwellian - Off-Equator LOS (x=6 R$_J$, z=0.5 R$_J$)', fontsize=14)
    ax.set_xlim(550, 2100)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('spectrum_double_off_equator.png', dpi=150)
    plt.close()
    print('Saved: spectrum_double_off_equator.png')

# Plot 6: Single vs Double Maxwellian comparison
if wav_s_eq is not None and wav_d_eq is not None:
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(wav_s_eq, spec_s_eq, 'C0-', linewidth=0.8, label='Single Maxwellian')
    ax.plot(wav_d_eq, spec_d_eq, 'C1-', linewidth=0.8, label='Double Maxwellian')
    ax.set_xlabel('Wavelength [Å]', fontsize=12)
    ax.set_ylabel('Brightness [R/Å]', fontsize=12)
    ax.set_title('Single vs Double Maxwellian - Equatorial LOS', fontsize=14)
    ax.set_xlim(550, 2100)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('spectrum_single_vs_double.png', dpi=150)
    plt.close()
    print('Saved: spectrum_single_vs_double.png')

print('\nAll plots generated successfully.')
)";
    
    file.close();
    std::cout << "Saved: plot_spectra.py" << std::endl;
}

/**
 * @brief Print strongest emission lines.
 */
void print_strongest_lines(const std::vector<ipt::EmissionLine>& lines,
                           const ipt::JovianUVEmissionRaytracer& raytracer,
                           size_t n_lines = 10) {
    if (lines.empty()) {
        std::cout << "  No emission lines found." << std::endl;
        return;
    }
    
    std::vector<ipt::EmissionLine> sorted_lines = lines;
    std::sort(sorted_lines.begin(), sorted_lines.end(),
              [](const ipt::EmissionLine& a, const ipt::EmissionLine& b) {
                  return a.brightness > b.brightness;
              });
    
    std::cout << "\n  Strongest emission lines:" << std::endl;
    std::cout << "  " << std::setw(12) << "Wavelength" << std::setw(15) << "Brightness" 
              << "  Species" << std::endl;
    std::cout << "  " << std::setw(12) << "[Å]" << std::setw(15) << "[R]" << std::endl;
    
    size_t count = std::min(n_lines, sorted_lines.size());
    for (size_t i = 0; i < count; ++i) {
        const auto& line = sorted_lines[i];
        std::cout << "  " << std::fixed << std::setprecision(2) << std::setw(12) << line.wavelength
                  << std::scientific << std::setprecision(4) << std::setw(15) << line.brightness
                  << "  " << raytracer.get_ion_name(line.species) << std::endl;
    }
}

int main() {
    std::cout << "======================================================================\n";
    std::cout << "Jovian UV Emission Raytracer - Line-of-Sight Integration (C++)\n";
    std::cout << "======================================================================\n\n";
    
    std::cout << "This example demonstrates UV emission calculations through\n";
    std::cout << "the Io Plasma Torus using proper line-of-sight integration\n";
    std::cout << "with species-by-species interpolation from CHIANTI tables.\n\n";
    
    std::cout << "Loading plasma model and emission tables...\n";
    
    try {
        // Initialize raytracer
        ipt::JovianUVEmissionRaytracer raytracer;
        
        // Define test line-of-sight geometries
        std::array<double, 3> equatorial_pos = {6.0, -20.0, 0.0};
        std::array<double, 3> off_equator_pos = {6.0, -20.0, 0.5};
        std::array<double, 3> direction = {0.0, 1.0, 0.0};
        
        // Calculation parameters
        std::pair<double, double> wavelength_range = {550.0, 2100.0};
        double bin_width = 1.0;
        double fwhm = 6.0;
        double ds = 0.1;
        
        // =====================================================================
        // Test Case 1: Equatorial Line of Sight (Single Maxwellian)
        // =====================================================================
        std::cout << "\n======================================================================\n";
        std::cout << "Test Case 1: Equatorial Line of Sight (Single Maxwellian)\n";
        std::cout << "----------------------------------------------------------------------\n";
        std::cout << "Starting position: (6.0, -20.0, 0.0) R_J\n";
        std::cout << "Direction vector: (0.0, 1.0, 0.0)\n\n";
        
        auto start_time = std::chrono::high_resolution_clock::now();
        auto result_single_eq = raytracer.calculate_spectrum_single(
            equatorial_pos, direction, wavelength_range, bin_width, fwhm, ds);
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        // Find peak
        double peak_brightness = 0.0;
        double peak_wavelength = 0.0;
        for (size_t i = 0; i < result_single_eq.spectrum.size(); ++i) {
            if (result_single_eq.spectrum[i] > peak_brightness) {
                peak_brightness = result_single_eq.spectrum[i];
                peak_wavelength = result_single_eq.wave_bins[i];
            }
        }
        
        std::cout << "\nSingle Maxwellian Results:\n";
        std::cout << "  Total brightness: " << std::fixed << std::setprecision(2) 
                  << result_single_eq.total_brightness << " Rayleighs\n";
        std::cout << "  Peak brightness: " << std::scientific << std::setprecision(2) 
                  << peak_brightness << " R/Å at " << std::fixed << std::setprecision(2) 
                  << peak_wavelength << " Å\n";
        std::cout << "  Calculation time: " << duration.count() << " ms\n";
        
        print_strongest_lines(result_single_eq.lines, raytracer, 10);
        
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
        
        peak_brightness = 0.0;
        peak_wavelength = 0.0;
        for (size_t i = 0; i < result_single_off.spectrum.size(); ++i) {
            if (result_single_off.spectrum[i] > peak_brightness) {
                peak_brightness = result_single_off.spectrum[i];
                peak_wavelength = result_single_off.wave_bins[i];
            }
        }
        
        std::cout << "\nSingle Maxwellian Results:\n";
        std::cout << "  Total brightness: " << std::fixed << std::setprecision(2) 
                  << result_single_off.total_brightness << " Rayleighs\n";
        std::cout << "  Peak brightness: " << std::scientific << std::setprecision(2) 
                  << peak_brightness << " R/Å at " << std::fixed << std::setprecision(2) 
                  << peak_wavelength << " Å\n";
        std::cout << "  Calculation time: " << duration.count() << " ms\n";
        
        // =====================================================================
        // Test Case 3: Equatorial LOS (Double Maxwellian)
        // =====================================================================
        ipt::SpectrumResult result_double_eq;
        if (raytracer.is_double_maxwellian_loaded()) {
            std::cout << "\n======================================================================\n";
            std::cout << "Test Case 3: Equatorial LOS (Double Maxwellian)\n";
            std::cout << "----------------------------------------------------------------------\n";
            
            start_time = std::chrono::high_resolution_clock::now();
            result_double_eq = raytracer.calculate_spectrum_double(
                equatorial_pos, direction, wavelength_range, bin_width, fwhm, ds);
            end_time = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
            
            peak_brightness = 0.0;
            peak_wavelength = 0.0;
            for (size_t i = 0; i < result_double_eq.spectrum.size(); ++i) {
                if (result_double_eq.spectrum[i] > peak_brightness) {
                    peak_brightness = result_double_eq.spectrum[i];
                    peak_wavelength = result_double_eq.wave_bins[i];
                }
            }
            
            double enhancement = (result_single_eq.total_brightness > 0) ? 
                result_double_eq.total_brightness / result_single_eq.total_brightness : 0.0;
            
            std::cout << "\nDouble Maxwellian Results:\n";
            std::cout << "  Total brightness: " << std::fixed << std::setprecision(2) 
                      << result_double_eq.total_brightness << " Rayleighs\n";
            std::cout << "  Peak brightness: " << std::scientific << std::setprecision(2) 
                      << peak_brightness << " R/Å at " << std::fixed << std::setprecision(2) 
                      << peak_wavelength << " Å\n";
            std::cout << "  Enhancement factor: " << std::fixed << std::setprecision(2) 
                      << enhancement << "x\n";
            std::cout << "  Calculation time: " << duration.count() << " ms\n";
        }
        
        // =====================================================================
        // Test Case 4: Off-Equator LOS (Double Maxwellian)
        // =====================================================================
        ipt::SpectrumResult result_double_off;
        if (raytracer.is_double_maxwellian_loaded()) {
            std::cout << "\n======================================================================\n";
            std::cout << "Test Case 4: Off-Equator LOS (Double Maxwellian)\n";
            std::cout << "----------------------------------------------------------------------\n";
            
            start_time = std::chrono::high_resolution_clock::now();
            result_double_off = raytracer.calculate_spectrum_double(
                off_equator_pos, direction, wavelength_range, bin_width, fwhm, ds);
            end_time = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
            
            peak_brightness = 0.0;
            peak_wavelength = 0.0;
            for (size_t i = 0; i < result_double_off.spectrum.size(); ++i) {
                if (result_double_off.spectrum[i] > peak_brightness) {
                    peak_brightness = result_double_off.spectrum[i];
                    peak_wavelength = result_double_off.wave_bins[i];
                }
            }
            
            double enhancement = (result_single_off.total_brightness > 0) ? 
                result_double_off.total_brightness / result_single_off.total_brightness : 0.0;
            
            std::cout << "\nDouble Maxwellian Results:\n";
            std::cout << "  Total brightness: " << std::fixed << std::setprecision(2) 
                      << result_double_off.total_brightness << " Rayleighs\n";
            std::cout << "  Peak brightness: " << std::scientific << std::setprecision(2) 
                      << peak_brightness << " R/Å at " << std::fixed << std::setprecision(2) 
                      << peak_wavelength << " Å\n";
            std::cout << "  Enhancement factor: " << std::fixed << std::setprecision(2) 
                      << enhancement << "x\n";
            std::cout << "  Calculation time: " << duration.count() << " ms\n";
        }
        
        // =====================================================================
        // Save Data Files
        // =====================================================================
        std::cout << "\n======================================================================\n";
        std::cout << "Saving Data Files\n";
        std::cout << "----------------------------------------------------------------------\n";
        
        save_spectrum("spectrum_single_equatorial.dat", 
                      result_single_eq.wave_bins, result_single_eq.spectrum,
                      "Single Maxwellian equatorial spectrum (x=6 R_J, z=0)");
        
        save_spectrum("spectrum_single_off_equator.dat",
                      result_single_off.wave_bins, result_single_off.spectrum,
                      "Single Maxwellian off-equator spectrum (x=6 R_J, z=0.5 R_J)");
        
        if (raytracer.is_double_maxwellian_loaded()) {
            save_spectrum("spectrum_double_equatorial.dat",
                          result_double_eq.wave_bins, result_double_eq.spectrum,
                          "Double Maxwellian equatorial spectrum (x=6 R_J, z=0)");
            
            save_spectrum("spectrum_double_off_equator.dat",
                          result_double_off.wave_bins, result_double_off.spectrum,
                          "Double Maxwellian off-equator spectrum (x=6 R_J, z=0.5 R_J)");
        }
        
        save_line_list("lines_single_equatorial.dat", result_single_eq.lines, raytracer);
        
        generate_plot_script();
        
        // =====================================================================
        // Summary
        // =====================================================================
        std::cout << "\n======================================================================\n";
        std::cout << "Summary\n";
        std::cout << "----------------------------------------------------------------------\n";
        std::cout << "All calculations completed successfully.\n";
        std::cout << "\nTo generate plots, run: python3 plot_spectra.py\n";
        std::cout << "======================================================================\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
