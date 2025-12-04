/**
 * @file IPT_emiss_MOP_community_code.hpp
 * @brief Io Plasma Torus (IPT) UV/Optical Emission Line-of-Sight Integration Module
 *
 * This header provides comprehensive functionality for calculating UV and optical emission
 * spectra from the Io Plasma Torus (IPT) using ray tracing through a 3D plasma model with
 * CHIANTI atomic emission tables. Supports both single and double Maxwellian electron
 * distributions with proper species-by-species interpolation and vectorized computations.
 *
 * @section physical_model PHYSICAL MODEL
 * - Plasma torus extends from ~5-10 R_J in cylindrical radius
 * - Peak emission near 6 R_J at centrifugal equator
 * - Scale height ~0.5-1 R_J
 * - Based on IPT Isotropic Model interpolated to rectilinear grid
 * - Per-ion photon emission rates from CHIANTI 11.0.2 atomic database
 * - Line-of-sight integration using Simpson's rule
 * - Optically thin plasma approximation (valid for IPT)
 *
 * @section computational COMPUTATIONAL APPROACH
 * - Ray tracing through 3D plasma model with trilinear interpolation
 * - Bilinear interpolation of single Maxwellian emission tables (2D: Te, ne)
 * - Quadrilinear interpolation of double Maxwellian emission tables (4D: Tec, Teh, ne, feh)
 * - Vectorized emission rate interpolation for all lines simultaneously
 * - Per-ion photon emission rates multiplied by ion density and 1e-6 Rayleigh factor
 * - Simpson's rule integration over line of sight
 * - Analytic ERF-based Gaussian convolution for instrument response
 *
 * @section units COORDINATE SYSTEMS AND UNITS
 * - Positions: Jupiter radii [R_J]
 * - Temperatures: electron volts [eV]
 * - Densities: particles per cubic centimeter [cm^-3]
 * - Wavelengths: Angstroms [Å]
 * - Emission rates: photons per second per ion [photons s^-1 ion^-1]
 * - Brightnesses: Rayleighs [R], where 1 R = 10^6 photons s^-1 cm^-2 (4π sr)^-1
 *
 * @section wavelengths APPLICABLE WAVELENGTH RANGES
 * - UV instruments: 550-2100 Å (JUICE-UVS, Europa-UVS, HST/STIS)
 * - Optical instruments: 3000-10000 Å (ground-based telescopes, HST optical)
 *
 * @section references REFERENCES
 * - CHIANTI database: Dere et al. 1997; Del Zanna et al. 2020; Dufresne et al. 2024
 * - IPT observations: Steffl et al. 2004a,b; Thomas et al. 2004; Bagenal & Delamere 2011
 * - Emission modeling: Nerney et al. 2017, 2020, 2022, 2025a, 2025b
 * - Electron distributions: Meyer-Vernet & Moncuquet 1989; Moncuquet et al. 2002
 *
 * @author Edward (Eddie) G. Nerney
 * @institution Laboratory for Atmospheric and Space Physics, University of Colorado Boulder
 * @license Open source for academic and research use
 * @version 2.0
 * @date November 2025
 *
 * @section acknowledgment CHIANTI ACKNOWLEDGMENT
 * CHIANTI is a collaborative project involving George Mason University, the University
 * of Michigan (USA), University of Cambridge (UK), and NASA Goddard Space Flight Center (USA).
 */

#ifndef IPT_EMISS_MOP_COMMUNITY_CODE_HPP
#define IPT_EMISS_MOP_COMMUNITY_CODE_HPP

#include <vector>
#include <array>
#include <string>
#include <unordered_map>
#include <memory>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <filesystem>

// =============================================================================
// PHYSICAL CONSTANTS
// =============================================================================

namespace ipt {

/// Jupiter radius in kilometers
constexpr double R_J = 71492.0;

/// Jupiter radius in centimeters
constexpr double R_J_CM = 7.1492e9;

/// Conversion factor for integration to Rayleighs: 1e-6 * R_J_CM
constexpr double RAYLEIGH_FACTOR = 1.0e-6 * R_J_CM;

/// Natural log of 2 for FWHM conversion
constexpr double LN2 = 0.693147180559945309417232121458;

/// Square root of 2
constexpr double SQRT2 = 1.41421356237309504880168872420969807857;

/// Square root of 2 * ln(2): sqrt(2 * 0.693147...) = 1.17741002251547469
constexpr double SQRT_2LN2 = 1.17741002251547469101156932645969;

/// FWHM to sigma conversion factor: 1 / (2 * sqrt(2 * ln(2)))
/// Pre-computed value since std::sqrt is not constexpr in C++17
constexpr double FWHM_TO_SIGMA = 0.42466090014400952;

// =============================================================================
// FORWARD DECLARATIONS
// =============================================================================

class JovianUVEmissionRaytracer;

// =============================================================================
// DATA STRUCTURES
// =============================================================================

/**
 * @struct EmissionLine
 * @brief Represents a single emission line with wavelength, brightness, and species.
 */
struct EmissionLine {
    double wavelength;      ///< Line wavelength in Angstroms
    double brightness;      ///< Integrated brightness in Rayleighs
    std::string species;    ///< Species identifier (e.g., "SP", "S2P", "OP")
    
    EmissionLine(double w, double b, const std::string& s)
        : wavelength(w), brightness(b), species(s) {}
};

/**
 * @struct SpectrumResult
 * @brief Contains the results of a spectrum calculation.
 */
struct SpectrumResult {
    std::vector<double> wave_bins;      ///< Wavelength bin centers [Angstroms]
    std::vector<double> spectrum;       ///< Convolved spectrum [Rayleighs/Angstrom]
    std::vector<EmissionLine> lines;    ///< Individual emission lines
    double total_brightness;            ///< Total integrated brightness [Rayleighs]
};

/**
 * @struct PlasmaModel
 * @brief Contains the 3D plasma model data and coordinate axes.
 */
struct PlasmaModel {
    // Coordinate axes in Jupiter radii
    std::vector<double> x_axis;
    std::vector<double> y_axis;
    std::vector<double> z_axis;
    
    // Grid dimensions
    size_t nx, ny, nz;
    
    // Electron density fields [cm^-3] stored as flat arrays
    std::vector<double> nec;        ///< Cold/core electron density
    std::vector<double> neh;        ///< Hot electron density
    std::vector<double> ne_total;   ///< Total electron density
    
    // Ion density fields [cm^-3]
    std::vector<double> nsp;        ///< S+ density
    std::vector<double> ns2p;       ///< S++ density
    std::vector<double> ns3p;       ///< S+++ density
    std::vector<double> nop;        ///< O+ density
    std::vector<double> no2p;       ///< O++ density
    
    // Temperature fields [eV]
    std::vector<double> Tec;        ///< Cold electron temperature
    std::vector<double> Teh;        ///< Hot electron temperature
    
    // Hot electron fraction
    std::vector<double> feh;        ///< Hot electron fraction (ne_h / ne_total)
    
    /**
     * @brief Get linear index from 3D indices.
     * 
     * Uses C-order (row-major) indexing to match HDF5/numpy storage format.
     * For array with shape (nx, ny, nz), the flat index for position (i, j, k) is:
     * flat_index = i * (ny * nz) + j * nz + k
     * 
     * @param i X index (slowest varying)
     * @param j Y index
     * @param k Z index (fastest varying)
     * @return Linear index into flat arrays
     */
    inline size_t index(size_t i, size_t j, size_t k) const {
        return (i * ny + j) * nz + k;
    }
};

/**
 * @struct EmissionTableSingle
 * @brief Single Maxwellian emission table data for one species.
 */
struct EmissionTableSingle {
    std::vector<double> wavelengths;    ///< Emission line wavelengths [Angstroms]
    std::vector<double> emissivities;   ///< Emission rates [photons/s/ion], shape: (n_lines, n_T, n_n)
    size_t n_lines;                     ///< Number of emission lines
    size_t n_T;                         ///< Number of temperature grid points
    size_t n_n;                         ///< Number of density grid points
    
    /**
     * @brief Get emissivity at (line_idx, T_idx, n_idx)
     */
    inline double get(size_t line_idx, size_t T_idx, size_t n_idx) const {
        return emissivities[line_idx * n_T * n_n + T_idx * n_n + n_idx];
    }
};

/**
 * @struct EmissionTableDouble
 * @brief Double Maxwellian emission table data for one species.
 */
struct EmissionTableDouble {
    std::vector<double> wavelengths;    ///< Emission line wavelengths [Angstroms]
    std::vector<double> emissivities;   ///< Emission rates, shape: (n_lines, n_Tec, n_Teh, n_ne, n_feh)
    size_t n_lines;                     ///< Number of emission lines
    size_t n_Tec;                       ///< Number of core temperature grid points
    size_t n_Teh;                       ///< Number of hot temperature grid points
    size_t n_ne;                        ///< Number of density grid points
    size_t n_feh;                       ///< Number of hot fraction grid points
    
    /**
     * @brief Get emissivity at (line_idx, Tec_idx, Teh_idx, ne_idx, feh_idx)
     */
    inline double get(size_t line_idx, size_t Tec_idx, size_t Teh_idx, 
                      size_t ne_idx, size_t feh_idx) const {
        size_t stride_feh = 1;
        size_t stride_ne = n_feh;
        size_t stride_Teh = n_ne * n_feh;
        size_t stride_Tec = n_Teh * n_ne * n_feh;
        size_t stride_line = n_Tec * n_Teh * n_ne * n_feh;
        return emissivities[line_idx * stride_line + Tec_idx * stride_Tec + 
                           Teh_idx * stride_Teh + ne_idx * stride_ne + feh_idx * stride_feh];
    }
};

// =============================================================================
// UTILITY FUNCTIONS
// =============================================================================

/**
 * @brief Check if a double value is valid (finite and not NaN).
 * 
 * This function provides a robust check that works correctly regardless
 * of compiler optimization flags like -ffast-math.
 * 
 * @param x Value to check
 * @return true if value is finite and not NaN
 */
inline bool is_valid_double(double x) {
    // Use direct comparison - a NaN is never equal to itself
    // This works even with -ffast-math on most compilers
    return (x == x) && (x != std::numeric_limits<double>::infinity()) 
                    && (x != -std::numeric_limits<double>::infinity());
}

/**
 * @brief Scalar error function for single values.
 * 
 * Uses Abramowitz and Stegun approximation (equation 7.1.26) with maximum error < 1.5e-7.
 * 
 * @param x Input value
 * @return Error function value erf(x)
 */
inline double erf_scalar(double x) {
    // Constants for Abramowitz and Stegun approximation (equation 7.1.26)
    constexpr double a1 =  0.254829592;
    constexpr double a2 = -0.284496736;
    constexpr double a3 =  1.421413741;
    constexpr double a4 = -1.453152027;
    constexpr double a5 =  1.061405429;
    constexpr double p  =  0.3275911;
    
    double sign = (x >= 0.0) ? 1.0 : -1.0;
    double abs_x = std::fabs(x);
    
    // Clamp to avoid overflow in exp
    if (abs_x > 6.0) {
        return sign * 1.0;
    }
    
    double t = 1.0 / (1.0 + p * abs_x);
    double t2 = t * t;
    double t3 = t2 * t;
    double t4 = t3 * t;
    double t5 = t4 * t;
    
    double y = 1.0 - (a1*t + a2*t2 + a3*t3 + a4*t4 + a5*t5) * std::exp(-abs_x * abs_x);
    
    return sign * y;
}

/**
 * @brief Simpson's rule integration for uniformly spaced data.
 * 
 * Implements composite Simpson's 1/3 rule matching scipy.integrate.simpson behavior.
 * For odd number of points (even intervals), uses pure Simpson's rule.
 * For even number of points (odd intervals), uses Simpson's rule plus trapezoidal for last interval.
 * 
 * @param y Function values at each point
 * @param dx Uniform spacing between x points
 * @return Integral value
 */
inline double simpson_integrate_uniform(const std::vector<double>& y, double dx) {
    size_t n = y.size();
    
    if (n < 2) return 0.0;
    if (n == 2) {
        // Trapezoidal rule for single interval
        return 0.5 * (y[0] + y[1]) * dx;
    }
    
    double result = 0.0;
    
    // For odd number of points (even number of intervals), use pure Simpson's
    // For even number of points (odd number of intervals), handle last interval separately
    size_t n_simpson = (n % 2 == 1) ? n : n - 1;
    
    if (n_simpson >= 3) {
        // Simpson's 1/3 rule: integral = (dx/3) * [y0 + 4*y1 + 2*y2 + 4*y3 + ... + yn]
        result = y[0] + y[n_simpson - 1];
        
        for (size_t i = 1; i < n_simpson - 1; i += 2) {
            result += 4.0 * y[i];
        }
        for (size_t i = 2; i < n_simpson - 1; i += 2) {
            result += 2.0 * y[i];
        }
        
        result *= dx / 3.0;
    }
    
    // Handle remaining interval with trapezoidal rule if n is even
    if (n % 2 == 0) {
        result += 0.5 * (y[n - 2] + y[n - 1]) * dx;
    }
    
    return result;
}

/**
 * @brief Binary search to find bracketing indices in sorted array.
 * 
 * Returns the index i such that arr[i-1] <= val < arr[i], clamped to valid range.
 * 
 * @param arr Sorted array
 * @param val Value to search for
 * @return Upper bracket index (clamped to [1, arr.size()-1])
 */
inline size_t searchsorted(const std::vector<double>& arr, double val) {
    if (arr.empty()) return 0;
    auto it = std::lower_bound(arr.begin(), arr.end(), val);
    size_t idx = static_cast<size_t>(std::distance(arr.begin(), it));
    if (idx == 0) idx = 1;
    if (idx >= arr.size()) idx = arr.size() - 1;
    return idx;
}

// =============================================================================
// MAIN RAYTRACER CLASS
// =============================================================================

/**
 * @class JovianUVEmissionRaytracer
 * @brief Main class for calculating UV and optical emission through the Jovian plasma torus.
 *
 * This class provides a complete workflow for line-of-sight integrated emission
 * calculations using ray tracing through a 3D plasma model combined with CHIANTI
 * atomic emission tables. Supports both single and double Maxwellian electron
 * distributions with optimized vectorized computations.
 *
 * @section usage Usage Example
 * @code
 * JovianUVEmissionRaytracer raytracer;
 * 
 * // Define line of sight
 * std::array<double, 3> start_pos = {6.0, -20.0, 0.0};
 * std::array<double, 3> direction = {0.0, 1.0, 0.0};
 * 
 * // Calculate spectrum
 * auto result = raytracer.calculate_spectrum_single(start_pos, direction);
 * @endcode
 */
class JovianUVEmissionRaytracer {
public:
    // =========================================================================
    // CONSTRUCTORS AND INITIALIZATION
    // =========================================================================
    
    /**
     * @brief Default constructor with automatic file path detection.
     */
    JovianUVEmissionRaytracer();
    
    /**
     * @brief Constructor with explicit file paths.
     */
    JovianUVEmissionRaytracer(const std::string& plasma_file,
                              const std::string& emission_file_single,
                              const std::string& emission_file_double = "");
    
    // =========================================================================
    // DATA LOADING METHODS
    // =========================================================================
    
    void load_plasma_model(const std::string& filename);
    void load_emission_tables_single(const std::string& filename);
    void load_emission_tables_double(const std::string& filename);
    
    // =========================================================================
    // RAY TRACING METHODS
    // =========================================================================
    
    std::tuple<std::vector<double>, std::vector<std::array<double, 3>>>
    trace_ray(const std::array<double, 3>& start_pos,
              const std::array<double, 3>& direction,
              double ds = 0.1,
              double max_distance = 40.0) const;
    
    // =========================================================================
    // INTERPOLATION METHODS
    // =========================================================================
    
    double interpolate_3d(const std::vector<double>& field,
                          double x, double y, double z) const;
    
    std::vector<double> interpolate_3d_vectorized(
        const std::vector<double>& field,
        const std::vector<std::array<double, 3>>& positions) const;
    
    std::tuple<std::vector<double>, std::vector<std::vector<double>>>
    interpolate_emission_rates_single_vectorized(
        const std::vector<double>& Te_arr,
        const std::vector<double>& ne_arr,
        const std::string& species_key,
        double min_wav = 550.0,
        double max_wav = 2100.0) const;
    
    std::tuple<std::vector<double>, std::vector<std::vector<double>>>
    interpolate_emission_rates_double_vectorized(
        const std::vector<double>& Tec_arr,
        const std::vector<double>& Teh_arr,
        const std::vector<double>& ne_arr,
        const std::vector<double>& feh_arr,
        const std::string& species_key,
        double min_wav = 550.0,
        double max_wav = 2100.0) const;
    
    // =========================================================================
    // LINE-OF-SIGHT INTEGRATION METHODS
    // =========================================================================
    
    std::tuple<std::vector<double>, std::vector<double>>
    integrate_species_emission_single(
        const std::vector<double>& s_values,
        const std::string& species_key,
        const std::vector<double>& Te_los,
        const std::vector<double>& ne_los,
        const std::vector<double>& n_ion_los,
        double min_wav = 550.0,
        double max_wav = 2100.0) const;
    
    std::tuple<std::vector<double>, std::vector<double>>
    integrate_species_emission_double(
        const std::vector<double>& s_values,
        const std::string& species_key,
        const std::vector<double>& Tec_los,
        const std::vector<double>& Teh_los,
        const std::vector<double>& ne_los,
        const std::vector<double>& feh_los,
        const std::vector<double>& n_ion_los,
        double min_wav = 550.0,
        double max_wav = 2100.0) const;
    
    // =========================================================================
    // INSTRUMENT RESPONSE CONVOLUTION
    // =========================================================================
    
    std::vector<double> convolve_spectrum_erf(
        const std::vector<double>& wavelength_grid,
        double bin_width,
        const std::vector<double>& line_wavelengths,
        const std::vector<double>& line_brightnesses,
        double fwhm = 6.0) const;
    
    // =========================================================================
    // MAIN SPECTRUM CALCULATION METHODS
    // =========================================================================
    
    SpectrumResult calculate_spectrum_single(
        const std::array<double, 3>& slit_pos_vec,
        const std::array<double, 3>& norm_vec,
        std::pair<double, double> wavelength_range = {550.0, 2100.0},
        double bin_width = 1.0,
        double fwhm = 6.0,
        double ds = 0.1) const;
    
    SpectrumResult calculate_spectrum_double(
        const std::array<double, 3>& slit_pos_vec,
        const std::array<double, 3>& norm_vec,
        std::pair<double, double> wavelength_range = {550.0, 2100.0},
        double bin_width = 1.0,
        double fwhm = 6.0,
        double ds = 0.1) const;
    
    // =========================================================================
    // ACCESSORS
    // =========================================================================
    
    bool is_single_maxwellian_loaded() const { return single_maxwellian_loaded_; }
    bool is_double_maxwellian_loaded() const { return double_maxwellian_loaded_; }
    std::string get_ion_name(const std::string& species_key) const;
    const std::vector<std::string>& get_default_species() const { return default_species_; }
    const PlasmaModel& get_plasma_model() const { return plasma_model_; }
    
private:
    PlasmaModel plasma_model_;
    
    // Single Maxwellian emission data
    std::vector<double> temp_arr_;
    std::vector<double> dens_arr_;
    std::vector<double> log_temp_;
    std::vector<double> log_dens_;
    std::unordered_map<std::string, EmissionTableSingle> emiss_single_;
    
    // Double Maxwellian emission data
    std::vector<double> tec_arr_;
    std::vector<double> teh_arr_;
    std::vector<double> ne_arr_;
    std::vector<double> feh_arr_;
    std::vector<double> log_tec_;
    std::vector<double> log_teh_;
    std::vector<double> log_ne_;
    std::vector<double> log_feh_;
    std::unordered_map<std::string, EmissionTableDouble> emiss_double_;
    
    bool single_maxwellian_loaded_ = false;
    bool double_maxwellian_loaded_ = false;
    
    std::vector<std::string> default_species_ = {"SP", "S2P", "S3P", "OP", "O2P"};
    std::unordered_map<std::string, std::string> ion_names_ = {
        {"SP", "S II"}, {"S2P", "S III"}, {"S3P", "S IV"},
        {"S4P", "S V"}, {"OP", "O II"}, {"O2P", "O III"}
    };
    
    std::string get_default_plasma_path() const;
    std::string get_default_emission_path_single() const;
    std::string get_default_emission_path_double() const;
    const std::vector<double>& get_ion_field(const std::string& species_key) const;
};

/**
 * @brief Standalone function for ERF-based spectrum convolution.
 */
std::vector<double> simulate_ipt_spectrum_rayleighs_erf_form(
    const std::vector<double>& wavelength_grid,
    double bin_width,
    const std::vector<double>& line_wavelengths,
    const std::vector<double>& line_brightnesses,
    double fwhm = 6.0);

} // namespace ipt

#endif // IPT_EMISS_MOP_COMMUNITY_CODE_HPP
