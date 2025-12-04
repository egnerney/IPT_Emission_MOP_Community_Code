/**
 * @file IPT_emiss_MOP_community_code.cpp
 * @brief Implementation of Io Plasma Torus UV/Optical Emission Line-of-Sight Integration
 *
 * @author Edward (Eddie) G. Nerney
 * @institution Laboratory for Atmospheric and Space Physics, University of Colorado Boulder
 * @license Open source for academic and research use
 * @version 2.0
 * @date November 2025
 */

#include "IPT_emiss_MOP_community_code.hpp"
#include <H5Cpp.h>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <sstream>
#include <fstream>

namespace ipt {

// =============================================================================
// CONSTRUCTORS
// =============================================================================

JovianUVEmissionRaytracer::JovianUVEmissionRaytracer() {
    std::cout << "Initializing Jovian UV/Optical Emission Raytracer..." << std::endl;
    
    std::string plasma_file = get_default_plasma_path();
    std::string emission_file_single = get_default_emission_path_single();
    std::string emission_file_double = get_default_emission_path_double();
    
    load_plasma_model(plasma_file);
    load_emission_tables_single(emission_file_single);
    
    if (!emission_file_double.empty()) {
        try {
            load_emission_tables_double(emission_file_double);
        } catch (const std::exception& e) {
            std::cerr << "Warning: Could not load double Maxwellian tables: " 
                      << e.what() << std::endl;
        }
    }
    
    std::cout << "Initialization complete." << std::endl;
}

JovianUVEmissionRaytracer::JovianUVEmissionRaytracer(
    const std::string& plasma_file,
    const std::string& emission_file_single,
    const std::string& emission_file_double) {
    
    std::cout << "Initializing Jovian UV/Optical Emission Raytracer..." << std::endl;
    
    load_plasma_model(plasma_file);
    load_emission_tables_single(emission_file_single);
    
    if (!emission_file_double.empty()) {
        try {
            load_emission_tables_double(emission_file_double);
        } catch (const std::exception& e) {
            std::cerr << "Warning: Could not load double Maxwellian tables: " 
                      << e.what() << std::endl;
        }
    }
    
    std::cout << "Initialization complete." << std::endl;
}

// =============================================================================
// PATH HELPERS
// =============================================================================

std::string JovianUVEmissionRaytracer::get_default_plasma_path() const {
    std::filesystem::path current = std::filesystem::current_path();
    std::filesystem::path plasma_path = current / ".." / ".." / "3D_Torus_Model" / 
                                        "jovian_plasma_interpolated_381x381x231.h5";
    
    if (std::filesystem::exists(plasma_path)) {
        return std::filesystem::canonical(plasma_path).string();
    }
    
    throw std::runtime_error("Plasma model not found at expected location.");
}

std::string JovianUVEmissionRaytracer::get_default_emission_path_single() const {
    std::filesystem::path current = std::filesystem::current_path();
    std::filesystem::path emiss_path = current / ".." / ".." / "Emiss_tables" / 
                                       "CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5";
    
    if (std::filesystem::exists(emiss_path)) {
        return std::filesystem::canonical(emiss_path).string();
    }
    
    throw std::runtime_error("Single Maxwellian emission tables not found.");
}

std::string JovianUVEmissionRaytracer::get_default_emission_path_double() const {
    std::filesystem::path current = std::filesystem::current_path();
    std::filesystem::path emiss_path = current / ".." / ".." / "Emiss_tables" / 
                                       "CHIANTI_11.0.2_emiss_tables_double_maxwellian_16x8x12x10.h5";
    
    if (std::filesystem::exists(emiss_path)) {
        return std::filesystem::canonical(emiss_path).string();
    }
    
    return "";
}

// =============================================================================
// DATA LOADING METHODS
// =============================================================================

void JovianUVEmissionRaytracer::load_plasma_model(const std::string& filename) {
    std::cout << "Loading plasma model from " << filename << "..." << std::endl;
    
    try {
        H5::H5File file(filename, H5F_ACC_RDONLY);
        
        // Helper to read 1D double array
        auto read_1d = [&file](const std::string& path) -> std::vector<double> {
            H5::DataSet dataset = file.openDataSet(path);
            H5::DataSpace dataspace = dataset.getSpace();
            hsize_t dims;
            dataspace.getSimpleExtentDims(&dims);
            std::vector<double> data(dims);
            dataset.read(data.data(), H5::PredType::NATIVE_DOUBLE);
            return data;
        };
        
        // Helper to read 3D double array
        auto read_3d = [&file](const std::string& path) -> std::vector<double> {
            H5::DataSet dataset = file.openDataSet(path);
            H5::DataSpace dataspace = dataset.getSpace();
            hsize_t dims[3];
            dataspace.getSimpleExtentDims(dims);
            std::vector<double> data(dims[0] * dims[1] * dims[2]);
            dataset.read(data.data(), H5::PredType::NATIVE_DOUBLE);
            return data;
        };
        
        // Load coordinate axes
        plasma_model_.x_axis = read_1d("coordinates/x");
        plasma_model_.y_axis = read_1d("coordinates/y");
        plasma_model_.z_axis = read_1d("coordinates/z");
        
        plasma_model_.nx = plasma_model_.x_axis.size();
        plasma_model_.ny = plasma_model_.y_axis.size();
        plasma_model_.nz = plasma_model_.z_axis.size();
        
        // Load electron density fields [cm^-3]
        plasma_model_.nec = read_3d("data/ne_c");
        plasma_model_.neh = read_3d("data/ne_h");
        
        // Load ion density fields [cm^-3]
        plasma_model_.nsp = read_3d("data/nsp");
        plasma_model_.ns2p = read_3d("data/ns2p");
        plasma_model_.ns3p = read_3d("data/ns3p");
        plasma_model_.nop = read_3d("data/nop");
        plasma_model_.no2p = read_3d("data/no2p");
        
        // Load temperature fields [eV]
        plasma_model_.Tec = read_3d("data/Te_c");
        plasma_model_.Teh = read_3d("data/Te_h");
        
        file.close();
        
        // Calculate derived quantities and clean data
        size_t total_size = plasma_model_.nx * plasma_model_.ny * plasma_model_.nz;
        plasma_model_.ne_total.resize(total_size);
        plasma_model_.feh.resize(total_size);
        
        // Helper to clean a single value - replace NaN/Inf/negative with 0
        auto clean_value = [](double val) -> double {
            if (!is_valid_double(val) || val < 0.0) return 0.0;
            return val;
        };
        
        for (size_t i = 0; i < total_size; ++i) {
            // Clean electron densities
            plasma_model_.nec[i] = clean_value(plasma_model_.nec[i]);
            plasma_model_.neh[i] = clean_value(plasma_model_.neh[i]);
            
            // Calculate total and hot fraction
            double ne_total = plasma_model_.nec[i] + plasma_model_.neh[i];
            plasma_model_.ne_total[i] = ne_total;
            
            // Hot fraction with safe division
            if (ne_total > 1.0e-10) {
                plasma_model_.feh[i] = plasma_model_.neh[i] / ne_total;
            } else {
                plasma_model_.feh[i] = 0.0;
            }
            
            // Clean ion densities
            plasma_model_.nsp[i] = clean_value(plasma_model_.nsp[i]);
            plasma_model_.ns2p[i] = clean_value(plasma_model_.ns2p[i]);
            plasma_model_.ns3p[i] = clean_value(plasma_model_.ns3p[i]);
            plasma_model_.nop[i] = clean_value(plasma_model_.nop[i]);
            plasma_model_.no2p[i] = clean_value(plasma_model_.no2p[i]);
            
            // Clean temperatures (must be positive for valid plasma)
            plasma_model_.Tec[i] = clean_value(plasma_model_.Tec[i]);
            plasma_model_.Teh[i] = clean_value(plasma_model_.Teh[i]);
        }
        
        std::cout << "Plasma model loaded: grid shape " 
                  << plasma_model_.nx << "x" << plasma_model_.ny << "x" << plasma_model_.nz << std::endl;
        std::cout << "  X range: " << plasma_model_.x_axis.front() << " to " 
                  << plasma_model_.x_axis.back() << " R_J" << std::endl;
        std::cout << "  Y range: " << plasma_model_.y_axis.front() << " to " 
                  << plasma_model_.y_axis.back() << " R_J" << std::endl;
        std::cout << "  Z range: " << plasma_model_.z_axis.front() << " to " 
                  << plasma_model_.z_axis.back() << " R_J" << std::endl;
        
    } catch (const H5::Exception& e) {
        throw std::runtime_error("Failed to load plasma model: " + std::string(e.getCDetailMsg()));
    }
}

void JovianUVEmissionRaytracer::load_emission_tables_single(const std::string& filename) {
    std::cout << "Loading single Maxwellian emission tables from " << filename << "..." << std::endl;
    
    try {
        H5::H5File file(filename, H5F_ACC_RDONLY);
        
        // Load temperature and density grids
        {
            H5::DataSet ds_T = file.openDataSet("T");
            H5::DataSpace space_T = ds_T.getSpace();
            hsize_t dim_T;
            space_T.getSimpleExtentDims(&dim_T);
            temp_arr_.resize(dim_T);
            ds_T.read(temp_arr_.data(), H5::PredType::NATIVE_DOUBLE);
        }
        
        {
            H5::DataSet ds_n = file.openDataSet("n");
            H5::DataSpace space_n = ds_n.getSpace();
            hsize_t dim_n;
            space_n.getSimpleExtentDims(&dim_n);
            dens_arr_.resize(dim_n);
            ds_n.read(dens_arr_.data(), H5::PredType::NATIVE_DOUBLE);
        }
        
        // Create log-space grids for interpolation
        log_temp_.resize(temp_arr_.size());
        log_dens_.resize(dens_arr_.size());
        for (size_t i = 0; i < temp_arr_.size(); ++i) {
            log_temp_[i] = std::log10(temp_arr_[i]);
        }
        for (size_t i = 0; i < dens_arr_.size(); ++i) {
            log_dens_[i] = std::log10(dens_arr_[i]);
        }
        
        // Load emission array: (n_T, n_n, n_lines)
        H5::DataSet ds_emiss = file.openDataSet("emiss");
        H5::DataSpace space_emiss = ds_emiss.getSpace();
        hsize_t dims_emiss[3];
        space_emiss.getSimpleExtentDims(dims_emiss);
        
        size_t n_T = dims_emiss[0];
        size_t n_n = dims_emiss[1];
        size_t n_lines_total = dims_emiss[2];
        
        std::vector<double> emissivity_all(n_T * n_n * n_lines_total);
        ds_emiss.read(emissivity_all.data(), H5::PredType::NATIVE_DOUBLE);
        
        // Clean emission data - replace NaN/Inf/negative with 0
        for (auto& val : emissivity_all) {
            if (!is_valid_double(val) || val < 0.0) val = 0.0;
        }
        
        // Load wavelengths
        H5::DataSet ds_wav = file.openDataSet("wavelength");
        H5::DataSpace space_wav = ds_wav.getSpace();
        hsize_t dim_wav;
        space_wav.getSimpleExtentDims(&dim_wav);
        std::vector<double> wavelength_all(dim_wav);
        ds_wav.read(wavelength_all.data(), H5::PredType::NATIVE_DOUBLE);
        
        // Load species labels
        H5::DataSet ds_species = file.openDataSet("species");
        H5::DataSpace space_species = ds_species.getSpace();
        hsize_t dim_species;
        space_species.getSimpleExtentDims(&dim_species);
        
        H5::DataType dtype_species = ds_species.getDataType();
        size_t str_size = dtype_species.getSize();
        std::vector<char> species_buffer(dim_species * str_size);
        ds_species.read(species_buffer.data(), dtype_species);
        
        std::vector<std::string> species_all;
        for (size_t i = 0; i < dim_species; ++i) {
            std::string s(species_buffer.data() + i * str_size, str_size);
            // Trim whitespace and null characters
            s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
                return !std::isspace(ch) && ch != '\0';
            }).base(), s.end());
            species_all.push_back(s);
        }
        
        file.close();
        
        // Get unique species
        std::vector<std::string> unique_species;
        for (const auto& s : species_all) {
            if (std::find(unique_species.begin(), unique_species.end(), s) == unique_species.end()) {
                unique_species.push_back(s);
            }
        }
        
        // Organize by species and transpose to (n_lines, n_T, n_n)
        for (const auto& species_key : unique_species) {
            std::vector<size_t> indices;
            for (size_t i = 0; i < species_all.size(); ++i) {
                if (species_all[i] == species_key) {
                    indices.push_back(i);
                }
            }
            
            EmissionTableSingle table;
            table.n_lines = indices.size();
            table.n_T = n_T;
            table.n_n = n_n;
            table.wavelengths.resize(indices.size());
            table.emissivities.resize(indices.size() * n_T * n_n, 0.0);
            
            for (size_t li = 0; li < indices.size(); ++li) {
                size_t orig_idx = indices[li];
                table.wavelengths[li] = wavelength_all[orig_idx];
                
                // Transpose: original (n_T, n_n, n_lines) -> target (n_lines, n_T, n_n)
                for (size_t ti = 0; ti < n_T; ++ti) {
                    for (size_t ni = 0; ni < n_n; ++ni) {
                        size_t src_idx = ti * n_n * n_lines_total + ni * n_lines_total + orig_idx;
                        size_t dst_idx = li * n_T * n_n + ti * n_n + ni;
                        table.emissivities[dst_idx] = emissivity_all[src_idx];
                    }
                }
            }
            
            emiss_single_[species_key] = std::move(table);
        }
        
        single_maxwellian_loaded_ = true;
        
        std::cout << "Single Maxwellian tables loaded successfully:" << std::endl;
        std::cout << "  Temperature range: " << temp_arr_.front() << " - " 
                  << temp_arr_.back() << " eV" << std::endl;
        std::cout << "  Density range: " << dens_arr_.front() << " - " 
                  << dens_arr_.back() << " cm^-3" << std::endl;
        std::cout << "  Grid size: " << temp_arr_.size() << " x " << dens_arr_.size() << std::endl;
        std::cout << "  Species: ";
        for (size_t i = 0; i < unique_species.size(); ++i) {
            std::cout << unique_species[i];
            if (i < unique_species.size() - 1) std::cout << ", ";
        }
        std::cout << std::endl;
        
    } catch (const H5::Exception& e) {
        throw std::runtime_error("Failed to load single Maxwellian tables: " + 
                                 std::string(e.getCDetailMsg()));
    }
}

void JovianUVEmissionRaytracer::load_emission_tables_double(const std::string& filename) {
    if (filename.empty()) {
        double_maxwellian_loaded_ = false;
        return;
    }
    
    std::cout << "Loading double Maxwellian emission tables from " << filename << "..." << std::endl;
    
    try {
        H5::H5File file(filename, H5F_ACC_RDONLY);
        
        auto read_1d = [&file](const std::string& path) -> std::vector<double> {
            H5::DataSet ds = file.openDataSet(path);
            H5::DataSpace space = ds.getSpace();
            hsize_t dim;
            space.getSimpleExtentDims(&dim);
            std::vector<double> data(dim);
            ds.read(data.data(), H5::PredType::NATIVE_DOUBLE);
            return data;
        };
        
        // Load parameter grids
        tec_arr_ = read_1d("T_cold");
        teh_arr_ = read_1d("T_hot");
        ne_arr_ = read_1d("n");
        feh_arr_ = read_1d("feh");
        
        // Create log-space grids
        log_tec_.resize(tec_arr_.size());
        log_teh_.resize(teh_arr_.size());
        log_ne_.resize(ne_arr_.size());
        log_feh_.resize(feh_arr_.size());
        
        for (size_t i = 0; i < tec_arr_.size(); ++i) log_tec_[i] = std::log10(tec_arr_[i]);
        for (size_t i = 0; i < teh_arr_.size(); ++i) log_teh_[i] = std::log10(teh_arr_[i]);
        for (size_t i = 0; i < ne_arr_.size(); ++i) log_ne_[i] = std::log10(ne_arr_[i]);
        for (size_t i = 0; i < feh_arr_.size(); ++i) log_feh_[i] = std::log10(feh_arr_[i]);
        
        // Load emission array: (n_Tc, n_Th, n_n, n_feh, n_lines)
        H5::DataSet ds_emiss = file.openDataSet("emiss");
        H5::DataSpace space_emiss = ds_emiss.getSpace();
        hsize_t dims_emiss[5];
        space_emiss.getSimpleExtentDims(dims_emiss);
        
        size_t n_Tec = dims_emiss[0];
        size_t n_Teh = dims_emiss[1];
        size_t n_ne = dims_emiss[2];
        size_t n_feh = dims_emiss[3];
        size_t n_lines_total = dims_emiss[4];
        
        std::vector<double> emissivity_all(n_Tec * n_Teh * n_ne * n_feh * n_lines_total);
        ds_emiss.read(emissivity_all.data(), H5::PredType::NATIVE_DOUBLE);
        
        // Clean emission data
        for (auto& val : emissivity_all) {
            if (!is_valid_double(val) || val < 0.0) val = 0.0;
        }
        
        // Load wavelengths
        std::vector<double> wavelength_all = read_1d("wavelength");
        
        // Load species labels
        H5::DataSet ds_species = file.openDataSet("species");
        H5::DataSpace space_species = ds_species.getSpace();
        hsize_t dim_species;
        space_species.getSimpleExtentDims(&dim_species);
        
        H5::DataType dtype_species = ds_species.getDataType();
        size_t str_size = dtype_species.getSize();
        std::vector<char> species_buffer(dim_species * str_size);
        ds_species.read(species_buffer.data(), dtype_species);
        
        std::vector<std::string> species_all;
        for (size_t i = 0; i < dim_species; ++i) {
            std::string s(species_buffer.data() + i * str_size, str_size);
            s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
                return !std::isspace(ch) && ch != '\0';
            }).base(), s.end());
            species_all.push_back(s);
        }
        
        file.close();
        
        // Get unique species
        std::vector<std::string> unique_species;
        for (const auto& s : species_all) {
            if (std::find(unique_species.begin(), unique_species.end(), s) == unique_species.end()) {
                unique_species.push_back(s);
            }
        }
        
        // Source strides for (n_Tec, n_Teh, n_ne, n_feh, n_lines)
        size_t src_stride_line = 1;
        size_t src_stride_feh = n_lines_total;
        size_t src_stride_ne = n_feh * n_lines_total;
        size_t src_stride_Teh = n_ne * n_feh * n_lines_total;
        size_t src_stride_Tec = n_Teh * n_ne * n_feh * n_lines_total;
        
        // Organize by species and transpose to (n_lines, n_Tec, n_Teh, n_ne, n_feh)
        for (const auto& species_key : unique_species) {
            std::vector<size_t> indices;
            for (size_t i = 0; i < species_all.size(); ++i) {
                if (species_all[i] == species_key) {
                    indices.push_back(i);
                }
            }
            
            EmissionTableDouble table;
            table.n_lines = indices.size();
            table.n_Tec = n_Tec;
            table.n_Teh = n_Teh;
            table.n_ne = n_ne;
            table.n_feh = n_feh;
            table.wavelengths.resize(indices.size());
            table.emissivities.resize(indices.size() * n_Tec * n_Teh * n_ne * n_feh, 0.0);
            
            // Dest strides for (n_lines, n_Tec, n_Teh, n_ne, n_feh)
            size_t dst_stride_feh = 1;
            size_t dst_stride_ne = n_feh;
            size_t dst_stride_Teh = n_ne * n_feh;
            size_t dst_stride_Tec = n_Teh * n_ne * n_feh;
            size_t dst_stride_line = n_Tec * n_Teh * n_ne * n_feh;
            
            for (size_t li = 0; li < indices.size(); ++li) {
                size_t orig_idx = indices[li];
                table.wavelengths[li] = wavelength_all[orig_idx];
                
                for (size_t i_Tec = 0; i_Tec < n_Tec; ++i_Tec) {
                    for (size_t i_Teh = 0; i_Teh < n_Teh; ++i_Teh) {
                        for (size_t i_ne = 0; i_ne < n_ne; ++i_ne) {
                            for (size_t i_feh = 0; i_feh < n_feh; ++i_feh) {
                                size_t src_idx = i_Tec * src_stride_Tec + 
                                                 i_Teh * src_stride_Teh +
                                                 i_ne * src_stride_ne +
                                                 i_feh * src_stride_feh +
                                                 orig_idx * src_stride_line;
                                size_t dst_idx = li * dst_stride_line +
                                                 i_Tec * dst_stride_Tec +
                                                 i_Teh * dst_stride_Teh +
                                                 i_ne * dst_stride_ne +
                                                 i_feh * dst_stride_feh;
                                table.emissivities[dst_idx] = emissivity_all[src_idx];
                            }
                        }
                    }
                }
            }
            
            emiss_double_[species_key] = std::move(table);
        }
        
        double_maxwellian_loaded_ = true;
        
        std::cout << "Double Maxwellian tables loaded successfully:" << std::endl;
        std::cout << "  Core temperature range: " << tec_arr_.front() << " - " 
                  << tec_arr_.back() << " eV" << std::endl;
        std::cout << "  Hot temperature range: " << teh_arr_.front() << " - " 
                  << teh_arr_.back() << " eV" << std::endl;
        std::cout << "  Density range: " << ne_arr_.front() << " - " 
                  << ne_arr_.back() << " cm^-3" << std::endl;
        std::cout << "  Hot fraction range: " << std::scientific << std::setprecision(6)
                  << feh_arr_.front() << " - " << std::fixed << std::setprecision(4)
                  << feh_arr_.back() << std::endl;
        std::cout << std::fixed << std::setprecision(1);
        std::cout << "  Grid size: " << tec_arr_.size() << " x " << teh_arr_.size() 
                  << " x " << ne_arr_.size() << " x " << feh_arr_.size() << std::endl;
        std::cout << "  Species: ";
        for (size_t i = 0; i < unique_species.size(); ++i) {
            std::cout << unique_species[i];
            if (i < unique_species.size() - 1) std::cout << ", ";
        }
        std::cout << std::endl;
        
    } catch (const H5::Exception& e) {
        throw std::runtime_error("Failed to load double Maxwellian tables: " + 
                                 std::string(e.getCDetailMsg()));
    }
}

// =============================================================================
// RAY TRACING
// =============================================================================

std::tuple<std::vector<double>, std::vector<std::array<double, 3>>>
JovianUVEmissionRaytracer::trace_ray(
    const std::array<double, 3>& start_pos,
    const std::array<double, 3>& direction,
    double ds,
    double max_distance) const {
    
    // Normalize direction vector
    double norm = std::sqrt(direction[0]*direction[0] + 
                           direction[1]*direction[1] + 
                           direction[2]*direction[2]);
    if (norm < 1e-15) norm = 1.0;
    
    std::array<double, 3> dir_norm = {
        direction[0] / norm,
        direction[1] / norm,
        direction[2] / norm
    };
    
    // Create ray points with uniform spacing
    size_t n_points = static_cast<size_t>(max_distance / ds) + 1;
    
    std::vector<double> s_values(n_points);
    std::vector<std::array<double, 3>> positions(n_points);
    
    for (size_t i = 0; i < n_points; ++i) {
        double s = static_cast<double>(i) * ds;
        s_values[i] = s;
        positions[i] = {
            start_pos[0] + s * dir_norm[0],
            start_pos[1] + s * dir_norm[1],
            start_pos[2] + s * dir_norm[2]
        };
    }
    
    return {s_values, positions};
}

// =============================================================================
// INTERPOLATION METHODS
// =============================================================================

double JovianUVEmissionRaytracer::interpolate_3d(
    const std::vector<double>& field,
    double x, double y, double z) const {
    
    const auto& x_ax = plasma_model_.x_axis;
    const auto& y_ax = plasma_model_.y_axis;
    const auto& z_ax = plasma_model_.z_axis;
    
    // Check bounds - return 0 for points outside the grid
    if (x < x_ax.front() || x > x_ax.back() ||
        y < y_ax.front() || y > y_ax.back() ||
        z < z_ax.front() || z > z_ax.back()) {
        return 0.0;
    }
    
    // Find bracketing indices
    size_t ix = searchsorted(x_ax, x);
    size_t iy = searchsorted(y_ax, y);
    size_t iz = searchsorted(z_ax, z);
    
    size_t ix0 = ix - 1, ix1 = ix;
    size_t iy0 = iy - 1, iy1 = iy;
    size_t iz0 = iz - 1, iz1 = iz;
    
    // Compute interpolation weights with safe division
    double dx = x_ax[ix1] - x_ax[ix0];
    double dy = y_ax[iy1] - y_ax[iy0];
    double dz = z_ax[iz1] - z_ax[iz0];
    
    double wx = (std::fabs(dx) > 1e-15) ? (x - x_ax[ix0]) / dx : 0.0;
    double wy = (std::fabs(dy) > 1e-15) ? (y - y_ax[iy0]) / dy : 0.0;
    double wz = (std::fabs(dz) > 1e-15) ? (z - z_ax[iz0]) / dz : 0.0;
    
    // Clamp weights to [0,1]
    wx = std::max(0.0, std::min(1.0, wx));
    wy = std::max(0.0, std::min(1.0, wy));
    wz = std::max(0.0, std::min(1.0, wz));
    
    // Get corner values with bounds checking
    auto safe_get = [&](size_t i, size_t j, size_t k) -> double {
        if (i >= plasma_model_.nx || j >= plasma_model_.ny || k >= plasma_model_.nz) {
            return 0.0;
        }
        double val = field[plasma_model_.index(i, j, k)];
        return is_valid_double(val) ? val : 0.0;
    };
    
    // Trilinear interpolation
    double c000 = safe_get(ix0, iy0, iz0);
    double c001 = safe_get(ix0, iy0, iz1);
    double c010 = safe_get(ix0, iy1, iz0);
    double c011 = safe_get(ix0, iy1, iz1);
    double c100 = safe_get(ix1, iy0, iz0);
    double c101 = safe_get(ix1, iy0, iz1);
    double c110 = safe_get(ix1, iy1, iz0);
    double c111 = safe_get(ix1, iy1, iz1);
    
    double c00 = c000 * (1.0 - wx) + c100 * wx;
    double c01 = c001 * (1.0 - wx) + c101 * wx;
    double c10 = c010 * (1.0 - wx) + c110 * wx;
    double c11 = c011 * (1.0 - wx) + c111 * wx;
    
    double c0 = c00 * (1.0 - wy) + c10 * wy;
    double c1 = c01 * (1.0 - wy) + c11 * wy;
    
    double result = c0 * (1.0 - wz) + c1 * wz;
    
    return is_valid_double(result) ? result : 0.0;
}

std::vector<double> JovianUVEmissionRaytracer::interpolate_3d_vectorized(
    const std::vector<double>& field,
    const std::vector<std::array<double, 3>>& positions) const {
    
    std::vector<double> result(positions.size());
    for (size_t i = 0; i < positions.size(); ++i) {
        result[i] = interpolate_3d(field, positions[i][0], positions[i][1], positions[i][2]);
    }
    return result;
}

const std::vector<double>& JovianUVEmissionRaytracer::get_ion_field(
    const std::string& species_key) const {
    
    if (species_key == "SP") return plasma_model_.nsp;
    if (species_key == "S2P") return plasma_model_.ns2p;
    if (species_key == "S3P") return plasma_model_.ns3p;
    if (species_key == "OP") return plasma_model_.nop;
    if (species_key == "O2P") return plasma_model_.no2p;
    
    throw std::invalid_argument("Unknown species key: " + species_key);
}

std::tuple<std::vector<double>, std::vector<std::vector<double>>>
JovianUVEmissionRaytracer::interpolate_emission_rates_single_vectorized(
    const std::vector<double>& Te_arr,
    const std::vector<double>& ne_arr,
    const std::string& species_key,
    double min_wav,
    double max_wav) const {
    
    auto it = emiss_single_.find(species_key);
    if (it == emiss_single_.end()) {
        return {{}, {}};
    }
    
    const EmissionTableSingle& table = it->second;
    
    // Filter wavelengths to range
    std::vector<size_t> wav_indices;
    std::vector<double> wavelengths;
    for (size_t i = 0; i < table.wavelengths.size(); ++i) {
        double wl = table.wavelengths[i];
        if (is_valid_double(wl) && wl >= min_wav && wl <= max_wav) {
            wav_indices.push_back(i);
            wavelengths.push_back(wl);
        }
    }
    
    if (wavelengths.empty()) {
        return {{}, {}};
    }
    
    size_t n_los = Te_arr.size();
    size_t n_lines = wavelengths.size();
    
    // Initialize output with zeros
    std::vector<std::vector<double>> emission_rates(n_los, std::vector<double>(n_lines, 0.0));
    
    // Get bounds from emission table
    double T_min = temp_arr_.front();
    double T_max = temp_arr_.back();
    double n_min = dens_arr_.front();
    double n_max = dens_arr_.back();
    
    // Process each LOS point
    for (size_t los = 0; los < n_los; ++los) {
        double Te = Te_arr[los];
        double ne = ne_arr[los];
        
        // Skip invalid or out-of-range values
        if (!is_valid_double(Te) || !is_valid_double(ne)) continue;
        if (Te < T_min || Te > T_max) continue;
        if (ne < n_min || ne > n_max) continue;
        
        // Convert to log space
        double log_T = std::log10(Te);
        double log_n = std::log10(ne);
        
        // Find bracketing indices
        size_t i_T = searchsorted(log_temp_, log_T);
        size_t i_n = searchsorted(log_dens_, log_n);
        
        size_t i_T0 = i_T - 1, i_T1 = i_T;
        size_t i_n0 = i_n - 1, i_n1 = i_n;
        
        // Compute interpolation weights
        double dT = log_temp_[i_T1] - log_temp_[i_T0];
        double dn = log_dens_[i_n1] - log_dens_[i_n0];
        
        double w_T = (std::fabs(dT) > 1e-15) ? (log_T - log_temp_[i_T0]) / dT : 0.0;
        double w_n = (std::fabs(dn) > 1e-15) ? (log_n - log_dens_[i_n0]) / dn : 0.0;
        
        // Clamp weights to [0,1]
        w_T = std::max(0.0, std::min(1.0, w_T));
        w_n = std::max(0.0, std::min(1.0, w_n));
        
        // Bilinear interpolation weights
        double w00 = (1.0 - w_T) * (1.0 - w_n);
        double w01 = (1.0 - w_T) * w_n;
        double w10 = w_T * (1.0 - w_n);
        double w11 = w_T * w_n;
        
        // Interpolate for each line
        for (size_t li = 0; li < n_lines; ++li) {
            size_t orig_idx = wav_indices[li];
            
            double E_00 = table.get(orig_idx, i_T0, i_n0);
            double E_01 = table.get(orig_idx, i_T0, i_n1);
            double E_10 = table.get(orig_idx, i_T1, i_n0);
            double E_11 = table.get(orig_idx, i_T1, i_n1);
            
            // Check for valid corner values
            if (!is_valid_double(E_00)) E_00 = 0.0;
            if (!is_valid_double(E_01)) E_01 = 0.0;
            if (!is_valid_double(E_10)) E_10 = 0.0;
            if (!is_valid_double(E_11)) E_11 = 0.0;
            
            double result = w00 * E_00 + w01 * E_01 + w10 * E_10 + w11 * E_11;
            emission_rates[los][li] = is_valid_double(result) ? result : 0.0;
        }
    }
    
    return {wavelengths, emission_rates};
}

std::tuple<std::vector<double>, std::vector<std::vector<double>>>
JovianUVEmissionRaytracer::interpolate_emission_rates_double_vectorized(
    const std::vector<double>& Tec_arr,
    const std::vector<double>& Teh_arr,
    const std::vector<double>& ne_arr,
    const std::vector<double>& feh_arr,
    const std::string& species_key,
    double min_wav,
    double max_wav) const {
    
    if (!double_maxwellian_loaded_) {
        return interpolate_emission_rates_single_vectorized(Tec_arr, ne_arr, species_key, 
                                                            min_wav, max_wav);
    }
    
    auto it = emiss_double_.find(species_key);
    if (it == emiss_double_.end()) {
        return {{}, {}};
    }
    
    const EmissionTableDouble& table = it->second;
    
    // Filter wavelengths
    std::vector<size_t> wav_indices;
    std::vector<double> wavelengths;
    for (size_t i = 0; i < table.wavelengths.size(); ++i) {
        double wl = table.wavelengths[i];
        if (is_valid_double(wl) && wl >= min_wav && wl <= max_wav) {
            wav_indices.push_back(i);
            wavelengths.push_back(wl);
        }
    }
    
    if (wavelengths.empty()) {
        return {{}, {}};
    }
    
    size_t n_los = Tec_arr.size();
    size_t n_lines = wavelengths.size();
    
    std::vector<std::vector<double>> emission_rates(n_los, std::vector<double>(n_lines, 0.0));
    
    // Get bounds from double maxwellian table
    double Tec_min = tec_arr_.front();
    double Tec_max = tec_arr_.back();
    double Teh_min = teh_arr_.front();
    double Teh_max = teh_arr_.back();
    double ne_min_d = ne_arr_.front();
    double ne_max_d = ne_arr_.back();
    double feh_min = feh_arr_.front();
    double feh_max = feh_arr_.back();
    
    for (size_t los = 0; los < n_los; ++los) {
        double Tec = Tec_arr[los];
        double Teh = Teh_arr[los];
        double ne = ne_arr[los];
        double feh = feh_arr[los];
        
        // Check validity
        if (!is_valid_double(Tec) || !is_valid_double(Teh) || 
            !is_valid_double(ne) || !is_valid_double(feh)) {
            continue;
        }
        
        // Use single Maxwellian for out-of-range values
        bool out_of_range = (Tec < Tec_min || Tec > Tec_max ||
                            Teh < Teh_min || Teh > Teh_max ||
                            ne < ne_min_d || ne > ne_max_d ||
                            feh < feh_min || feh > feh_max);
        
        if (out_of_range) {
            // Fall back to single Maxwellian using core temperature
            auto [wav_s, emiss_s] = interpolate_emission_rates_single_vectorized(
                {Tec}, {ne}, species_key, min_wav, max_wav);
            if (!wav_s.empty() && !emiss_s.empty()) {
                emission_rates[los] = emiss_s[0];
            }
            continue;
        }
        
        // Convert to log space
        double log_Tec = std::log10(Tec);
        double log_Teh = std::log10(Teh);
        double log_ne = std::log10(ne);
        double log_feh = std::log10(feh);
        
        // Find bracketing indices
        size_t i_Tec = searchsorted(log_tec_, log_Tec);
        size_t i_Teh = searchsorted(log_teh_, log_Teh);
        size_t i_ne = searchsorted(log_ne_, log_ne);
        size_t i_feh = searchsorted(log_feh_, log_feh);
        
        size_t i0_Tec = i_Tec - 1, i1_Tec = i_Tec;
        size_t i0_Teh = i_Teh - 1, i1_Teh = i_Teh;
        size_t i0_ne = i_ne - 1, i1_ne = i_ne;
        size_t i0_feh = i_feh - 1, i1_feh = i_feh;
        
        // Compute interpolation weights
        double dTec = log_tec_[i1_Tec] - log_tec_[i0_Tec];
        double dTeh = log_teh_[i1_Teh] - log_teh_[i0_Teh];
        double dne = log_ne_[i1_ne] - log_ne_[i0_ne];
        double dfeh = log_feh_[i1_feh] - log_feh_[i0_feh];
        
        double w_Tec = (std::fabs(dTec) > 1e-15) ? (log_Tec - log_tec_[i0_Tec]) / dTec : 0.0;
        double w_Teh = (std::fabs(dTeh) > 1e-15) ? (log_Teh - log_teh_[i0_Teh]) / dTeh : 0.0;
        double w_ne = (std::fabs(dne) > 1e-15) ? (log_ne - log_ne_[i0_ne]) / dne : 0.0;
        double w_feh = (std::fabs(dfeh) > 1e-15) ? (log_feh - log_feh_[i0_feh]) / dfeh : 0.0;
        
        // Clamp weights
        w_Tec = std::max(0.0, std::min(1.0, w_Tec));
        w_Teh = std::max(0.0, std::min(1.0, w_Teh));
        w_ne = std::max(0.0, std::min(1.0, w_ne));
        w_feh = std::max(0.0, std::min(1.0, w_feh));
        
        // Quadrilinear interpolation for each line
        for (size_t li = 0; li < n_lines; ++li) {
            size_t orig_idx = wav_indices[li];
            double emiss_val = 0.0;
            
            // Loop over 16 corners of 4D hypercube
            for (int b_Tec = 0; b_Tec < 2; ++b_Tec) {
                for (int b_Teh = 0; b_Teh < 2; ++b_Teh) {
                    for (int b_ne = 0; b_ne < 2; ++b_ne) {
                        for (int b_feh = 0; b_feh < 2; ++b_feh) {
                            size_t idx_Tec = (b_Tec == 0) ? i0_Tec : i1_Tec;
                            size_t idx_Teh = (b_Teh == 0) ? i0_Teh : i1_Teh;
                            size_t idx_ne = (b_ne == 0) ? i0_ne : i1_ne;
                            size_t idx_feh = (b_feh == 0) ? i0_feh : i1_feh;
                            
                            double wt_Tec = (b_Tec == 0) ? (1.0 - w_Tec) : w_Tec;
                            double wt_Teh = (b_Teh == 0) ? (1.0 - w_Teh) : w_Teh;
                            double wt_ne = (b_ne == 0) ? (1.0 - w_ne) : w_ne;
                            double wt_feh = (b_feh == 0) ? (1.0 - w_feh) : w_feh;
                            
                            double weight = wt_Tec * wt_Teh * wt_ne * wt_feh;
                            double corner_val = table.get(orig_idx, idx_Tec, idx_Teh, idx_ne, idx_feh);
                            
                            if (is_valid_double(corner_val)) {
                                emiss_val += weight * corner_val;
                            }
                        }
                    }
                }
            }
            
            emission_rates[los][li] = is_valid_double(emiss_val) ? emiss_val : 0.0;
        }
    }
    
    return {wavelengths, emission_rates};
}

// =============================================================================
// LINE-OF-SIGHT INTEGRATION
// =============================================================================

std::tuple<std::vector<double>, std::vector<double>>
JovianUVEmissionRaytracer::integrate_species_emission_single(
    const std::vector<double>& s_values,
    const std::string& species_key,
    const std::vector<double>& Te_los,
    const std::vector<double>& ne_los,
    const std::vector<double>& n_ion_los,
    double min_wav,
    double max_wav) const {
    
    auto [wavelengths, emission_rates] = interpolate_emission_rates_single_vectorized(
        Te_los, ne_los, species_key, min_wav, max_wav);
    
    if (wavelengths.empty()) {
        return {{}, {}};
    }
    
    size_t n_los = s_values.size();
    size_t n_lines = wavelengths.size();
    
    // Calculate step size
    double ds = (n_los > 1) ? (s_values[1] - s_values[0]) : 0.1;
    
    // Integrate each line along LOS
    std::vector<double> brightnesses(n_lines, 0.0);
    
    for (size_t li = 0; li < n_lines; ++li) {
        // Build integrand for this line: emission_rate * n_ion
        std::vector<double> integrand(n_los);
        for (size_t i = 0; i < n_los; ++i) {
            double emiss = emission_rates[i][li];
            double n_ion = n_ion_los[i];
            
            // Only include valid, positive values
            if (is_valid_double(emiss) && is_valid_double(n_ion) && emiss > 0 && n_ion > 0) {
                integrand[i] = emiss * n_ion;
            } else {
                integrand[i] = 0.0;
            }
        }
        
        // Integrate using Simpson's rule
        double integral = simpson_integrate_uniform(integrand, ds);
        
        // Convert to Rayleighs
        double brightness = RAYLEIGH_FACTOR * integral;
        brightnesses[li] = is_valid_double(brightness) ? brightness : 0.0;
    }
    
    return {wavelengths, brightnesses};
}

std::tuple<std::vector<double>, std::vector<double>>
JovianUVEmissionRaytracer::integrate_species_emission_double(
    const std::vector<double>& s_values,
    const std::string& species_key,
    const std::vector<double>& Tec_los,
    const std::vector<double>& Teh_los,
    const std::vector<double>& ne_los,
    const std::vector<double>& feh_los,
    const std::vector<double>& n_ion_los,
    double min_wav,
    double max_wav) const {
    
    auto [wavelengths, emission_rates] = interpolate_emission_rates_double_vectorized(
        Tec_los, Teh_los, ne_los, feh_los, species_key, min_wav, max_wav);
    
    if (wavelengths.empty()) {
        return {{}, {}};
    }
    
    size_t n_los = s_values.size();
    size_t n_lines = wavelengths.size();
    
    double ds = (n_los > 1) ? (s_values[1] - s_values[0]) : 0.1;
    
    std::vector<double> brightnesses(n_lines, 0.0);
    
    for (size_t li = 0; li < n_lines; ++li) {
        std::vector<double> integrand(n_los);
        for (size_t i = 0; i < n_los; ++i) {
            double emiss = emission_rates[i][li];
            double n_ion = n_ion_los[i];
            
            if (is_valid_double(emiss) && is_valid_double(n_ion) && emiss > 0 && n_ion > 0) {
                integrand[i] = emiss * n_ion;
            } else {
                integrand[i] = 0.0;
            }
        }
        
        double integral = simpson_integrate_uniform(integrand, ds);
        double brightness = RAYLEIGH_FACTOR * integral;
        brightnesses[li] = is_valid_double(brightness) ? brightness : 0.0;
    }
    
    return {wavelengths, brightnesses};
}

// =============================================================================
// INSTRUMENT RESPONSE CONVOLUTION
// =============================================================================

std::vector<double> JovianUVEmissionRaytracer::convolve_spectrum_erf(
    const std::vector<double>& wavelength_grid,
    double bin_width,
    const std::vector<double>& line_wavelengths,
    const std::vector<double>& line_brightnesses,
    double fwhm) const {
    
    size_t n_bins = wavelength_grid.size();
    std::vector<double> spectrum(n_bins, 0.0);
    
    if (line_wavelengths.empty() || line_brightnesses.empty() || fwhm <= 0 || bin_width <= 0) {
        return spectrum;
    }
    
    // Convert FWHM to Gaussian sigma
    double sigma = fwhm * FWHM_TO_SIGMA;
    double sigma_sqrt2 = sigma * SQRT2;
    
    if (sigma_sqrt2 < 1e-15) {
        return spectrum;
    }
    
    size_t n_lines = line_wavelengths.size();
    
    // For each bin, calculate contribution from all lines using ERF-based integration
    for (size_t bi = 0; bi < n_bins; ++bi) {
        double bin_center = wavelength_grid[bi];
        double bin_lo = bin_center - bin_width / 2.0;
        double bin_hi = bin_center + bin_width / 2.0;
        
        double bin_sum = 0.0;
        for (size_t li = 0; li < n_lines; ++li) {
            double brightness = line_brightnesses[li];
            
            // Skip invalid or non-positive brightnesses
            if (!is_valid_double(brightness) || brightness <= 0.0) {
                continue;
            }
            
            double line_wav = line_wavelengths[li];
            if (!is_valid_double(line_wav)) {
                continue;
            }
            
            // ERF-based Gaussian integration:
            // I_bin = I_line * 0.5 * [erf((λ_hi - λ_0)/(σ√2)) - erf((λ_lo - λ_0)/(σ√2))]
            double erf_arg_lo = (bin_lo - line_wav) / sigma_sqrt2;
            double erf_arg_hi = (bin_hi - line_wav) / sigma_sqrt2;
            
            double bin_contrib = 0.5 * (erf_scalar(erf_arg_hi) - erf_scalar(erf_arg_lo));
            bin_sum += brightness * bin_contrib;
        }
        
        // Divide by bin width to get Rayleighs per Angstrom
        spectrum[bi] = bin_sum / bin_width;
    }
    
    return spectrum;
}

// =============================================================================
// MAIN SPECTRUM CALCULATION METHODS
// =============================================================================

SpectrumResult JovianUVEmissionRaytracer::calculate_spectrum_single(
    const std::array<double, 3>& slit_pos_vec,
    const std::array<double, 3>& norm_vec,
    std::pair<double, double> wavelength_range,
    double bin_width,
    double fwhm,
    double ds) const {
    
    std::cout << "Calculating single Maxwellian spectrum:" << std::endl;
    std::cout << "  LOS start: [" << std::fixed << std::setprecision(1) 
              << slit_pos_vec[0] << ", " << slit_pos_vec[1] << ", " 
              << slit_pos_vec[2] << "] R_J" << std::endl;
    std::cout << std::setprecision(2);
    std::cout << "  Direction: [" << norm_vec[0] << ", " << norm_vec[1] << ", " 
              << norm_vec[2] << "]" << std::endl;
    
    double min_wav = wavelength_range.first;
    double max_wav = wavelength_range.second;
    
    // Create output wavelength grid matching Python np.arange behavior
    size_t n_wav = static_cast<size_t>((max_wav - min_wav) / bin_width) + 1;
    std::vector<double> wave_bins(n_wav);
    for (size_t i = 0; i < n_wav; ++i) {
        wave_bins[i] = min_wav + static_cast<double>(i) * bin_width;
    }
    
    // Trace ray through plasma model
    auto [s_values, positions] = trace_ray(slit_pos_vec, norm_vec, ds);
    
    // Interpolate plasma parameters along LOS (use cold electron parameters for single Maxwellian)
    std::vector<double> ne_los = interpolate_3d_vectorized(plasma_model_.nec, positions);
    std::vector<double> Te_los = interpolate_3d_vectorized(plasma_model_.Tec, positions);
    
    // Collect all emission lines from all species
    std::vector<double> all_wavelengths;
    std::vector<double> all_brightnesses;
    std::vector<std::string> all_species;
    
    for (const auto& species_key : default_species_) {
        const std::vector<double>& ion_field = get_ion_field(species_key);
        std::vector<double> n_ion_los = interpolate_3d_vectorized(ion_field, positions);
        
        auto [wavelengths, brightnesses] = integrate_species_emission_single(
            s_values, species_key, Te_los, ne_los, n_ion_los, min_wav, max_wav);
        
        if (!wavelengths.empty()) {
            for (size_t i = 0; i < wavelengths.size(); ++i) {
                if (is_valid_double(brightnesses[i]) && brightnesses[i] > 0) {
                    all_wavelengths.push_back(wavelengths[i]);
                    all_brightnesses.push_back(brightnesses[i]);
                    all_species.push_back(species_key);
                }
            }
        }
    }
    
    // Create line list
    std::vector<EmissionLine> line_list;
    for (size_t i = 0; i < all_wavelengths.size(); ++i) {
        if (all_brightnesses[i] > 1e-10) {
            line_list.emplace_back(all_wavelengths[i], all_brightnesses[i], all_species[i]);
        }
    }
    
    std::cout << "  Processed " << all_wavelengths.size() << " lines, " 
              << line_list.size() << " with significant brightness" << std::endl;
    
    // Convolve with instrument response
    std::vector<double> spectrum = convolve_spectrum_erf(wave_bins, bin_width,
                                                          all_wavelengths, all_brightnesses, fwhm);
    
    // Calculate total brightness by integrating spectrum
    double total = simpson_integrate_uniform(spectrum, bin_width);
    if (!is_valid_double(total)) total = 0.0;
    
    return {wave_bins, spectrum, line_list, total};
}

SpectrumResult JovianUVEmissionRaytracer::calculate_spectrum_double(
    const std::array<double, 3>& slit_pos_vec,
    const std::array<double, 3>& norm_vec,
    std::pair<double, double> wavelength_range,
    double bin_width,
    double fwhm,
    double ds) const {
    
    if (!double_maxwellian_loaded_) {
        throw std::runtime_error("Double Maxwellian tables not loaded.");
    }
    
    std::cout << "Calculating double Maxwellian spectrum:" << std::endl;
    std::cout << "  LOS start: [" << std::fixed << std::setprecision(1) 
              << slit_pos_vec[0] << ", " << slit_pos_vec[1] << ", " 
              << slit_pos_vec[2] << "] R_J" << std::endl;
    std::cout << std::setprecision(2);
    std::cout << "  Direction: [" << norm_vec[0] << ", " << norm_vec[1] << ", " 
              << norm_vec[2] << "]" << std::endl;
    
    double min_wav = wavelength_range.first;
    double max_wav = wavelength_range.second;
    
    // Create wavelength grid
    size_t n_wav = static_cast<size_t>((max_wav - min_wav) / bin_width) + 1;
    std::vector<double> wave_bins(n_wav);
    for (size_t i = 0; i < n_wav; ++i) {
        wave_bins[i] = min_wav + static_cast<double>(i) * bin_width;
    }
    
    // Trace ray
    auto [s_values, positions] = trace_ray(slit_pos_vec, norm_vec, ds);
    
    // Interpolate plasma parameters (use total electron density for double Maxwellian)
    std::vector<double> ne_los = interpolate_3d_vectorized(plasma_model_.ne_total, positions);
    std::vector<double> feh_los = interpolate_3d_vectorized(plasma_model_.feh, positions);
    std::vector<double> Tec_los = interpolate_3d_vectorized(plasma_model_.Tec, positions);
    std::vector<double> Teh_los = interpolate_3d_vectorized(plasma_model_.Teh, positions);
    
    // Collect emission lines
    std::vector<double> all_wavelengths;
    std::vector<double> all_brightnesses;
    std::vector<std::string> all_species;
    
    for (const auto& species_key : default_species_) {
        const std::vector<double>& ion_field = get_ion_field(species_key);
        std::vector<double> n_ion_los = interpolate_3d_vectorized(ion_field, positions);
        
        auto [wavelengths, brightnesses] = integrate_species_emission_double(
            s_values, species_key, Tec_los, Teh_los, ne_los, feh_los,
            n_ion_los, min_wav, max_wav);
        
        if (!wavelengths.empty()) {
            for (size_t i = 0; i < wavelengths.size(); ++i) {
                if (is_valid_double(brightnesses[i]) && brightnesses[i] > 0) {
                    all_wavelengths.push_back(wavelengths[i]);
                    all_brightnesses.push_back(brightnesses[i]);
                    all_species.push_back(species_key);
                }
            }
        }
    }
    
    // Create line list
    std::vector<EmissionLine> line_list;
    for (size_t i = 0; i < all_wavelengths.size(); ++i) {
        if (all_brightnesses[i] > 1e-10) {
            line_list.emplace_back(all_wavelengths[i], all_brightnesses[i], all_species[i]);
        }
    }
    
    std::cout << "  Processed " << all_wavelengths.size() << " lines, " 
              << line_list.size() << " with significant brightness" << std::endl;
    
    // Convolve
    std::vector<double> spectrum = convolve_spectrum_erf(wave_bins, bin_width,
                                                          all_wavelengths, all_brightnesses, fwhm);
    
    // Calculate total
    double total = simpson_integrate_uniform(spectrum, bin_width);
    if (!is_valid_double(total)) total = 0.0;
    
    return {wave_bins, spectrum, line_list, total};
}

// =============================================================================
// ACCESSORS
// =============================================================================

std::string JovianUVEmissionRaytracer::get_ion_name(const std::string& species_key) const {
    auto it = ion_names_.find(species_key);
    return (it != ion_names_.end()) ? it->second : species_key;
}

// =============================================================================
// STANDALONE HELPER FUNCTIONS
// =============================================================================

std::vector<double> simulate_ipt_spectrum_rayleighs_erf_form(
    const std::vector<double>& wavelength_grid,
    double bin_width,
    const std::vector<double>& line_wavelengths,
    const std::vector<double>& line_brightnesses,
    double fwhm) {
    
    size_t n_bins = wavelength_grid.size();
    std::vector<double> spectrum(n_bins, 0.0);
    
    if (line_wavelengths.empty() || line_brightnesses.empty() || fwhm <= 0 || bin_width <= 0) {
        return spectrum;
    }
    
    double sigma = fwhm * FWHM_TO_SIGMA;
    double sigma_sqrt2 = sigma * SQRT2;
    
    if (sigma_sqrt2 < 1e-15) return spectrum;
    
    size_t n_lines = line_wavelengths.size();
    
    for (size_t bi = 0; bi < n_bins; ++bi) {
        double bin_lo = wavelength_grid[bi] - bin_width / 2.0;
        double bin_hi = wavelength_grid[bi] + bin_width / 2.0;
        
        double bin_sum = 0.0;
        for (size_t li = 0; li < n_lines; ++li) {
            if (!is_valid_double(line_brightnesses[li]) || line_brightnesses[li] <= 0) continue;
            if (!is_valid_double(line_wavelengths[li])) continue;
            
            double erf_arg_lo = (bin_lo - line_wavelengths[li]) / sigma_sqrt2;
            double erf_arg_hi = (bin_hi - line_wavelengths[li]) / sigma_sqrt2;
            
            double bin_contrib = 0.5 * (erf_scalar(erf_arg_hi) - erf_scalar(erf_arg_lo));
            bin_sum += line_brightnesses[li] * bin_contrib;
        }
        
        spectrum[bi] = bin_sum / bin_width;
    }
    
    return spectrum;
}

} // namespace ipt
