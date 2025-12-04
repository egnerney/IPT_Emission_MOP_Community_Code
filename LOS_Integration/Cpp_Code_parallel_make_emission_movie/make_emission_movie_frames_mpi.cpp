/**
 * @file make_emission_movie_frames_mpi.cpp
 * @brief MPI-Parallelized Hybrid Optical+UV Emission Movie Frame Generator
 *
 * This program generates emission data for movie frames showing three specific
 * emission lines from the Jovian plasma torus in a stacked 3-panel view:
 *
 *     Top panel:    S+ 6731 Å (optical)
 *     Middle panel: S++ 680 Å (UV)  
 *     Bottom panel: O+ 833 Å (UV)
 *
 * @section parallelization MPI PARALLELIZATION
 * - Master process (rank 0) coordinates work distribution and I/O
 * - Worker processes compute emission for assigned grid points
 * - Each frame is parallelized over all available MPI ranks
 * - Shared memory avoided by having each rank load data independently
 *
 * @section output OUTPUT
 * Binary files in movie_data/ directory:
 *   - frame_XXXX.bin: Raw emission data for each frame
 *   - colorbar_limits.txt: Fixed colorbar limits from prescan
 *   - grid_info.txt: Grid dimensions and parameters
 *
 * A separate Python script (generate_frames_from_binary.py) reads these
 * files and generates PNG frames.
 *
 * @section movie_params MOVIE PARAMETERS
 * - Frame rate: 15 fps
 * - Duration: 5 seconds
 * - Total frames: 75 (covering 360° rotation)
 * - CML step: 4.8° per frame
 *
 * @section spectral_params SPECTRAL PARAMETERS
 * Each emission line uses optimized spectral windows:
 *
 *   S+ 6731 Å (Optical):
 *     - FWHM: 1.2 Å (narrow optical spectrometer)
 *     - Bin width: 0.6 Å
 *     - Window: ±10 Å → range 6721-6741 Å
 *
 *   S++ 680 Å (UV):
 *     - FWHM: 6.0 Å (broader UV)
 *     - Bin width: 1.0 Å
 *     - Window: ±10 Å → range 670-690 Å
 *
 *   O+ 833 Å (UV):
 *     - FWHM: 6.0 Å
 *     - Bin width: 1.0 Å
 *     - Window: ±10 Å → range 823-843 Å
 *
 * @section geometry VIEWING GEOMETRY
 * For each frame at CML angle θ (degrees):
 *   - CML = 0°:   Observer at X = -R_obs, looking in +X direction
 *   - CML = 90°:  Observer at Y = +R_obs, looking in -Y direction
 *   - CML = 180°: Observer at X = +R_obs, looking in -X direction
 *   - CML = 270°: Observer at Y = -R_obs, looking in +Y direction
 *
 * Ray direction: [cos(θ), -sin(θ), 0]
 * Ray starting position:
 *   start_x = -R_obs·cos(θ) - ρ·sin(θ)
 *   start_y =  R_obs·sin(θ) - ρ·cos(θ)
 *   start_z = z
 *
 * @section grid GRID CONFIGURATION
 * - ρ position: 221 points from -8 to +8 R_J (non-uniform sampling)
 * - z position: 131 points from -2 to +2 R_J (non-uniform sampling)
 * - Fine sampling near torus edges and equatorial plane
 * - Total per frame: 221 × 131 = 28,951 calculations
 *
 * @section colorbar COLORBAR LIMITS
 * - S+ 6731 Å: Automatic (prescan 4 CML angles with 6% headroom)
 * - S++ 680 Å: Fixed at 0-185 R
 * - O+ 833 Å: Fixed at 0-175 R
 *
 * @section usage USAGE
 * Compile:
 *   make
 *
 * Run with MPI:
 *   mpirun -np 8 ./make_emission_movie_frames_mpi
 *
 * Generate frames:
 *   python3 generate_frames_from_binary.py
 *
 * Create movie:
 *   python3 make_movies_ffmpeg.py
 *
 * @author Edward (Eddie) G. Nerney
 * @institution Laboratory for Atmospheric and Space Physics, University of Colorado Boulder
 * @date November 2025
 * @license Open source for academic and research use
 */

#include "../Cpp_Code/IPT_emiss_MOP_community_code.hpp"
#include <mpi.h>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <filesystem>
#include <sstream>

// =============================================================================
// CONFIGURATION CONSTANTS
// =============================================================================

namespace config {

// Movie parameters
constexpr int FRAME_RATE = 15;           // frames per second
constexpr int DURATION = 5;              // seconds
constexpr int N_FRAMES = FRAME_RATE * DURATION;  // 75 frames

// Ray tracing parameters
constexpr double R_OBS = 20.0;           // Observer distance [R_J]
constexpr double DS = 0.1;               // Integration step [R_J]

// S+ 6731 Å (Optical)
constexpr double S1P_6731_WAV = 6731.0;
constexpr double S1P_6731_FWHM = 1.2;
constexpr double S1P_6731_BIN_WIDTH = 0.6;
constexpr double S1P_6731_MIN = 6721.0;
constexpr double S1P_6731_MAX = 6741.0;

// S++ 680 Å (UV)
constexpr double S2P_680_WAV = 680.0;
constexpr double S2P_680_FWHM = 6.0;
constexpr double S2P_680_BIN_WIDTH = 1.0;
constexpr double S2P_680_MIN = 670.0;
constexpr double S2P_680_MAX = 690.0;

// O+ 833 Å (UV)
constexpr double OP_833_WAV = 833.0;
constexpr double OP_833_FWHM = 6.0;
constexpr double OP_833_BIN_WIDTH = 1.0;
constexpr double OP_833_MIN = 823.0;
constexpr double OP_833_MAX = 843.0;

// UV colorbar limits (fixed)
constexpr double S2P_680_VMAX = 185.0;
constexpr double OP_833_VMAX = 175.0;

// Prescan CML angles
const std::vector<double> PRESCAN_ANGLES = {0.0, 90.0, 180.0, 270.0};
constexpr double HEADROOM_FACTOR = 1.06;  // 6% headroom

// Output directory
const std::string OUTPUT_DIR = "movie_data";

}  // namespace config

// =============================================================================
// GRID GENERATION
// =============================================================================

/**
 * @brief Create a non-uniform grid from segment specifications.
 *
 * Matches the Python create_nonuniform_grid() function exactly.
 *
 * @param segments Vector of (start, stop, step) tuples
 * @return Concatenated grid points
 */
std::vector<double> create_nonuniform_grid(
    const std::vector<std::tuple<double, double, double>>& segments) {
    
    std::vector<double> grid;
    
    for (size_t seg_idx = 0; seg_idx < segments.size(); ++seg_idx) {
        double start = std::get<0>(segments[seg_idx]);
        double stop = std::get<1>(segments[seg_idx]);
        double step = std::get<2>(segments[seg_idx]);
        
        // For intermediate segments, don't include endpoint
        // For last segment, include endpoint
        double end_val = (seg_idx < segments.size() - 1) ? stop : stop + step / 2.0;
        
        for (double val = start; val < end_val; val += step) {
            grid.push_back(val);
        }
    }
    
    return grid;
}

/**
 * @brief Create the rho grid (radial positions).
 *
 * Non-uniform sampling with fine resolution near torus edges (4.5-6 R_J).
 *
 * @return Vector of rho positions [R_J], 221 points
 */
std::vector<double> create_rho_grid() {
    std::vector<std::tuple<double, double, double>> segments = {
        {-8.0, -6.0, 0.1},      // Outer edge
        {-6.0, -4.5, 0.025},    // Transition region (fine)
        {-4.5, -3.0, 0.1},      // Intermediate
        {-3.0, 3.0, 0.2},       // Central torus (coarse)
        {3.0, 4.5, 0.1},        // Intermediate
        {4.5, 6.0, 0.025},      // Transition region (fine)
        {6.0, 8.0, 0.1}         // Outer edge
    };
    return create_nonuniform_grid(segments);
}

/**
 * @brief Create the z grid (vertical positions).
 *
 * Non-uniform sampling with very fine resolution near equatorial plane.
 *
 * @return Vector of z positions [R_J], 131 points
 */
std::vector<double> create_z_grid() {
    std::vector<std::tuple<double, double, double>> segments = {
        {-2.0, -0.5, 0.1},      // Lower region
        {-0.5, 0.5, 0.01},      // Equatorial plane (very fine)
        {0.5, 2.0, 0.1}         // Upper region
    };
    return create_nonuniform_grid(segments);
}

// =============================================================================
// GEOMETRY FUNCTIONS
// =============================================================================

/**
 * @brief Compute ray starting position and direction for given sky coordinates and CML.
 *
 * Matches the Python compute_ray_geometry() function exactly.
 *
 * @param rho_sky Projected radial position in observer's sky plane [R_J]
 * @param z_sky Height above centrifugal equator [R_J]
 * @param cml_deg Central Meridian Longitude [degrees]
 * @param R_obs Observer distance from Jupiter center [R_J]
 * @param[out] slit_pos Ray starting position [x, y, z]
 * @param[out] norm_vec Ray direction vector (normalized)
 */
void compute_ray_geometry(double rho_sky, double z_sky, double cml_deg, double R_obs,
                          std::array<double, 3>& slit_pos, 
                          std::array<double, 3>& norm_vec) {
    double theta = cml_deg * M_PI / 180.0;
    
    // Viewing direction (toward Jupiter center)
    norm_vec[0] = std::cos(theta);
    norm_vec[1] = -std::sin(theta);
    norm_vec[2] = 0.0;
    
    // Ray starting position
    slit_pos[0] = -R_obs * std::cos(theta) - rho_sky * std::sin(theta);
    slit_pos[1] =  R_obs * std::sin(theta) - rho_sky * std::cos(theta);
    slit_pos[2] = z_sky;
}

// =============================================================================
// EMISSION CALCULATION
// =============================================================================

/**
 * @brief Calculate total emission for a single line of sight.
 *
 * Computes the integrated spectrum brightness for all three emission lines
 * at a single (rho, z) grid point for a given CML angle.
 *
 * @param raytracer The raytracer instance
 * @param rho_sky Projected radial position [R_J]
 * @param z_sky Height above equator [R_J]
 * @param cml_deg Central Meridian Longitude [degrees]
 * @param[out] s1p_6731 S+ 6731 Å emission [Rayleighs]
 * @param[out] s2p_680 S++ 680 Å emission [Rayleighs]
 * @param[out] op_833 O+ 833 Å emission [Rayleighs]
 */
void calculate_emission_point(const ipt::JovianUVEmissionRaytracer& raytracer,
                              double rho_sky, double z_sky, double cml_deg,
                              double& s1p_6731, double& s2p_680, double& op_833) {
    
    std::array<double, 3> slit_pos, norm_vec;
    compute_ray_geometry(rho_sky, z_sky, cml_deg, config::R_OBS, slit_pos, norm_vec);

    // Suppress raytracer diagnostic output
    std::streambuf* cout_buf = std::cout.rdbuf();
    std::ostringstream null_stream;
    std::cout.rdbuf(null_stream.rdbuf());
    
    try {
        // S+ 6731 Å (Optical)
        auto result_6731 = raytracer.calculate_spectrum_single(
            slit_pos, norm_vec,
            {config::S1P_6731_MIN, config::S1P_6731_MAX},
            config::S1P_6731_BIN_WIDTH,
            config::S1P_6731_FWHM,
            config::DS
        );
        s1p_6731 = result_6731.total_brightness;
        
        // S++ 680 Å (UV)
        auto result_680 = raytracer.calculate_spectrum_single(
            slit_pos, norm_vec,
            {config::S2P_680_MIN, config::S2P_680_MAX},
            config::S2P_680_BIN_WIDTH,
            config::S2P_680_FWHM,
            config::DS
        );
        s2p_680 = result_680.total_brightness;
        
        // O+ 833 Å (UV)
        auto result_833 = raytracer.calculate_spectrum_single(
            slit_pos, norm_vec,
            {config::OP_833_MIN, config::OP_833_MAX},
            config::OP_833_BIN_WIDTH,
            config::OP_833_FWHM,
            config::DS
        );
        op_833 = result_833.total_brightness;

	// Restore stdout
        std::cout.rdbuf(cout_buf);
        
    } catch (const std::exception&) {
        std::cout.rdbuf(cout_buf);  // Restore stdout on error
        s1p_6731 = 0.0;
        s2p_680 = 0.0;
        op_833 = 0.0;
    }
}

// =============================================================================
// FILE I/O
// =============================================================================

/**
 * @brief Save emission data for a single frame to binary file.
 *
 * File format:
 *   - Header: int32 nrho, int32 nz, float64 cml_deg
 *   - Data: nz × nrho float64 arrays for each emission line
 *
 * @param filename Output filename
 * @param s1p_6731_data S+ 6731 Å emission array [nz × nrho]
 * @param s2p_680_data S++ 680 Å emission array [nz × nrho]
 * @param op_833_data O+ 833 Å emission array [nz × nrho]
 * @param nrho Number of rho grid points
 * @param nz Number of z grid points
 * @param cml_deg CML angle for this frame
 */
void save_frame_binary(const std::string& filename,
                       const std::vector<double>& s1p_6731_data,
                       const std::vector<double>& s2p_680_data,
                       const std::vector<double>& op_833_data,
                       int nrho, int nz, double cml_deg) {
    
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << " for writing." << std::endl;
        return;
    }
    
    // Write header
    file.write(reinterpret_cast<const char*>(&nrho), sizeof(int));
    file.write(reinterpret_cast<const char*>(&nz), sizeof(int));
    file.write(reinterpret_cast<const char*>(&cml_deg), sizeof(double));
    
    // Write emission data
    file.write(reinterpret_cast<const char*>(s1p_6731_data.data()), 
               s1p_6731_data.size() * sizeof(double));
    file.write(reinterpret_cast<const char*>(s2p_680_data.data()),
               s2p_680_data.size() * sizeof(double));
    file.write(reinterpret_cast<const char*>(op_833_data.data()),
               op_833_data.size() * sizeof(double));
    
    file.close();
}

/**
 * @brief Save grid information to text file.
 */
void save_grid_info(const std::string& filename,
                    const std::vector<double>& rho_grid,
                    const std::vector<double>& z_grid) {
    
    std::ofstream file(filename);
    file << "# Grid information for emission movie frames\n";
    file << "nrho = " << rho_grid.size() << "\n";
    file << "nz = " << z_grid.size() << "\n";
    file << "n_frames = " << config::N_FRAMES << "\n";
    file << "frame_rate = " << config::FRAME_RATE << "\n";
    file << "duration = " << config::DURATION << "\n";
    file << "\n# rho_grid [R_J]\n";
    for (double rho : rho_grid) {
        file << std::fixed << std::setprecision(6) << rho << "\n";
    }
    file << "\n# z_grid [R_J]\n";
    for (double z : z_grid) {
        file << std::fixed << std::setprecision(6) << z << "\n";
    }
    file.close();
}

/**
 * @brief Save colorbar limits to text file.
 */
void save_colorbar_limits(const std::string& filename,
                          double s1p_6731_vmax, double s2p_680_vmax, double op_833_vmax) {
    
    std::ofstream file(filename);
    file << "# Colorbar limits for emission movie frames\n";
    file << "# Format: line_name vmin vmax\n";
    file << "s1p_6731 0.0 " << std::fixed << std::setprecision(2) << s1p_6731_vmax << "\n";
    file << "s2p_680 0.0 " << std::fixed << std::setprecision(2) << s2p_680_vmax << "\n";
    file << "op_833 0.0 " << std::fixed << std::setprecision(2) << op_833_vmax << "\n";
    file.close();
}

// =============================================================================
// MAIN FUNCTION
// =============================================================================

int main(int argc, char* argv[]) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    auto program_start = std::chrono::high_resolution_clock::now();
    
    // Only rank 0 prints header information
    if (world_rank == 0) {
        std::cout << "======================================================================\n";
        std::cout << "Io Plasma Torus Hybrid Emission Movie Frame Generator (MPI C++)\n";
        std::cout << "Optical (S+ 6731 Å) + UV (S++ 680 Å, O+ 833 Å)\n";
        std::cout << "Single Maxwellian Distribution\n";
        std::cout << "======================================================================\n\n";
        
        std::cout << "MPI Configuration: " << world_size << " processes\n\n";
    }
    
    // =========================================================================
    // INITIALIZE RAYTRACER (all ranks)
    // =========================================================================
    
    if (world_rank == 0) {
        std::cout << "Initializing raytracer on all ranks...\n" << std::flush;
    }
    
    // Suppress output from non-root ranks during initialization
    std::streambuf* cout_buf = std::cout.rdbuf();
    std::ostringstream null_stream;
    if (world_rank != 0) {
        std::cout.rdbuf(null_stream.rdbuf());
    }
    
    ipt::JovianUVEmissionRaytracer raytracer;
    
    // Restore cout for all ranks
    std::cout.rdbuf(cout_buf);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (world_rank == 0) {
        std::cout << "All ranks initialized successfully.\n\n";
    }
    
    // =========================================================================
    // SETUP GRIDS
    // =========================================================================
    
    std::vector<double> rho_grid = create_rho_grid();
    std::vector<double> z_grid = create_z_grid();
    
    int nrho = static_cast<int>(rho_grid.size());
    int nz = static_cast<int>(z_grid.size());
    int total_points = nrho * nz;
    
    if (world_rank == 0) {
        std::cout << "Grid configuration:\n";
        std::cout << "  Grid: " << nrho << " × " << nz << " = " << total_points << " points per frame\n";
        std::cout << "  ρ range: " << rho_grid.front() << " to " << rho_grid.back() << " R_J\n";
        std::cout << "  z range: " << z_grid.front() << " to " << z_grid.back() << " R_J\n";
        std::cout << "  Total frames: " << config::N_FRAMES << " (covering 360° rotation)\n";
        std::cout << "  CML step: " << std::fixed << std::setprecision(2) 
                  << 360.0 / config::N_FRAMES << "° per frame\n\n";
        
        std::cout << "Spectral parameters:\n";
        std::cout << "  S+ 6731 Å: FWHM=" << config::S1P_6731_FWHM << " Å, bin=" 
                  << config::S1P_6731_BIN_WIDTH << " Å, range=" << config::S1P_6731_MIN 
                  << "-" << config::S1P_6731_MAX << " Å\n";
        std::cout << "  S++ 680 Å: FWHM=" << config::S2P_680_FWHM << " Å, bin=" 
                  << config::S2P_680_BIN_WIDTH << " Å, range=" << config::S2P_680_MIN 
                  << "-" << config::S2P_680_MAX << " Å\n";
        std::cout << "  O+ 833 Å: FWHM=" << config::OP_833_FWHM << " Å, bin=" 
                  << config::OP_833_BIN_WIDTH << " Å, range=" << config::OP_833_MIN 
                  << "-" << config::OP_833_MAX << " Å\n";
        std::cout << "  Integration step: ds=" << config::DS << " R_J\n\n";
        
        // Create output directory
        std::filesystem::create_directories(config::OUTPUT_DIR);
        
        // Save grid info
        save_grid_info(config::OUTPUT_DIR + "/grid_info.txt", rho_grid, z_grid);
        std::cout << "Grid info saved to " << config::OUTPUT_DIR << "/grid_info.txt\n\n";
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    // =========================================================================
    // PRESCAN FOR OPTICAL COLORBAR LIMIT
    // =========================================================================
    
    double s1p_6731_vmax = 0.0;
    
    if (world_rank == 0) {
        std::cout << "======================================================================\n";
        std::cout << "Pre-scanning " << config::PRESCAN_ANGLES.size() 
                  << " CML angles for optical colorbar limit...\n";
        std::cout << "======================================================================\n\n";
    }
    
    for (double prescan_cml : config::PRESCAN_ANGLES) {
        if (world_rank == 0) {
            std::cout << "Scanning CML = " << std::fixed << std::setprecision(0) 
                      << prescan_cml << "°..." << std::flush;
        }
        
        auto scan_start = std::chrono::high_resolution_clock::now();
        
        // Distribute work among ranks
        double local_max = 0.0;
        
        for (int idx = world_rank; idx < total_points; idx += world_size) {
            int i = idx / nz;  // rho index
            int j = idx % nz;  // z index
            
            double s1p_6731, s2p_680, op_833;
            calculate_emission_point(raytracer, rho_grid[i], z_grid[j], prescan_cml,
                                    s1p_6731, s2p_680, op_833);
            
            if (s1p_6731 > local_max) local_max = s1p_6731;
        }
        
        // Reduce to find global max
        double global_max;
        MPI_Reduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        
        if (world_rank == 0) {
            auto scan_end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(scan_end - scan_start);
            
            std::cout << " max = " << std::fixed << std::setprecision(1) << global_max 
                      << " R (" << duration.count() / 1000.0 << " s)\n";
            
            if (global_max > s1p_6731_vmax) {
                s1p_6731_vmax = global_max;
            }
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    // Apply headroom and broadcast
    if (world_rank == 0) {
        double raw_max = s1p_6731_vmax;
        s1p_6731_vmax *= config::HEADROOM_FACTOR;
        
        std::cout << "\nOptical colorbar limit: 0 - " << std::fixed << std::setprecision(1)
                  << s1p_6731_vmax << " R (6% headroom from max " << raw_max << " R)\n\n";
        
        // Save colorbar limits
        save_colorbar_limits(config::OUTPUT_DIR + "/colorbar_limits.txt",
                            s1p_6731_vmax, config::S2P_680_VMAX, config::OP_833_VMAX);
        std::cout << "Colorbar limits saved to " << config::OUTPUT_DIR << "/colorbar_limits.txt\n\n";
    }
    
    MPI_Bcast(&s1p_6731_vmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // =========================================================================
    // GENERATE ALL FRAMES
    // =========================================================================
    
    if (world_rank == 0) {
        std::cout << "======================================================================\n";
        std::cout << "Generating " << config::N_FRAMES << " frames...\n";
        std::cout << "======================================================================\n\n";
    }
    
    auto frames_start = std::chrono::high_resolution_clock::now();
    
    for (int frame_idx = 0; frame_idx < config::N_FRAMES; ++frame_idx) {
        double cml_deg = 360.0 * frame_idx / config::N_FRAMES;
        
        auto frame_start = std::chrono::high_resolution_clock::now();
        
        // Local storage for this rank's computed points
        std::vector<int> local_indices;
        std::vector<double> local_s1p_6731, local_s2p_680, local_op_833;
        
        // Distribute work among ranks
        for (int idx = world_rank; idx < total_points; idx += world_size) {
            int i = idx / nz;  // rho index
            int j = idx % nz;  // z index
            
            double s1p_6731, s2p_680, op_833;
            calculate_emission_point(raytracer, rho_grid[i], z_grid[j], cml_deg,
                                    s1p_6731, s2p_680, op_833);
            
            local_indices.push_back(idx);
            local_s1p_6731.push_back(s1p_6731);
            local_s2p_680.push_back(s2p_680);
            local_op_833.push_back(op_833);
        }
        
        // Gather results on rank 0
        std::vector<double> global_s1p_6731(total_points, 0.0);
        std::vector<double> global_s2p_680(total_points, 0.0);
        std::vector<double> global_op_833(total_points, 0.0);
        
        // Each rank sends its computed values
        int local_count = static_cast<int>(local_indices.size());
        
        if (world_rank == 0) {
            // Rank 0 fills its own values
            for (size_t k = 0; k < local_indices.size(); ++k) {
                int idx = local_indices[k];
                global_s1p_6731[idx] = local_s1p_6731[k];
                global_s2p_680[idx] = local_s2p_680[k];
                global_op_833[idx] = local_op_833[k];
            }
            
            // Receive from other ranks
            for (int src = 1; src < world_size; ++src) {
                int recv_count;
                MPI_Recv(&recv_count, 1, MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                std::vector<int> recv_indices(recv_count);
                std::vector<double> recv_s1p(recv_count), recv_s2p(recv_count), recv_op(recv_count);
                
                MPI_Recv(recv_indices.data(), recv_count, MPI_INT, src, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(recv_s1p.data(), recv_count, MPI_DOUBLE, src, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(recv_s2p.data(), recv_count, MPI_DOUBLE, src, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(recv_op.data(), recv_count, MPI_DOUBLE, src, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                for (int k = 0; k < recv_count; ++k) {
                    int idx = recv_indices[k];
                    global_s1p_6731[idx] = recv_s1p[k];
                    global_s2p_680[idx] = recv_s2p[k];
                    global_op_833[idx] = recv_op[k];
                }
            }
        } else {
            // Send to rank 0
            MPI_Send(&local_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(local_indices.data(), local_count, MPI_INT, 0, 1, MPI_COMM_WORLD);
            MPI_Send(local_s1p_6731.data(), local_count, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
            MPI_Send(local_s2p_680.data(), local_count, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
            MPI_Send(local_op_833.data(), local_count, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
        }
        
        // Rank 0 saves the frame
        if (world_rank == 0) {
            // Reshape to (nz, nrho) layout
            std::vector<double> s1p_6731_frame(total_points);
            std::vector<double> s2p_680_frame(total_points);
            std::vector<double> op_833_frame(total_points);
            
            // Convert from linear idx (rho-major) to (z, rho) layout
            for (int idx = 0; idx < total_points; ++idx) {
                int i = idx / nz;  // rho index
                int j = idx % nz;  // z index
                int out_idx = j * nrho + i;  // (z, rho) layout
                
                s1p_6731_frame[out_idx] = global_s1p_6731[idx];
                s2p_680_frame[out_idx] = global_s2p_680[idx];
                op_833_frame[out_idx] = global_op_833[idx];
            }
            
            // Save binary file
            std::ostringstream filename;
            filename << config::OUTPUT_DIR << "/frame_" << std::setw(4) << std::setfill('0') 
                     << frame_idx << ".bin";
            save_frame_binary(filename.str(), s1p_6731_frame, s2p_680_frame, op_833_frame,
                             nrho, nz, cml_deg);
            
            // Calculate max values for diagnostics
            double max_s1p = *std::max_element(s1p_6731_frame.begin(), s1p_6731_frame.end());
            double max_s2p = *std::max_element(s2p_680_frame.begin(), s2p_680_frame.end());
            double max_op = *std::max_element(op_833_frame.begin(), op_833_frame.end());
            
            auto frame_end = std::chrono::high_resolution_clock::now();
            auto frame_duration = std::chrono::duration_cast<std::chrono::milliseconds>(frame_end - frame_start);
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(frame_end - frames_start);
            
            double avg_time = static_cast<double>(elapsed.count()) / (frame_idx + 1);
            double remaining = avg_time * (config::N_FRAMES - frame_idx - 1);
            
            std::cout << "Frame " << std::setw(2) << frame_idx + 1 << "/" << config::N_FRAMES
                      << " | CML = " << std::fixed << std::setprecision(1) << std::setw(5) << cml_deg << "°"
                      << " | " << std::setprecision(1) << frame_duration.count() / 1000.0 << "s"
                      << " | S+ max=" << std::setw(5) << max_s1p << " R"
                      << " | S++ max=" << std::setw(5) << max_s2p << " R"
                      << " | O+ max=" << std::setw(5) << max_op << " R"
                      << " | ETA: " << std::setprecision(1) << remaining / 60.0 << " min\n" << std::flush;
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    // =========================================================================
    // SUMMARY
    // =========================================================================
    
    auto program_end = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(program_end - program_start);
    
    if (world_rank == 0) {
        std::cout << "\n======================================================================\n";
        std::cout << "Frame Data Generation Complete!\n";
        std::cout << "======================================================================\n\n";
        
        std::cout << "Total time: " << total_duration.count() / 60.0 << " minutes\n";
        std::cout << "Average time per frame: " << total_duration.count() / static_cast<double>(config::N_FRAMES) << " seconds\n\n";
        
        std::cout << "Output directory: " << config::OUTPUT_DIR << "/\n";
        std::cout << "  - " << config::N_FRAMES << " binary frame files (frame_XXXX.bin)\n";
        std::cout << "  - grid_info.txt\n";
        std::cout << "  - colorbar_limits.txt\n\n";
        
        std::cout << "Colorbar limits:\n";
        std::cout << "  S+ 6731 Å:  0 - " << std::fixed << std::setprecision(1) << s1p_6731_vmax << " R\n";
        std::cout << "  S++ 680 Å:  0 - " << config::S2P_680_VMAX << " R\n";
        std::cout << "  O+ 833 Å:   0 - " << config::OP_833_VMAX << " R\n\n";
        
        std::cout << "Next steps:\n";
        std::cout << "  1. Generate PNG frames: python3 generate_frames_from_binary.py\n";
        std::cout << "  2. Create movie: python3 make_movies_ffmpeg.py\n";
        std::cout << "======================================================================\n";
    }
    
    MPI_Finalize();
    return 0;
}
