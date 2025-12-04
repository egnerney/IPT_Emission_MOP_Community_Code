#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_emission_movie_frames_hybrid.py
====================================
Generate Hybrid Optical+UV Emission Movie Frames for Full Jovian Rotation

This script creates image frames showing three specific emission lines from the
Jovian plasma torus in a single stacked 3-panel view, combining one optical line
with two UV lines using optimized spectral windows:

    Top panel:    S+ 6731 Å (optical)
    Middle panel: S++ 680 Å (UV)
    Bottom panel: O+ 833 Å (UV)

MEMORY OPTIMIZATION:
Uses a subclass of JovianUVEmissionRaytracer that loads only single Maxwellian
emission tables (skipping double Maxwellian to save ~4 GB per worker):
  - Plasma model: ~1 GB as float64
  - Single Maxwellian tables: ~250 MB as float64
  - Per-worker memory: ~1.5 GB
  - Total with 23 workers: ~35 GB (fits in 64 GB RAM)

MOVIE OUTPUT:
A single set of stacked 3-panel frames is generated in:
  frames_stacked/  - Combined 3-panel view (S+ 6731, S++ 680, O+ 833)

A separate script (make_movies_ffmpeg_hybrid.py) combines these into an MP4 movie.

MOVIE PARAMETERS:
- Frame rate: 15 fps
- Duration: 5 seconds
- Total frames: 75 (covering 360° rotation)
- CML step: 4.8° per frame

SPECTRAL PARAMETERS:
Each emission line uses optimized spectral windows:

  S+ 6731 Å (Optical):
    - FWHM: 1.2 Å (narrow optical spectrometer)
    - Bin width: 0.6 Å
    - Window: ±10 Å (to safely capture emission line)
    - Wavelength range: 6721 - 6741 Å

  S++ 680 Å (UV):
    - FWHM: 6.0 Å (broader UV)
    - Bin width: 1.0 Å
    - Window: ±10 Å
    - Wavelength range: 670 - 690 Å

  O+ 833 Å (UV):
    - FWHM: 6.0 Å
    - Bin width: 1.0 Å
    - Window: ±10 Å
    - Wavelength range: 823 - 843 Å

VIEWING GEOMETRY:
For each frame at CML angle θ (degrees), the observer is positioned at longitude θ
relative to Jupiter's +X axis (counterclockwise when viewed from +Z/north):

    CML = 0°:   Observer at X = -R_obs, looking in +X direction
    CML = 90°:  Observer at Y = +R_obs, looking in -Y direction
    CML = 180°: Observer at X = +R_obs, looking in -X direction
    CML = 270°: Observer at Y = -R_obs, looking in +Y direction

COORDINATE TRANSFORMATION:
For observer at CML = θ viewing a sky position (ρ, z):
  - ρ: projected radial distance from spin axis [R_J]
  - z: height above centrifugal equator [R_J]

Ray direction: [cos(θ), -sin(θ), 0]

Ray starting position:
  start_x = -R_obs·cos(θ) - ρ·sin(θ)
  start_y =  R_obs·sin(θ) - ρ·cos(θ)
  start_z = z

GRID CONFIGURATION (NON-UNIFORM):
- ρ position: 221 points from -8 to +8 R_J (non-uniform sampling)
- z position: 131 points from -2 to +2 R_J (non-uniform sampling)
- Fine sampling near torus edges and equatorial plane
- Total per frame: 221 × 131 = 28,951 calculations

COLORBAR LIMITS:
- S+ 6731 Å: Automatic (prescan 4 CML angles with 6% headroom)
- S++ 680 Å: Manual (0 - 185 R)
- O+ 833 Å: Manual (0 - 175 R)

PHYSICS:
Uses the JovianUVEmissionRaytracer class from IPT_emiss_MOP_community_code.py:
1. Ray tracing through 3D plasma model (trilinear interpolation)
2. Emission rate lookup from CHIANTI tables (bilinear in log-log space)
3. Simpson's rule integration along line of sight
4. ERF-based Gaussian convolution for instrument response

AUTHOR: Edward (Eddie) G. Nerney
INSTITUTION: LASP, University of Colorado Boulder
DATE: November 2025
LICENSE: MIT
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import sys
from pathlib import Path
import time
import multiprocessing as mp
from multiprocessing import Value
import os
import shutil
import io

# =============================================================================
# CONFIGURATION
# =============================================================================

# Movie parameters
FRAME_RATE = 15        # frames per second
DURATION = 5           # seconds
N_FRAMES = FRAME_RATE * DURATION  # 75 frames for full rotation

# Grid configuration - Non-uniform optimized for torus structure
# Radial grid: 221 points from -8 to +8 R_J
RHO_SEGMENTS = [
    (-8.0, -6.0, 0.1),      # Outer edge
    (-6.0, -4.5, 0.025),    # Transition region (fine)
    (-4.5, -3.0, 0.1),      # Intermediate
    (-3.0, 3.0, 0.2),       # Central torus (coarse)
    (3.0, 4.5, 0.1),        # Intermediate
    (4.5, 6.0, 0.025),      # Transition region (fine)
    (6.0, 8.0, 0.1),        # Outer edge
]

# Vertical grid: 131 points from -2 to +2 R_J
Z_SEGMENTS = [
    (-2.0, -0.5, 0.1),      # Lower region
    (-0.5, 0.5, 0.01),      # Equatorial plane (very fine)
    (0.5, 2.0, 0.1),        # Upper region
]

# Ray tracing parameters
R_OBS = 20.0           # Observer distance from Jupiter [R_J]
DS = 0.1               # Integration step size [R_J]

# -------------------------------------------------------------------------
# Optical line: S+ 6731 Å
# -------------------------------------------------------------------------
S1P_6731_WAV = 6731.0          # Line center [Å]
S1P_6731_FWHM = 1.2            # Instrumental FWHM [Å] (narrow optical)
S1P_6731_BIN_WIDTH = 0.6       # Bin width [Å]
# Window: ±10 Å to safely capture emission line
S1P_6731_RANGE = (S1P_6731_WAV - 10.0, S1P_6731_WAV + 10.0)

# -------------------------------------------------------------------------
# UV line 1: S++ 680 Å
# -------------------------------------------------------------------------
S2P_680_WAV = 680.0            # Line center [Å]
S2P_680_FWHM = 6.0             # Instrumental FWHM [Å] (broader UV)
S2P_680_BIN_WIDTH = 1.0        # Bin width [Å]
# Window: ±10 Å
S2P_680_RANGE = (S2P_680_WAV - 10.0, S2P_680_WAV + 10.0)

# -------------------------------------------------------------------------
# UV line 2: O+ 833 Å
# -------------------------------------------------------------------------
OP_833_WAV = 833.0             # Line center [Å]
OP_833_FWHM = 6.0              # Instrumental FWHM [Å]
OP_833_BIN_WIDTH = 1.0         # Bin width [Å]
# Window: ±10 Å
OP_833_RANGE = (OP_833_WAV - 10.0, OP_833_WAV + 10.0)

# -------------------------------------------------------------------------
# Manual colorbar limits for UV lines (from UV example code)
# -------------------------------------------------------------------------
UV_LIMITS = {
    's2p_680': (0.0, 185.0),   # Rayleighs
    'op_833': (0.0, 175.0)     # Rayleighs
}

# Prescan CML angles for optical line colorbar limit
PRESCAN_ANGLES = [0.0, 90.0, 180.0, 270.0]

# Output directory
OUTPUT_BASE = Path("movie_frames")
DIR_STACKED = OUTPUT_BASE / "frames_stacked"

# =============================================================================
# PATH SETUP
# =============================================================================

SCRIPT_DIR = Path(__file__).resolve().parent
RAYTRACER_DIR = SCRIPT_DIR.parent / "IPT_emission_model"

# Add raytracer directory to Python path
if str(RAYTRACER_DIR) not in sys.path:
    sys.path.insert(0, str(RAYTRACER_DIR))


# =============================================================================
# SINGLE MAXWELLIAN RAYTRACER SUBCLASS
# =============================================================================

def create_single_maxwellian_raytracer_class():
    """
    Create a subclass of JovianUVEmissionRaytracer that skips double Maxwellian
    table loading to save memory (~4 GB per instance).
    
    Returns
    -------
    class
        SingleMaxwellianRaytracer class
    """
    from IPT_emiss_MOP_community_code import JovianUVEmissionRaytracer
    
    class SingleMaxwellianRaytracer(JovianUVEmissionRaytracer):
        """
        Memory-efficient raytracer that loads only single Maxwellian tables.
        
        Saves ~4 GB per instance by skipping double Maxwellian emission tables.
        Uses calculate_spectrum_single() for all calculations.
        """
        
        def _get_default_emission_path_double(self):
            """Override to return None, skipping double Maxwellian loading."""
            return None
    
    return SingleMaxwellianRaytracer


# =============================================================================
# GRID FUNCTIONS
# =============================================================================

def create_nonuniform_grid(segments):
    """
    Create a non-uniform grid from segment specifications.
    
    Parameters
    ----------
    segments : list of tuples
        Each tuple is (start, stop, step) defining a segment
    
    Returns
    -------
    ndarray
        Concatenated grid points
    """
    parts = []
    for i, (start, stop, step) in enumerate(segments):
        if i < len(segments) - 1:
            # Don't include endpoint for intermediate segments
            parts.append(np.arange(start, stop, step))
        else:
            # Include endpoint for last segment
            parts.append(np.arange(start, stop + step/2, step))
    
    return np.concatenate(parts)


# =============================================================================
# GEOMETRY FUNCTIONS
# =============================================================================

def compute_ray_geometry(rho_sky, z_sky, cml_deg, R_obs=20.0):
    """
    Compute ray starting position and direction for given sky coordinates and CML.
    
    The observer rotates around Jupiter's spin axis (Z). The Central Meridian
    Longitude (CML) is the sub-observer longitude - the longitude on Jupiter
    directly beneath the observer.
    
    CML Convention (longitude measured counterclockwise from +X when viewed from +Z):
        CML = 0°:   Observer at X = -R_obs, looking in +X direction
        CML = 90°:  Observer at Y = +R_obs, looking in -Y direction
        CML = 180°: Observer at X = +R_obs, looking in -X direction
        CML = 270°: Observer at Y = -R_obs, looking in +Y direction
    
    Parameters
    ----------
    rho_sky : float
        Projected radial position in observer's sky plane [R_J]
        (horizontal axis in output image, perpendicular to LOS and Z)
    z_sky : float
        Height above/below centrifugal equator [R_J]
        (vertical axis in output image)
    cml_deg : float
        Central Meridian Longitude [degrees]
    R_obs : float
        Observer distance from Jupiter center [R_J]
    
    Returns
    -------
    slit_pos : ndarray, shape (3,)
        Ray starting position [x, y, z] in R_J
    norm_vec : ndarray, shape (3,)
        Ray direction vector (normalized)
    
    Notes
    -----
    At CML = 270° with rho_sky = x0, z_sky = z0:
        slit_pos = [x0, -R_obs, z0]
        norm_vec = [0, 1, 0]
    This exactly matches the original single-frame geometry.
    """
    theta = np.radians(cml_deg)
    
    # Viewing direction (toward Jupiter center)
    norm_vec = np.array([np.cos(theta), -np.sin(theta), 0.0], dtype=np.float64)
    
    # Ray starting position
    start_x = -R_obs * np.cos(theta) - rho_sky * np.sin(theta)
    start_y =  R_obs * np.sin(theta) - rho_sky * np.cos(theta)
    start_z = z_sky
    
    slit_pos = np.array([start_x, start_y, start_z], dtype=np.float64)
    
    return slit_pos, norm_vec


# =============================================================================
# WORKER FUNCTIONS
# =============================================================================

# Global worker state
_worker_raytracer = None


def init_worker(print_flag):
    """
    Initialize worker process with its own raytracer instance.
    
    Each worker loads the plasma model and single Maxwellian emission tables.
    Memory usage per worker: ~1.5 GB
    
    Parameters
    ----------
    print_flag : multiprocessing.Value
        Shared flag to control diagnostic printing (first worker prints)
    """
    global _worker_raytracer
    
    pid = os.getpid()
    print(f"[Worker {pid}] Starting initialization...", flush=True)
    
    # First worker to grab flag prints, others suppress
    with print_flag.get_lock():
        should_print = not print_flag.value
        print_flag.value = True
    
    try:
        t0 = time.time()
        
        if not should_print:
            old_stdout = sys.stdout
            sys.stdout = io.StringIO()
        
        SingleMaxwellianRaytracer = create_single_maxwellian_raytracer_class()
        _worker_raytracer = SingleMaxwellianRaytracer()
        
        if not should_print:
            sys.stdout = old_stdout
        
        dt = time.time() - t0
        
        print(f"[Worker {pid}] Ready in {dt:.1f}s", flush=True)
        
    except Exception as e:
        if not should_print:
            sys.stdout = old_stdout
        print(f"[Worker {pid}] ERROR: {e}", flush=True)
        import traceback
        traceback.print_exc()
        raise


def calculate_point(args):
    """
    Calculate three emission lines for a single grid point.
    
    Computes S+ 6731 Å (optical), S++ 680 Å (UV), and O+ 833 Å (UV)
    using optimized spectral windows with different FWHM and bin widths.
    
    Parameters
    ----------
    args : tuple
        (i, j, rho_sky, z_sky, cml_deg) where:
        - i: rho index
        - j: z index
        - rho_sky: projected radial position [R_J]
        - z_sky: height above equator [R_J]
        - cml_deg: Central Meridian Longitude [degrees]
    
    Returns
    -------
    tuple
        (i, j, s1p_6731_emission, s2p_680_emission, op_833_emission)
        Emissions are in Rayleighs
    """
    global _worker_raytracer
    
    i, j, rho_sky, z_sky, cml_deg = args
    
    if _worker_raytracer is None:
        return i, j, 0.0, 0.0, 0.0
    
    try:
        # Compute ray geometry for this CML
        slit_pos, norm_vec = compute_ray_geometry(rho_sky, z_sky, cml_deg, R_OBS)
        
        # Suppress raytracer stdout messages
        old_stdout = sys.stdout
        sys.stdout = io.StringIO()
        
        try:
            # -----------------------------------------------------------------
            # S+ 6731 Å (Optical - narrow FWHM, fine bins)
            # -----------------------------------------------------------------
            wave_bins_6731, spectrum_6731, _ = _worker_raytracer.calculate_spectrum_single(
                slit_pos, norm_vec,
                wavelength_range=S1P_6731_RANGE,
                bin_width=S1P_6731_BIN_WIDTH,
                fwhm=S1P_6731_FWHM,
                ds=DS
            )
            s1p_6731_emission = simpson(spectrum_6731, x=wave_bins_6731)
            
            # -----------------------------------------------------------------
            # S++ 680 Å (UV - broader FWHM)
            # -----------------------------------------------------------------
            wave_bins_680, spectrum_680, _ = _worker_raytracer.calculate_spectrum_single(
                slit_pos, norm_vec,
                wavelength_range=S2P_680_RANGE,
                bin_width=S2P_680_BIN_WIDTH,
                fwhm=S2P_680_FWHM,
                ds=DS
            )
            s2p_680_emission = simpson(spectrum_680, x=wave_bins_680)
            
            # -----------------------------------------------------------------
            # O+ 833 Å (UV - broader FWHM)
            # -----------------------------------------------------------------
            wave_bins_833, spectrum_833, _ = _worker_raytracer.calculate_spectrum_single(
                slit_pos, norm_vec,
                wavelength_range=OP_833_RANGE,
                bin_width=OP_833_BIN_WIDTH,
                fwhm=OP_833_FWHM,
                ds=DS
            )
            op_833_emission = simpson(spectrum_833, x=wave_bins_833)
            
        finally:
            sys.stdout = old_stdout
        
        return i, j, s1p_6731_emission, s2p_680_emission, op_833_emission
        
    except Exception:
        return i, j, 0.0, 0.0, 0.0


# =============================================================================
# FRAME COMPUTATION
# =============================================================================

def compute_single_frame(cml_deg, rho_grid, z_grid, pool):
    """
    Compute emission maps for a single CML angle.
    
    Parameters
    ----------
    cml_deg : float
        Central Meridian Longitude [degrees]
    rho_grid : ndarray
        Array of rho positions [R_J]
    z_grid : ndarray
        Array of z positions [R_J]
    pool : multiprocessing.Pool
        Worker pool for parallel computation
    
    Returns
    -------
    tuple
        (s1p_6731_emission, s2p_680_emission, op_833_emission) arrays,
        each with shape (nz, nrho) in Rayleighs
    """
    nrho = len(rho_grid)
    nz = len(z_grid)
    
    s1p_6731_emission = np.zeros((nz, nrho), dtype=np.float64)
    s2p_680_emission = np.zeros((nz, nrho), dtype=np.float64)
    op_833_emission = np.zeros((nz, nrho), dtype=np.float64)
    
    # Prepare work items
    all_args = []
    for i, rho in enumerate(rho_grid):
        for j, z in enumerate(z_grid):
            all_args.append((i, j, rho, z, cml_deg))
    
    # Process with pool using unordered for efficiency
    for result in pool.imap_unordered(calculate_point, all_args, chunksize=50):
        i, j, s1p_6731, s2p_680, op_833 = result
        s1p_6731_emission[j, i] = s1p_6731
        s2p_680_emission[j, i] = s2p_680
        op_833_emission[j, i] = op_833
    
    return s1p_6731_emission, s2p_680_emission, op_833_emission


def prescan_optical_colorbar_limit(rho_grid, z_grid, pool, sample_cmls=None):
    """
    Pre-scan representative CML angles to determine optical line colorbar limit.
    
    Only the S+ 6731 Å line needs prescanning; UV lines use manual limits.
    
    Parameters
    ----------
    rho_grid : ndarray
        Array of rho positions [R_J]
    z_grid : ndarray
        Array of z positions [R_J]
    pool : multiprocessing.Pool
        Worker pool for parallel computation
    sample_cmls : list, optional
        CML angles to sample [degrees]. Default: [0, 90, 180, 270]
    
    Returns
    -------
    float
        Maximum emission value for S+ 6731 Å with 6% headroom
    """
    if sample_cmls is None:
        sample_cmls = PRESCAN_ANGLES
    
    print(f"\nPre-scanning {len(sample_cmls)} CML angles for optical colorbar limit...", flush=True)
    
    all_s1p_6731_max = []
    
    for cml in sample_cmls:
        print(f"  Scanning CML = {cml:.0f}°...", flush=True)
        t0 = time.time()
        s1p_6731, s2p_680, op_833 = compute_single_frame(cml, rho_grid, z_grid, pool)
        dt = time.time() - t0
        print(f"    Completed in {dt:.1f}s | S+ 6731 max = {s1p_6731.max():.1f} R", flush=True)
        all_s1p_6731_max.append(s1p_6731.max())
    
    # Apply 6% headroom to prevent saturation
    optical_vmax = 1.06 * max(all_s1p_6731_max)
    
    print(f"\nOptical colorbar limit determined:", flush=True)
    print(f"  S+ 6731 Å: 0 - {optical_vmax:.1f} R (6% headroom from max {max(all_s1p_6731_max):.1f} R)", flush=True)
    
    return optical_vmax


# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

def save_stacked_frame(s1p_6731_data, s2p_680_data, op_833_data, 
                       rho_grid, z_grid, cml_deg, output_path, limits):
    """
    Save a stacked 3-panel frame showing all three emission lines.
    
    Panel layout:
      - Top: S+ 6731 Å (optical)
      - Middle: S++ 680 Å (UV)
      - Bottom: O+ 833 Å (UV)
    
    Uses viridis colormap with smooth interpolation (gouraud shading).
    Consistent colorbar ranges across all frames.
    
    Parameters
    ----------
    s1p_6731_data : ndarray, shape (nz, nrho)
        S+ 6731 Å emission [Rayleighs]
    s2p_680_data : ndarray, shape (nz, nrho)
        S++ 680 Å emission [Rayleighs]
    op_833_data : ndarray, shape (nz, nrho)
        O+ 833 Å emission [Rayleighs]
    rho_grid : ndarray
        Array of rho positions [R_J]
    z_grid : ndarray
        Array of z positions [R_J]
    cml_deg : float
        Central Meridian Longitude [degrees]
    output_path : Path
        Output file path
    limits : dict
        Colorbar limits dictionary with keys 's1p_6731', 's2p_680', 'op_833'
    """
    RHO, Z = np.meshgrid(rho_grid, z_grid)
    
    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    
    # Define datasets with labels and limits
    datasets = [
        (s1p_6731_data, r'S$^{+}$ 6731 Å (Optical)', limits['s1p_6731']),
        (s2p_680_data, r'S$^{++}$ 680 Å (UV)', limits['s2p_680']),
        (op_833_data, r'O$^{+}$ 833 Å (UV)', limits['op_833'])
    ]
    
    for ax, (data, label, (vmin, vmax)) in zip(axes, datasets):
        cs = ax.pcolormesh(RHO, Z, data, cmap='viridis', vmin=vmin, vmax=vmax, shading='gouraud')
        ax.set_ylabel('Z (R$_J$)', fontsize=11)
        ax.set_title(label, fontsize=12)
        ax.set_aspect('equal')
        plt.colorbar(cs, ax=ax, label='Brightness (R)')
    
    axes[-1].set_xlabel(r'$\rho$ (R$_J$)', fontsize=12)
    
    fig.suptitle(f'Io Plasma Torus Emission | CML = {cml_deg:.1f}°', fontsize=14, y=0.98)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)


# =============================================================================
# MAIN FUNCTION
# =============================================================================

def main():
    """Generate all movie frames for Jovian plasma torus rotation."""
    
    print("=" * 70, flush=True)
    print("Io Plasma Torus Hybrid Emission Movie Frame Generator", flush=True)
    print("Optical (S+ 6731 Å) + UV (S++ 680 Å, O+ 833 Å)", flush=True)
    print("Using JovianUVEmissionRaytracer (Single Maxwellian)", flush=True)
    print("=" * 70, flush=True)
    print(flush=True)
    
    # =========================================================================
    # VERIFY RAYTRACER AVAILABILITY
    # =========================================================================
    
    print("Checking raytracer availability...", flush=True)
    
    try:
        old_stdout = sys.stdout
        sys.stdout = io.StringIO()
        try:
            SingleMaxwellianRaytracer = create_single_maxwellian_raytracer_class()
            test_raytracer = SingleMaxwellianRaytracer()
        finally:
            sys.stdout = old_stdout
        print("  Raytracer loaded successfully", flush=True)
        print(f"  Single Maxwellian loaded: {test_raytracer.single_maxwellian_loaded}", flush=True)
        print(f"  Double Maxwellian loaded: {test_raytracer.double_maxwellian_loaded}", flush=True)
        del test_raytracer
    except Exception as e:
        print(f"ERROR: Could not load raytracer: {e}", flush=True)
        print(f"  Expected location: {RAYTRACER_DIR}", flush=True)
        return
    
    # =========================================================================
    # SETUP OUTPUT DIRECTORY
    # =========================================================================
    
    print(flush=True)
    print("Setting up output directory...", flush=True)
    
    if DIR_STACKED.exists():
        shutil.rmtree(DIR_STACKED)
    DIR_STACKED.mkdir(parents=True, exist_ok=True)
    print(f"  Created: {DIR_STACKED}", flush=True)
    
    # =========================================================================
    # GRID SETUP
    # =========================================================================
    
    print(flush=True)
    print("Setting up calculation grid...", flush=True)
    
    rho_grid = create_nonuniform_grid(RHO_SEGMENTS)
    z_grid = create_nonuniform_grid(Z_SEGMENTS)
    
    nrho = len(rho_grid)
    nz = len(z_grid)
    total_per_frame = nrho * nz
    
    print(f"  Grid: {nrho} × {nz} = {total_per_frame} points per frame", flush=True)
    print(f"  ρ range: {rho_grid.min():.1f} to {rho_grid.max():.1f} R_J", flush=True)
    print(f"  z range: {z_grid.min():.1f} to {z_grid.max():.1f} R_J", flush=True)
    print(f"  Total frames: {N_FRAMES} (covering 360° rotation)", flush=True)
    print(f"  CML step: {360.0/N_FRAMES:.2f}° per frame", flush=True)
    
    # =========================================================================
    # SPECTRAL PARAMETERS
    # =========================================================================
    
    print(flush=True)
    print("Spectral parameters:", flush=True)
    print(f"  S+ 6731 Å: FWHM={S1P_6731_FWHM} Å, bin={S1P_6731_BIN_WIDTH} Å, range={S1P_6731_RANGE}", flush=True)
    print(f"  S++ 680 Å: FWHM={S2P_680_FWHM} Å, bin={S2P_680_BIN_WIDTH} Å, range={S2P_680_RANGE}", flush=True)
    print(f"  O+ 833 Å:  FWHM={OP_833_FWHM} Å, bin={OP_833_BIN_WIDTH} Å, range={OP_833_RANGE}", flush=True)
    print(f"  Integration step: ds={DS} R_J", flush=True)
    
    # =========================================================================
    # SYSTEM CONFIGURATION
    # =========================================================================
    
    print(flush=True)
    num_cores = mp.cpu_count()
    num_workers = max(1, num_cores - 1)
    
    print(f"System: {num_cores} CPU cores, using {num_workers} workers", flush=True)
    print(f"Estimated memory: {num_workers} × ~1.5 GB = ~{num_workers * 1.5:.0f} GB", flush=True)
    
    # =========================================================================
    # CREATE WORKER POOL
    # =========================================================================
    
    print(flush=True)
    print(f"Creating worker pool with {num_workers} workers...", flush=True)
    print("(Each worker loading plasma model and emission tables)", flush=True)
    
    ctx = mp.get_context('spawn')
    print_flag = ctx.Value('b', False)
    pool = ctx.Pool(processes=num_workers, initializer=init_worker, initargs=(print_flag,))
    
    # Wait for workers to initialize
    print("Waiting for workers to initialize...", flush=True)
    time.sleep(10)
    print("Pool created and workers initialized.", flush=True)
    
    # =========================================================================
    # DETERMINE COLORBAR LIMITS
    # =========================================================================
    
    # Pre-scan for optical line limit
    optical_vmax = prescan_optical_colorbar_limit(rho_grid, z_grid, pool)
    
    # Build complete limits dictionary
    limits = {
        's1p_6731': (0.0, optical_vmax),
        's2p_680': UV_LIMITS['s2p_680'],
        'op_833': UV_LIMITS['op_833']
    }
    
    print(flush=True)
    print("Final colorbar limits:", flush=True)
    print(f"  S+ 6731 Å:  {limits['s1p_6731'][0]:.1f} - {limits['s1p_6731'][1]:.1f} R (prescanned)", flush=True)
    print(f"  S++ 680 Å:  {limits['s2p_680'][0]:.1f} - {limits['s2p_680'][1]:.1f} R (manual)", flush=True)
    print(f"  O+ 833 Å:   {limits['op_833'][0]:.1f} - {limits['op_833'][1]:.1f} R (manual)", flush=True)
    
    # =========================================================================
    # GENERATE ALL FRAMES
    # =========================================================================
    
    print(flush=True)
    print("=" * 70, flush=True)
    print(f"Generating {N_FRAMES} stacked frames...", flush=True)
    print("=" * 70, flush=True)
    
    start_time = time.time()
    cml_values = np.linspace(0, 360, N_FRAMES, endpoint=False)
    
    for frame_idx, cml_deg in enumerate(cml_values):
        frame_start = time.time()
        
        print(flush=True)
        print(f"Frame {frame_idx + 1}/{N_FRAMES} | CML = {cml_deg:.1f}°", flush=True)
        
        # Compute emission maps
        s1p_6731_data, s2p_680_data, op_833_data = compute_single_frame(
            cml_deg, rho_grid, z_grid, pool
        )
        
        # Save stacked frame
        save_stacked_frame(
            s1p_6731_data, s2p_680_data, op_833_data, 
            rho_grid, z_grid, cml_deg,
            DIR_STACKED / f"frame_{frame_idx:04d}.png",
            limits
        )
        
        frame_time = time.time() - frame_start
        elapsed = time.time() - start_time
        remaining = (N_FRAMES - frame_idx - 1) * (elapsed / (frame_idx + 1))
        
        print(f"  Frame time: {frame_time:.1f}s | Elapsed: {elapsed/60:.1f}min | Remaining: ~{remaining/60:.1f}min", flush=True)
        print(f"  S+ 6731: max={s1p_6731_data.max():.1f} R | S++ 680: max={s2p_680_data.max():.1f} R | O+ 833: max={op_833_data.max():.1f} R", flush=True)
    
    # =========================================================================
    # CLEANUP
    # =========================================================================
    
    pool.close()
    pool.join()
    
    total_time = time.time() - start_time
    
    # =========================================================================
    # SUMMARY
    # =========================================================================
    
    print(flush=True)
    print("=" * 70, flush=True)
    print("Frame Generation Complete!", flush=True)
    print("=" * 70, flush=True)
    print(flush=True)
    
    print(f"Total computation time: {total_time/60:.1f} minutes ({total_time/3600:.2f} hours)", flush=True)
    print(f"Average time per frame: {total_time/N_FRAMES:.1f} seconds", flush=True)
    print(flush=True)
    
    print("Output directory:", flush=True)
    print(f"  {DIR_STACKED} ({len(list(DIR_STACKED.glob('*.png')))} frames)", flush=True)
    print(flush=True)
    
    print("Colorbar limits used:", flush=True)
    print(f"  S+ 6731 Å:  {limits['s1p_6731'][0]:.1f} - {limits['s1p_6731'][1]:.1f} R", flush=True)
    print(f"  S++ 680 Å:  {limits['s2p_680'][0]:.1f} - {limits['s2p_680'][1]:.1f} R", flush=True)
    print(f"  O+ 833 Å:   {limits['op_833'][0]:.1f} - {limits['op_833'][1]:.1f} R", flush=True)
    print(flush=True)
    
    print("Next step: Run make_movies_ffmpeg_hybrid.py to create MP4 movie.", flush=True)


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    mp.freeze_support()
    
    try:
        main()
    except KeyboardInterrupt:
        print("\nInterrupted by user.", flush=True)
    except Exception as e:
        print(f"\nUnexpected error: {e}", flush=True)
        import traceback
        traceback.print_exc()
