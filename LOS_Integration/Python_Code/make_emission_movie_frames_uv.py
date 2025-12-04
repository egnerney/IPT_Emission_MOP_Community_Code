#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_emission_movie_frames.py
=============================
Generate UV Emission Movie Frames for Full Jovian Rotation

This script creates image frames showing UV emission from the Jovian plasma torus
as viewed from different observer longitudes (Central Meridian Longitude, CML),
simulating a complete rotation of Jupiter as seen by a distant observer.

MEMORY OPTIMIZATION:
Uses a subclass of JovianUVEmissionRaytracer that loads only single Maxwellian
emission tables (skipping double Maxwellian to save ~4 GB per worker):
  - Plasma model: ~1 GB as float64
  - Single Maxwellian tables: ~250 MB as float64
  - Per-worker memory: ~1.5 GB
  - Total with 23 workers: ~35 GB (fits in 64 GB RAM)

MOVIE OUTPUTS:
Four sets of frames are generated in separate subdirectories:
1. frames_total_uv/  - Total UV emission (550-2100 Å)
2. frames_s2p/       - S++ emission (680 Å ± 3σ)
3. frames_op/        - O+ emission (833 Å ± 3σ)
4. frames_stacked/   - Combined 3-panel stacked view

A separate script (make_movies_ffmpeg.py) combines these into MP4 movies.

MOVIE PARAMETERS:
- Frame rate: 15 fps
- Duration: 5 seconds
- Total frames: 75 (covering 360° rotation)
- CML step: 4.8° per frame

VIEWING GEOMETRY:
For each frame at CML angle θ (degrees), the observer is positioned at longitude θ
relative to Jupiter's +X axis (counterclockwise when viewed from +Z/north):

    CML = 0°:   Observer at X = -R_obs, looking in +X direction
    CML = 90°:  Observer at Y = +R_obs, looking in -Y direction
    CML = 180°: Observer at X = +R_obs, looking in -X direction
    CML = 270°: Observer at Y = -R_obs, looking in +Y direction (ORIGINAL GEOMETRY)

COORDINATE TRANSFORMATION:
For observer at CML = θ viewing a sky position (ρ, z):
  - ρ: projected radial distance from spin axis [R_J]
  - z: height above centrifugal equator [R_J]

Ray direction: [cos(θ), -sin(θ), 0]

Ray starting position:
  start_x = -R_obs·cos(θ) - ρ·sin(θ)
  start_y =  R_obs·sin(θ) - ρ·cos(θ)
  start_z = z

At CML = 270° with ρ = x0, z = z0:
  start = [x0, -R_obs, z0], direction = [0, 1, 0]
  (This exactly matches the original single-frame geometry)

GRID CONFIGURATION:
- ρ position: -10 to +10 R_J in 0.1 R_J steps (201 points)
- z position: -2.5 to +2.5 R_J in 0.1 R_J steps (51 points)
- Total per frame: 201 × 51 = 10,251 calculations

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

# Grid configuration
RHO_MIN = -10.0        # R_J
RHO_MAX = 10.0         # R_J
RHO_STEP = 0.1         # R_J
Z_MIN = -2.5           # R_J
Z_MAX = 2.5            # R_J
Z_STEP = 0.1           # R_J

# Ray tracing parameters
R_OBS = 20.0           # Observer distance from Jupiter [R_J]
FWHM = 6.0             # Instrumental FWHM [Å]
DS = 0.1               # Integration step size [R_J]

# Emission line centers
S3_WAV = 680.0         # S++ (S III) line center [Å]
O2_WAV = 833.0         # O+ (O II) line center [Å]

# Wavelength range
WAV_MIN = 550.0        # Minimum wavelength [Å]
WAV_MAX = 2100.0       # Maximum wavelength [Å]

# Output directories
OUTPUT_BASE = Path("movie_frames")
DIR_TOTAL_UV = OUTPUT_BASE / "frames_total_uv"
DIR_S2P = OUTPUT_BASE / "frames_s2p"
DIR_OP = OUTPUT_BASE / "frames_op"
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
        CML = 270°: Observer at Y = -R_obs, looking in +Y direction (ORIGINAL GEOMETRY)
    
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
    Calculate emissions for a single grid point.
    
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
        (i, j, total_emission, s3_emission, o2_emission)
    """
    global _worker_raytracer
    
    i, j, rho_sky, z_sky, cml_deg = args
    
    if _worker_raytracer is None:
        return i, j, 0.0, 0.0, 0.0
    
    try:
        # Compute ray geometry for this CML
        slit_pos, norm_vec = compute_ray_geometry(rho_sky, z_sky, cml_deg, R_OBS)
        
        # Calculate spectrum using raytracer
        # Note: prints are suppressed by redirecting to /dev/null or NUL
        old_stdout = sys.stdout
        sys.stdout = io.StringIO()
        try:
            wave_bins, spectrum, line_list = _worker_raytracer.calculate_spectrum_single(
                slit_pos, norm_vec,
                wavelength_range=(WAV_MIN, WAV_MAX),
                bin_width=1.0,
                fwhm=FWHM,
                ds=DS
            )
        finally:
            sys.stdout = old_stdout
        
        # Total UV emission (integrated over wavelength)
        total_emission = simpson(spectrum, x=wave_bins)
        
        # Extract S++ and O+ emissions from line list
        sigma = FWHM / 2.35482
        s3_emission = 0.0
        o2_emission = 0.0
        
        for wav, brightness, species in line_list:
            if abs(wav - S3_WAV) <= 3.0 * sigma:
                s3_emission += brightness
            if abs(wav - O2_WAV) <= 3.0 * sigma:
                o2_emission += brightness
        
        return i, j, total_emission, s3_emission, o2_emission
        
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
        (total_emission, s3_emission, o2_emission) arrays, each shape (nz, nrho)
    """
    nrho = len(rho_grid)
    nz = len(z_grid)
    
    total_emission = np.zeros((nz, nrho), dtype=np.float64)
    s3_emission = np.zeros((nz, nrho), dtype=np.float64)
    o2_emission = np.zeros((nz, nrho), dtype=np.float64)
    
    # Prepare work items
    all_args = []
    for i, rho in enumerate(rho_grid):
        for j, z in enumerate(z_grid):
            all_args.append((i, j, rho, z, cml_deg))
    
    # Process with pool using unordered for efficiency
    for result in pool.imap_unordered(calculate_point, all_args, chunksize=50):
        i, j, total, s3, o2 = result
        total_emission[j, i] = total
        s3_emission[j, i] = s3
        o2_emission[j, i] = o2
    
    return total_emission, s3_emission, o2_emission


def prescan_colorbar_limits(rho_grid, z_grid, pool, sample_cmls=[22.5, 112.5, 202.5, 292.5]):#, sample_cmls=[67.5, 157.5, 247.5, 337.5]):#, sample_cmls=[45, 135, 225, 315]):# sample_cmls=[0, 90, 180, 270]):
    """
    Pre-scan representative CML angles to determine global colorbar limits.
    
    Parameters
    ----------
    rho_grid : ndarray
        Array of rho positions [R_J]
    z_grid : ndarray
        Array of z positions [R_J]
    pool : multiprocessing.Pool
        Worker pool for parallel computation
    sample_cmls : list
        CML angles to sample [degrees]
    
    Returns
    -------
    dict
        Dictionary with 'total', 's3', 'o2' keys, each containing (vmin, vmax)
    """
    print(f"\nPre-scanning {len(sample_cmls)} CML angles for colorbar limits...", flush=True)
    
    all_total = []
    all_s3 = []
    all_o2 = []
    
    for cml in sample_cmls:
        print(f"  Scanning CML = {cml:.0f}°...", flush=True)
        t0 = time.time()
        total, s3, o2 = compute_single_frame(cml, rho_grid, z_grid, pool)
        dt = time.time() - t0
        print(f"    Completed in {dt:.1f}s", flush=True)
        all_total.append(total)
        all_s3.append(s3)
        all_o2.append(o2)
    

    # 6% headroom from prescan to prevent saturation but not waste too much dynamic range
    limits = {
        'total': (0.0, 1.06 * max(arr.max() for arr in all_total)),
        's3': (0.0, 1.06 * max(arr.max() for arr in all_s3)),
        'o2': (0.0, 1.06 * max(arr.max() for arr in all_o2))
    }
    
    
    print(f"\nColorbar limits determined:", flush=True)
    print(f"  Total UV: 0 - {limits['total'][1]:.1f} R", flush=True)
    print(f"  S++:      0 - {limits['s3'][1]:.2f} R", flush=True)
    print(f"  O+:       0 - {limits['o2'][1]:.2f} R", flush=True)
    
    return limits


# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

def save_individual_frame(data, rho_grid, z_grid, cml_deg, output_path, 
                          title, vmin, vmax):
    """
    Save a single emission frame as PNG.
    
    Uses viridis colormap with 100 levels, no contour lines.
    Consistent colorbar range across all frames.
    
    Parameters
    ----------
    data : ndarray, shape (nz, nrho)
        Emission data [Rayleighs]
    rho_grid : ndarray
        Array of rho positions [R_J]
    z_grid : ndarray
        Array of z positions [R_J]
    cml_deg : float
        Central Meridian Longitude [degrees]
    output_path : Path
        Output file path
    title : str
        Plot title
    vmin, vmax : float
        Colorbar limits
    """
    RHO, Z = np.meshgrid(rho_grid, z_grid)
    
    fig, ax = plt.subplots(figsize=(10, 4))
    
    # Use contourf with 100 levels, no contour lines
    #levels = np.linspace(vmin, vmax, 256)
    #cs = ax.contourf(RHO, Z, data, levels=levels, cmap='viridis', vmin=vmin, vmax=vmax)
    cs = ax.pcolormesh(RHO, Z, data, cmap='viridis', vmin=vmin, vmax=vmax, shading='gouraud')
    
    ax.set_xlabel(r'$\rho$ (R$_J$)', fontsize=12)
    ax.set_ylabel('Z (R$_J$)', fontsize=12)
    ax.set_title(f'{title} | CML = {cml_deg:.1f}°', fontsize=14)
    ax.set_aspect('equal')
    
    plt.colorbar(cs, ax=ax, label='Brightness (R)')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)


def save_stacked_frame(total_data, s3_data, o2_data, rho_grid, z_grid, cml_deg,
                       output_path, limits):
    """
    Save a stacked 3-panel frame showing all emission types.
    
    Uses viridis colormap with 100 levels, no contour lines.
    Consistent colorbar ranges across all frames.
    
    Parameters
    ----------
    total_data : ndarray, shape (nz, nrho)
        Total UV emission [Rayleighs]
    s3_data : ndarray, shape (nz, nrho)
        S++ emission [Rayleighs]
    o2_data : ndarray, shape (nz, nrho)
        O+ emission [Rayleighs]
    rho_grid : ndarray
        Array of rho positions [R_J]
    z_grid : ndarray
        Array of z positions [R_J]
    cml_deg : float
        Central Meridian Longitude [degrees]
    output_path : Path
        Output file path
    limits : dict
        Colorbar limits dictionary
    """
    RHO, Z = np.meshgrid(rho_grid, z_grid)
    
    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    
    datasets = [
        (total_data, 'Total UV (550-2100 Å)', limits['total']),
        (s3_data, r'S$^{++}$ (680 Å ± 3$\sigma$)', limits['s3']),
        (o2_data, r'O$^{+}$ (833 Å ± 3$\sigma$)', limits['o2'])
    ]
    
    for ax, (data, label, (vmin, vmax)) in zip(axes, datasets):
        #levels = np.linspace(vmin, vmax, 256)
        #cs = ax.contourf(RHO, Z, data, levels=levels, cmap='viridis', vmin=vmin, vmax=vmax)
        cs = ax.pcolormesh(RHO, Z, data, cmap='viridis', vmin=vmin, vmax=vmax, shading='gouraud')
        ax.set_ylabel('Z (R$_J$)', fontsize=11)
        ax.set_title(label, fontsize=12)
        ax.set_aspect('equal')
        plt.colorbar(cs, ax=ax, label='Brightness (R)')
    
    axes[-1].set_xlabel(r'$\rho$ (R$_J$)', fontsize=12)
    
    fig.suptitle(f'Io Plasma Torus UV Emission | CML = {cml_deg:.1f}°', fontsize=14, y=0.98)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)


# =============================================================================
# MAIN FUNCTION
# =============================================================================

def main():
    """Generate all movie frames for Jovian plasma torus rotation."""
    
    print("=" * 70, flush=True)
    print("Io Plasma Torus UV Emission Movie Frame Generator", flush=True)
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
    # SETUP OUTPUT DIRECTORIES
    # =========================================================================
    
    print(flush=True)
    print("Setting up output directories...", flush=True)
    
    for dir_path in [DIR_TOTAL_UV, DIR_S2P, DIR_OP, DIR_STACKED]:
        if dir_path.exists():
            shutil.rmtree(dir_path)
        dir_path.mkdir(parents=True, exist_ok=True)
        print(f"  Created: {dir_path}", flush=True)
    
    # =========================================================================
    # GRID SETUP
    # =========================================================================
    
    print(flush=True)
    print("Setting up calculation grid...", flush=True)
    
    rho_grid = np.arange(RHO_MIN, RHO_MAX + RHO_STEP/2, RHO_STEP)
    z_grid = np.arange(Z_MIN, Z_MAX + Z_STEP/2, Z_STEP)
    
    nrho = len(rho_grid)
    nz = len(z_grid)
    total_per_frame = nrho * nz
    
    print(f"  Grid: {nrho} × {nz} = {total_per_frame} points per frame", flush=True)
    print(f"  Total frames: {N_FRAMES} (covering 360° rotation)", flush=True)
    print(f"  CML step: {360.0/N_FRAMES:.2f}° per frame", flush=True)
    
    # =========================================================================
    # SYSTEM CONFIGURATION
    # =========================================================================
    
    print(flush=True)
    num_cores = mp.cpu_count()
    num_workers = max(1, num_cores - 1) # 16 # this gives 24 - 1 = 23 on my machine for 24 threads even thouygh only have 16 physical cores 
    
    print(f"System: {num_cores} CPU cores, using {num_workers} workers", flush=True)
    print(f"Estimated memory: {num_workers} × ~1.5 GB = ~{num_workers * 1.5:.0f} GB", flush=True)
    
    # =========================================================================
    # CREATE WORKER POOL
    # =========================================================================
    
    print(flush=True)
    print(f"Creating worker pool with {num_workers} workers...", flush=True)
    print("(Each worker loading plasma model and emission tables)", flush=True)
    
    # Use spawn context for Windows compatibility
    #ctx = mp.get_context('spawn')
    #pool = ctx.Pool(processes=num_workers, initializer=init_worker)
    ctx = mp.get_context('spawn')
    print_flag = ctx.Value('b', False)  # shared boolean flag
    pool = ctx.Pool(processes=num_workers, initializer=init_worker, initargs=(print_flag,))
    
    # Wait for workers to initialize
    print("Waiting for workers to initialize...", flush=True)
    time.sleep(10)
    print("Pool created and workers initialized.", flush=True)
    
    # =========================================================================
    # PRE-SCAN FOR COLORBAR LIMITS
    # =========================================================================
    
    #limits = prescan_colorbar_limits(rho_grid, z_grid, pool)
    # Savetime no prescan
    # manual set limits for single maxwellian uv frames found from 0-337.5 for cml with 22.5 degree steps for 16 prescans with about 6% headroom set to nice even values to not staturate but not waste much dynamic range
    limits = {
    'total': (0.0, 2000.0), #rayleighs
    's3': (0.0, 185.0),     #rayleighs
    'o2': (0.0, 175.0)      #rayleighs
}
    
    
    # =========================================================================
    # GENERATE ALL FRAMES
    # =========================================================================
    
    print(flush=True)
    print("=" * 70, flush=True)
    print(f"Generating {N_FRAMES} frames...", flush=True)
    print("=" * 70, flush=True)
    
    start_time = time.time()
    cml_values = np.linspace(0, 360, N_FRAMES, endpoint=False)
    
    for frame_idx, cml_deg in enumerate(cml_values):
        frame_start = time.time()
        
        print(flush=True)
        print(f"Frame {frame_idx + 1}/{N_FRAMES} | CML = {cml_deg:.1f}°", flush=True)
        
        # Compute emission maps
        total_data, s3_data, o2_data = compute_single_frame(cml_deg, rho_grid, z_grid, pool)
        
        # Save individual frames
        save_individual_frame(
            total_data, rho_grid, z_grid, cml_deg,
            DIR_TOTAL_UV / f"frame_{frame_idx:04d}.png",
            "Total UV Emission (550-2100 Å)",
            limits['total'][0], limits['total'][1]
        )
        
        save_individual_frame(
            s3_data, rho_grid, z_grid, cml_deg,
            DIR_S2P / f"frame_{frame_idx:04d}.png",
            r"S$^{++}$ Emission (680 Å ± 3$\sigma$)",
            limits['s3'][0], limits['s3'][1]
        )
        
        save_individual_frame(
            o2_data, rho_grid, z_grid, cml_deg,
            DIR_OP / f"frame_{frame_idx:04d}.png",
            r"O$^{+}$ Emission (833 Å ± 3$\sigma$)",
            limits['o2'][0], limits['o2'][1]
        )
        
        # Save stacked frame
        save_stacked_frame(
            total_data, s3_data, o2_data, rho_grid, z_grid, cml_deg,
            DIR_STACKED / f"frame_{frame_idx:04d}.png",
            limits
        )
        
        frame_time = time.time() - frame_start
        elapsed = time.time() - start_time
        remaining = (N_FRAMES - frame_idx - 1) * (elapsed / (frame_idx + 1))
        
        print(f"  Frame time: {frame_time:.1f}s | Elapsed: {elapsed/60:.1f}min | Remaining: ~{remaining/60:.1f}min", flush=True)
        print(f"  Total UV: max={total_data.max():.1f} R | S++: max={s3_data.max():.2f} R | O+: max={o2_data.max():.2f} R", flush=True)
    
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
    
    print("Output directories:", flush=True)
    print(f"  {DIR_TOTAL_UV} ({len(list(DIR_TOTAL_UV.glob('*.png')))} frames)", flush=True)
    print(f"  {DIR_S2P} ({len(list(DIR_S2P.glob('*.png')))} frames)", flush=True)
    print(f"  {DIR_OP} ({len(list(DIR_OP.glob('*.png')))} frames)", flush=True)
    print(f"  {DIR_STACKED} ({len(list(DIR_STACKED.glob('*.png')))} frames)", flush=True)
    print(flush=True)
    
    print("Colorbar limits used:", flush=True)
    print(f"  Total UV: {limits['total'][0]:.1f} - {limits['total'][1]:.1f} R", flush=True)
    print(f"  S++:      {limits['s3'][0]:.2f} - {limits['s3'][1]:.2f} R", flush=True)
    print(f"  O+:       {limits['o2'][0]:.2f} - {limits['o2'][1]:.2f} R", flush=True)
    print(flush=True)
    
    print("Next step: Run make_movies_ffmpeg.py to create MP4 movies.", flush=True)


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