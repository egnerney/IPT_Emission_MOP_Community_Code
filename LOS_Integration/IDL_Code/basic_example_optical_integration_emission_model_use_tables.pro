;+
; basic_example_optical_integration_emission_model_use_tables.pro
;
; Example Script for Optical Emission Line-of-Sight Integration Through Io Plasma Torus
; ====================================================================================
;
; This script demonstrates optical emission calculations through the Io Plasma Torus
; using proper line-of-sight integration with ray tracing and 3D plasma model
; interpolation. Designed for ground-based telescope observations with realistic
; test geometries for both single and double Maxwellian electron distributions.
;
; Physical Model:
; ---------------
; The calculation follows these steps at each point along the line of sight:
;
; 1. Ray Tracing: Create ray from observer through plasma torus with configurable
;    step size (default 0.01 R_J for accurate integration)
;
; 2. 3D Plasma Interpolation: At each ray point, trilinear interpolation retrieves:
;    - Electron temperature Te [eV]
;    - Electron density ne [cm^-3]
;    - Ion number densities n_ion [cm^-3] for each species
;
; 3. Emissivity Interpolation: From pre-calculated CHIANTI tables:
;    - Single Maxwellian: Bilinear interpolation in log10(Te)-log10(ne) space
;    - Double Maxwellian: Quadlinear interpolation in (Tec, Teh, ne, feh) space
;    - Returns volume emissivity [photons s^-1 cm^-3] for all emission lines
;
; 4. Per-Ion Emission Rate: Convert volume emissivity to per-ion rate:
;    per_ion_rate = volume_emissivity / ne  [photons s^-1 ion^-1]
;
; 5. Brightness Integrand: Multiply by ion number density:
;    integrand = per_ion_rate * n_ion  [photons s^-1 cm^-3]
;
; 6. Line-of-Sight Integration: Simpson's rule integration:
;    brightness = integral(integrand) * R_J_cm * 1e-6  [Rayleighs]
;    where R_J_cm = 7.1492e9 cm and 1e-6 converts to Rayleighs
;
; 7. Instrumental Response: Convolve discrete line brightnesses with Gaussian
;    instrumental response using analytic ERF formulation:
;    spectrum = sum over lines of (brightness_i * gaussian_response_i)
;    Final output in R/Angstrom (brightness per wavelength)
;
; Test Cases:
; -----------
; 1. Equatorial line of sight (single Maxwellian)
; 2. Off-equator line of sight (single Maxwellian)
; 3. Equatorial line of sight (double Maxwellian)
; 4. Off-equator line of sight (double Maxwellian)
;
; All lines of sight start at x=6 R_J and look through the plasma torus
; in the +y direction, sampling the peak emission region near Io's orbit.
;
; REQUIRED DATA FILES:
; - Plasma Model: ../../3D_Torus_Model/jovian_plasma_interpolated_381x381x231.h5
; - Single Maxwellian: ../../Emiss_tables/CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5
; - Double Maxwellian: ../../Emiss_tables/CHIANTI_11.0.2_emiss_tables_double_maxwellian_16x8x12x10.h5
;
; OUTPUT:
; - Optical emission spectra for single and double Maxwellian cases
; - Comparison plots showing hot electron enhancement
; - Analysis of key optical diagnostic lines ([S II], [O II], [O III], [S III])
; - Publication-quality plots saved as PNG files
;
; WAVELENGTH RANGE:
; - Optical: 3000-10000 Angstroms (ground-based telescopes)
;
; Key Optical Diagnostic Lines:
; - [O II] 3726/3729 Angstrom doublet (density diagnostic)
; - [S II] 4069/4076 Angstrom (auroral lines)
; - [O III] 4363 Angstrom (temperature diagnostic)
; - [O III] 4959/5007 Angstrom (nebular lines)
; - [S III] 6312 Angstrom (temperature diagnostic)
; - [S II] 6716/6731 Angstrom (density diagnostic doublet)
; - [S III] 9069 Angstrom (near-IR forbidden line)
;
; AUTHOR: Edward (Eddie) G. Nerney
; INSTITUTION: Laboratory for Atmospheric and Space Physics, University of Colorado Boulder
; VERSION: 1.0
; DATE: November 2025
;
; ACKNOWLEDGMENTS:
; This code uses emission tables calculated from the CHIANTI atomic database
; (Dere et al. 1997, 2019). Please cite CHIANTI when using this code.
;-

; Load the IPT emission raytracer library
@IPT_emiss_MOP_community_code

PRO basic_example_optical_integration_emission_model_use_tables
  ;+
  ; NAME:
  ;   basic_example_optical_integration_emission_model_use_tables
  ;
  ; PURPOSE:
  ;   Main procedure demonstrating optical emission calculations through the Io Plasma
  ;   Torus using proper line-of-sight integration with ray tracing and 3D plasma
  ;   model interpolation. Calculates per-ion photon emission rates and integrates
  ;   over the line of sight to obtain brightness in Rayleighs for ground-based
  ;   telescope observations.
  ;
  ; CALLING SEQUENCE:
  ;   basic_example_optical_integration_emission_model_use_tables
  ;
  ; INPUTS:
  ;   None (data files are loaded from standard directories)
  ;
  ; OUTPUTS:
  ;   Creates plots showing simulated optical spectra in Rayleighs/Angstrom
  ;   Saves PNG files of the spectra
  ;   Prints summary of key optical line diagnostics
  ;
  ; REQUIRED FILES:
  ;   Plasma model HDF5 file
  ;   Single Maxwellian emission table HDF5 file
  ;   Double Maxwellian emission table HDF5 file
  ;
  ; MODIFICATION HISTORY:
  ;   Written by E. G. Nerney, November 2025
  ;   Version 2.0: Proper LOS integration with per-ion emission rates
  ;-

  ; Record start time for performance measurement
  start_time = SYSTIME(/SECONDS)

  PRINT, '======================================================================'
  PRINT, 'JOVIAN OPTICAL EMISSION RAYTRACER - LINE-OF-SIGHT INTEGRATION'
  PRINT, 'Io Plasma Torus Optical Emission Model for Ground-Based Telescopes'
  PRINT, 'Version 2.0 - Proper Per-Ion Emission Rate Integration'
  PRINT, '======================================================================'
  PRINT, ''
  PRINT, 'This example demonstrates optical emission calculations through the'
  PRINT, 'Io Plasma Torus using proper line-of-sight integration with:'
  PRINT, '  - Ray tracing through 3D plasma model'
  PRINT, '  - Trilinear interpolation of Te, ne, n_ion at each LOS point'
  PRINT, '  - Per-ion photon emission rates from CHIANTI tables'
  PRINT, '  - Simpson''s rule integration for brightness in Rayleighs'
  PRINT, '  - Gaussian instrumental response via analytic ERF'
  PRINT, ''

  ; ============================================================================
  ; LOAD DATA FILES
  ; ============================================================================

  PRINT, '======================================================================'
  PRINT, 'LOADING PLASMA MODEL AND EMISSION TABLES'
  PRINT, '======================================================================'
  PRINT, ''

  ; Load 3D plasma model
  ; The plasma model contains electron temperature, electron density, and
  ; ion number densities on a 3D Cartesian grid in Jupiter-centered coordinates
  plasma_filename = 'jovian_plasma_interpolated_381x381x231.h5'
  plasma_file = FILEPATH(plasma_filename, ROOT_DIR='..', SUBDIRECTORY=['..', '3D_Torus_Model'])
  plasma_model = IPT_LOAD_PLASMA_MODEL(plasma_file)
  PRINT, ''

  ; Load single Maxwellian emission tables
  ; Tables contain volume emissivity [photons/s/cm^3] as function of Te, ne
  ; for each emission line of each ion species
  emission_filename_single = 'CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5'
  emission_file_single = FILEPATH(emission_filename_single, ROOT_DIR='..', SUBDIRECTORY=['..', 'Emiss_tables'])
  tables_single = IPT_LOAD_EMISSION_TABLES_SINGLE(emission_file_single)
  PRINT, ''

  ; Load double Maxwellian emission tables
  ; Tables contain volume emissivity as function of Tec, Teh, ne_total, feh
  ; for modeling suprathermal electron populations
  emission_filename_double = 'CHIANTI_11.0.2_emiss_tables_double_maxwellian_16x8x12x10.h5'
  emission_file_double = FILEPATH(emission_filename_double, ROOT_DIR='..', SUBDIRECTORY=['..', 'Emiss_tables'])
  tables_double = IPT_LOAD_EMISSION_TABLES_DOUBLE(emission_file_double)
  PRINT, ''

  ; ============================================================================
  ; INSTRUMENT PARAMETERS - GROUND-BASED OPTICAL TELESCOPE
  ; ============================================================================
  ; These parameters are typical for ground-based optical spectroscopy
  ; of the Io Plasma Torus (e.g., Keck HIRES, VLT UVES, Subaru HDS)

  wavelength_range = [3000.0d, 10000.0d]  ; Angstroms (ground-based optical + near-IR)
  bin_width = 0.6d                         ; Angstroms (typical for R~3000 spectrograph)
  fwhm = 1.2d                              ; Angstroms (typical ground-based seeing-limited)
  ds = 0.1d                               ; Integration step size in R_J (714.92 km)

  PRINT, 'Instrument Parameters:'
  PRINT, '  Wavelength range: ', STRTRIM(STRING(wavelength_range[0], FORMAT='(F7.1)'),2), $
    ' - ', STRTRIM(STRING(wavelength_range[1], FORMAT='(F8.1)'),2), ' Angstrom'
  PRINT, '  Spectral bin width: ', STRTRIM(STRING(bin_width, FORMAT='(F5.2)'),2), ' Angstrom'
  PRINT, '  Instrumental FWHM: ', STRTRIM(STRING(fwhm, FORMAT='(F5.2)'),2), ' Angstrom'
  PRINT, '  LOS integration step: ', STRTRIM(STRING(ds, FORMAT='(F6.4)'),2), ' R_J'
  PRINT, ''

  ; ============================================================================
  ; TEST CASE 1: EQUATORIAL LINE OF SIGHT (SINGLE MAXWELLIAN)
  ; ============================================================================
  ; This geometry samples the densest part of the Io Plasma Torus,
  ; looking through the centrifugal equator where plasma density peaks

  PRINT, ''
  PRINT, '======================================================================'
  PRINT, 'TEST CASE 1: EQUATORIAL LINE OF SIGHT (SINGLE MAXWELLIAN)'
  PRINT, '----------------------------------------------------------------------'

  ; Define line of sight geometry
  ; Start position: 6 R_J from Jupiter, 20 R_J in -y direction, on equator
  ; Direction: Looking in +y direction through the torus
  slit_pos_vec1 = [6.0d, -20.0d, 0.0d]   ; Start at x=6 R_J, equator
  norm_vec1 = [0.0d, 1.0d, 0.0d]         ; Look in +y direction

  PRINT, 'Starting position: (', STRTRIM(STRING(slit_pos_vec1[0], FORMAT='(F5.1)'),2), ', ', $
    STRTRIM(STRING(slit_pos_vec1[1], FORMAT='(F6.1)'),2), ', ', $
    STRTRIM(STRING(slit_pos_vec1[2], FORMAT='(F5.1)'),2), ') R_J'
  PRINT, 'Direction vector: (', STRTRIM(STRING(norm_vec1[0], FORMAT='(F4.1)'),2), ', ', $
    STRTRIM(STRING(norm_vec1[1], FORMAT='(F4.1)'),2), ', ', $
    STRTRIM(STRING(norm_vec1[2], FORMAT='(F4.1)'),2), ')'
  PRINT, ''
  PRINT, 'Calculating optical spectrum using single Maxwellian...'
  PRINT, '  (Ray tracing with trilinear plasma interpolation)'
  PRINT, '  (Per-ion emission rates with Simpson''s rule integration)'

  ; Calculate spectrum with proper LOS integration
  ; This procedure:
  ; 1. Creates ray through plasma model
  ; 2. Interpolates Te, ne, n_ion at each point
  ; 3. Calculates per-ion emission rates from tables
  ; 4. Integrates brightness along LOS using Simpson's rule
  ; 5. Convolves with instrumental response
  IPT_CALCULATE_SPECTRUM_SINGLE, plasma_model, tables_single, $
    slit_pos_vec1, norm_vec1, wave_bins1_single, spectrum1_single, lines1_single, $
    WAVELENGTH_RANGE=wavelength_range, BIN_WIDTH=bin_width, FWHM=fwhm, DS=ds

  ; Calculate total integrated brightness
  total_brightness1_single = IPT_SIMPSON_INTEGRATE(wave_bins1_single, spectrum1_single)
  peak_idx1_single = WHERE(spectrum1_single EQ MAX(spectrum1_single))

  PRINT, ''
  PRINT, 'Single Maxwellian Results:'
  PRINT, '  Total integrated brightness: ', STRTRIM(STRING(total_brightness1_single, FORMAT='(F12.1)'),2), ' Rayleighs'
  PRINT, '  Peak spectral brightness: ', STRTRIM(STRING(MAX(spectrum1_single), FORMAT='(F10.2)'),2), $
    ' R/Angstrom at ', STRTRIM(STRING(wave_bins1_single[peak_idx1_single[0]], FORMAT='(F8.1)'),2), ' Angstrom'
  PRINT, '  Number of emission lines: ', STRTRIM(STRING(N_ELEMENTS(lines1_single.wavelength)),2)

  ; ============================================================================
  ; TEST CASE 2: OFF-EQUATOR LINE OF SIGHT (SINGLE MAXWELLIAN)
  ; ============================================================================
  ; This geometry samples plasma above the centrifugal equator,
  ; where ion densities decrease due to scale height effects

  PRINT, ''
  PRINT, '======================================================================'
  PRINT, 'TEST CASE 2: OFF-EQUATOR LINE OF SIGHT (SINGLE MAXWELLIAN)'
  PRINT, '----------------------------------------------------------------------'

  ; Define off-equator line of sight
  ; 0.5 R_J above the centrifugal equator
  slit_pos_vec2 = [6.0d, -20.0d, 0.5d]   ; 0.5 R_J above equator
  norm_vec2 = [0.0d, 1.0d, 0.0d]         ; Same direction

  PRINT, 'Starting position: (', STRTRIM(STRING(slit_pos_vec2[0], FORMAT='(F5.1)'),2), ', ', $
    STRTRIM(STRING(slit_pos_vec2[1], FORMAT='(F6.1)'),2), ', ', $
    STRTRIM(STRING(slit_pos_vec2[2], FORMAT='(F5.1)'),2), ') R_J'
  PRINT, 'Direction vector: (', STRTRIM(STRING(norm_vec2[0], FORMAT='(F4.1)'),2), ', ', $
    STRTRIM(STRING(norm_vec2[1], FORMAT='(F4.1)'),2), ', ', $
    STRTRIM(STRING(norm_vec2[2], FORMAT='(F4.1)'),2), ')'
  PRINT, ''
  PRINT, 'Calculating optical spectrum using single Maxwellian...'

  IPT_CALCULATE_SPECTRUM_SINGLE, plasma_model, tables_single, $
    slit_pos_vec2, norm_vec2, wave_bins2_single, spectrum2_single, lines2_single, $
    WAVELENGTH_RANGE=wavelength_range, BIN_WIDTH=bin_width, FWHM=fwhm, DS=ds

  total_brightness2_single = IPT_SIMPSON_INTEGRATE(wave_bins2_single, spectrum2_single)

  PRINT, ''
  PRINT, 'Single Maxwellian Results:'
  PRINT, '  Total integrated brightness: ', STRTRIM(STRING(total_brightness2_single, FORMAT='(F12.1)'),2), ' Rayleighs'
  PRINT, '  Equator/Off-equator ratio: ', $
    STRTRIM(STRING(total_brightness1_single/total_brightness2_single, FORMAT='(F6.2)'),2)

  ; ============================================================================
  ; TEST CASE 3: EQUATORIAL LINE OF SIGHT (DOUBLE MAXWELLIAN)
  ; ============================================================================
  ; Double Maxwellian includes hot electron population that enhances
  ; emission from higher ionization states and collisionally excited lines

  PRINT, ''
  PRINT, '======================================================================'
  PRINT, 'TEST CASE 3: EQUATORIAL LINE OF SIGHT (DOUBLE MAXWELLIAN)'
  PRINT, '----------------------------------------------------------------------'

  PRINT, 'Starting position: (', STRTRIM(STRING(slit_pos_vec1[0], FORMAT='(F5.1)'),2), ', ', $
    STRTRIM(STRING(slit_pos_vec1[1], FORMAT='(F6.1)'),2), ', ', $
    STRTRIM(STRING(slit_pos_vec1[2], FORMAT='(F5.1)'),2), ') R_J'
  PRINT, 'Direction vector: (', STRTRIM(STRING(norm_vec1[0], FORMAT='(F4.1)'),2), ', ', $
    STRTRIM(STRING(norm_vec1[1], FORMAT='(F4.1)'),2), ', ', $
    STRTRIM(STRING(norm_vec1[2], FORMAT='(F4.1)'),2), ')'
  PRINT, ''
  PRINT, 'Calculating optical spectrum using double Maxwellian...'
  PRINT, '  (Includes hot electron population from wave-particle heating)'
  PRINT, '  (Automatic fallback to single Maxwellian for low feh regions)'

  ; Calculate spectrum with double Maxwellian tables
  ; Uses 4D interpolation when hot electron fraction is significant
  ; Falls back to single Maxwellian for regions with low feh
  IPT_CALCULATE_SPECTRUM_DOUBLE, plasma_model, tables_single, tables_double, $
    slit_pos_vec1, norm_vec1, wave_bins1_double, spectrum1_double, lines1_double, $
    WAVELENGTH_RANGE=wavelength_range, BIN_WIDTH=bin_width, FWHM=fwhm, DS=ds

  total_brightness1_double = IPT_SIMPSON_INTEGRATE(wave_bins1_double, spectrum1_double)
  peak_idx1_double = WHERE(spectrum1_double EQ MAX(spectrum1_double))

  PRINT, ''
  PRINT, 'Double Maxwellian Results:'
  PRINT, '  Total integrated brightness: ', STRTRIM(STRING(total_brightness1_double, FORMAT='(F12.1)'),2), ' Rayleighs'
  PRINT, '  Peak spectral brightness: ', STRTRIM(STRING(MAX(spectrum1_double), FORMAT='(F10.2)'),2), $
    ' R/Angstrom at ', STRTRIM(STRING(wave_bins1_double[peak_idx1_double[0]], FORMAT='(F8.1)'),2), ' Angstrom'
  PRINT, '  Hot electron enhancement factor: ', $
    STRTRIM(STRING(total_brightness1_double/total_brightness1_single, FORMAT='(F6.2)'),2)

  ; ============================================================================
  ; TEST CASE 4: OFF-EQUATOR LINE OF SIGHT (DOUBLE MAXWELLIAN)
  ; ============================================================================

  PRINT, ''
  PRINT, '======================================================================'
  PRINT, 'TEST CASE 4: OFF-EQUATOR LINE OF SIGHT (DOUBLE MAXWELLIAN)'
  PRINT, '----------------------------------------------------------------------'

  PRINT, 'Starting position: (', STRTRIM(STRING(slit_pos_vec2[0], FORMAT='(F5.1)'),2), ', ', $
    STRTRIM(STRING(slit_pos_vec2[1], FORMAT='(F6.1)'),2), ', ', $
    STRTRIM(STRING(slit_pos_vec2[2], FORMAT='(F5.1)'),2), ') R_J'
  PRINT, 'Direction vector: (', STRTRIM(STRING(norm_vec2[0], FORMAT='(F4.1)'),2), ', ', $
    STRTRIM(STRING(norm_vec2[1], FORMAT='(F4.1)'),2), ', ', $
    STRTRIM(STRING(norm_vec2[2], FORMAT='(F4.1)'),2), ')'
  PRINT, ''
  PRINT, 'Calculating optical spectrum using double Maxwellian...'

  IPT_CALCULATE_SPECTRUM_DOUBLE, plasma_model, tables_single, tables_double, $
    slit_pos_vec2, norm_vec2, wave_bins2_double, spectrum2_double, lines2_double, $
    WAVELENGTH_RANGE=wavelength_range, BIN_WIDTH=bin_width, FWHM=fwhm, DS=ds

  total_brightness2_double = IPT_SIMPSON_INTEGRATE(wave_bins2_double, spectrum2_double)

  PRINT, ''
  PRINT, 'Double Maxwellian Results:'
  PRINT, '  Total integrated brightness: ', STRTRIM(STRING(total_brightness2_double, FORMAT='(F12.1)'),2), ' Rayleighs'
  PRINT, '  Hot electron enhancement factor: ', $
    STRTRIM(STRING(total_brightness2_double/total_brightness2_single, FORMAT='(F6.2)'),2)

  ; ============================================================================
  ; SUMMARY TABLE
  ; ============================================================================

  PRINT, ''
  PRINT, '======================================================================'
  PRINT, 'BRIGHTNESS SUMMARY'
  PRINT, '======================================================================'
  PRINT, ''
  PRINT, STRING('Geometry', FORMAT='(A25)'), $
    STRING('Single [R]', FORMAT='(A15)'), $
    STRING('Double [R]', FORMAT='(A15)'), $
    STRING('Enhancement', FORMAT='(A12)')
  PRINT, REPLICATE('-', 67)
  PRINT, STRING('Equatorial (z=0)', FORMAT='(A25)'), $
    STRING(total_brightness1_single, FORMAT='(F15.1)'), $
    STRING(total_brightness1_double, FORMAT='(F15.1)'), $
    STRING(total_brightness1_double/total_brightness1_single, FORMAT='(F12.2)')
  PRINT, STRING('Off-equator (z=0.5 R_J)', FORMAT='(A25)'), $
    STRING(total_brightness2_single, FORMAT='(F15.1)'), $
    STRING(total_brightness2_double, FORMAT='(F15.1)'), $
    STRING(total_brightness2_double/total_brightness2_single, FORMAT='(F12.2)')
  PRINT, REPLICATE('=', 67)

  ; ============================================================================
  ; MATPLOTLIB COLORS FOR PLOTS
  ; ============================================================================
  ; Define colors matching matplotlib's default color cycle for consistency
  ; with Python-based analysis tools and publications
  mp_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', $
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

  ; ============================================================================
  ; PLOT 1: EQUATORIAL SINGLE VS DOUBLE MAXWELLIAN
  ; ============================================================================

  PRINT, ''
  PRINT, '======================================================================'
  PRINT, 'CREATING PUBLICATION-QUALITY PLOTS'
  PRINT, '======================================================================'

  p1 = PLOT(wave_bins1_single, spectrum1_single, $
    COLOR=mp_colors[0], THICK=2, $
    XTITLE='Wavelength [Angstrom]', $
    YTITLE='Brightness [R/Angstrom]', $
    TITLE='IPT Optical Emission: Equatorial LOS (Proper Integration)', $
    XRANGE=wavelength_range, $
    YRANGE=[0, MAX([MAX(spectrum1_single), MAX(spectrum1_double)])*1.1], $
    FONT_SIZE=12, DIMENSIONS=[1200, 600], $
    NAME='Single Maxwellian')

  p1_double = PLOT(wave_bins1_double, spectrum1_double, /OVERPLOT, $
    COLOR=mp_colors[1], THICK=2, TRANSPARENCY=30, $
    NAME='Double Maxwellian')

  leg1 = LEGEND(TARGET=[p1, p1_double], $
    /RELATIVE, POSITION=[0.85, 0.85], $
    HORIZONTAL_ALIGNMENT='right', VERTICAL_ALIGNMENT='top', $
    /AUTO_TEXT_COLOR, FONT_SIZE=11)

  p1.XGRIDSTYLE = 1
  p1.YGRIDSTYLE = 1

  p1.SAVE, 'ipt_optical_emission_integration_equatorial_comparison.png', RESOLUTION=300
  PRINT, '  Saved: ipt_optical_emission_integration_equatorial_comparison.png'

  ; ============================================================================
  ; PLOT 2: OFF-EQUATOR SINGLE VS DOUBLE MAXWELLIAN
  ; ============================================================================

  p2 = PLOT(wave_bins2_single, spectrum2_single, $
    COLOR=mp_colors[0], THICK=2, $
    XTITLE='Wavelength [Angstrom]', $
    YTITLE='Brightness [R/Angstrom]', $
    TITLE='IPT Optical Emission: Off-Equator LOS (z = 0.5 R_J)', $
    XRANGE=wavelength_range, $
    YRANGE=[0, MAX([MAX(spectrum2_single), MAX(spectrum2_double)])*1.1], $
    FONT_SIZE=12, DIMENSIONS=[1200, 600], $
    NAME='Single Maxwellian')

  p2_double = PLOT(wave_bins2_double, spectrum2_double, /OVERPLOT, $
    COLOR=mp_colors[1], THICK=2, TRANSPARENCY=30, $
    NAME='Double Maxwellian')

  leg2 = LEGEND(TARGET=[p2, p2_double], $
    /RELATIVE, POSITION=[0.85, 0.85], $
    HORIZONTAL_ALIGNMENT='right', VERTICAL_ALIGNMENT='top', $
    /AUTO_TEXT_COLOR, FONT_SIZE=11)

  p2.XGRIDSTYLE = 1
  p2.YGRIDSTYLE = 1

  p2.SAVE, 'ipt_optical_emission_integration_offequator_comparison.png', RESOLUTION=300
  PRINT, '  Saved: ipt_optical_emission_integration_offequator_comparison.png'

  ; ============================================================================
  ; PLOT 3: FOUR-PANEL COMPARISON
  ; ============================================================================

  w = WINDOW(DIMENSIONS=[1400, 1400])

  ; Panel 1: Equatorial single Maxwellian
  p3a = PLOT(wave_bins1_single, spectrum1_single, $
    COLOR=mp_colors[0], THICK=2, $
    XTITLE='Wavelength [Angstrom]', $
    YTITLE='Brightness [R/Angstrom]', $
    TITLE='Equatorial LOS: Single Maxwellian', $
    XRANGE=wavelength_range, $
    YRANGE=[0, MAX(spectrum1_single)*1.1], $
    FONT_SIZE=10, LAYOUT=[2,2,1], /CURRENT)
  p3a.XGRIDSTYLE = 1
  p3a.YGRIDSTYLE = 1

  ; Panel 2: Equatorial double Maxwellian
  p3b = PLOT(wave_bins1_double, spectrum1_double, $
    COLOR=mp_colors[1], THICK=2, $
    XTITLE='Wavelength [Angstrom]', $
    YTITLE='Brightness [R/Angstrom]', $
    TITLE='Equatorial LOS: Double Maxwellian', $
    XRANGE=wavelength_range, $
    YRANGE=[0, MAX(spectrum1_double)*1.1], $
    FONT_SIZE=10, LAYOUT=[2,2,2], /CURRENT)
  p3b.XGRIDSTYLE = 1
  p3b.YGRIDSTYLE = 1

  ; Panel 3: Off-equator single Maxwellian
  p3c = PLOT(wave_bins2_single, spectrum2_single, $
    COLOR=mp_colors[0], THICK=2, $
    XTITLE='Wavelength [Angstrom]', $
    YTITLE='Brightness [R/Angstrom]', $
    TITLE='Off-Equator LOS (z=0.5 R_J): Single Maxwellian', $
    XRANGE=wavelength_range, $
    YRANGE=[0, MAX(spectrum2_single)*1.1], $
    FONT_SIZE=10, LAYOUT=[2,2,3], /CURRENT)
  p3c.XGRIDSTYLE = 1
  p3c.YGRIDSTYLE = 1

  ; Panel 4: Off-equator double Maxwellian
  p3d = PLOT(wave_bins2_double, spectrum2_double, $
    COLOR=mp_colors[1], THICK=2, $
    XTITLE='Wavelength [Angstrom]', $
    YTITLE='Brightness [R/Angstrom]', $
    TITLE='Off-Equator LOS (z=0.5 R_J): Double Maxwellian', $
    XRANGE=wavelength_range, $
    YRANGE=[0, MAX(spectrum2_double)*1.1], $
    FONT_SIZE=10, LAYOUT=[2,2,4], /CURRENT)
  p3d.XGRIDSTYLE = 1
  p3d.YGRIDSTYLE = 1

  w.SAVE, 'ipt_optical_emission_integration_four_panel.png', RESOLUTION=300
  PRINT, '  Saved: ipt_optical_emission_integration_four_panel.png'

  ; ============================================================================
  ; PLOT 4: [S II] 6716/6731 DOUBLET REGION (DENSITY DIAGNOSTIC)
  ; ============================================================================
  ; This doublet is one of the most important density diagnostics for
  ; ground-based observations of the Io Plasma Torus

  mask_sii = WHERE(wave_bins1_single GE 6680 AND wave_bins1_single LE 6760, count_sii)

  IF count_sii GT 0 THEN BEGIN
    p4 = PLOT(wave_bins1_single[mask_sii], spectrum1_single[mask_sii], $
      COLOR=mp_colors[0], THICK=2, $
      XTITLE='Wavelength [Angstrom]', $
      YTITLE='Brightness [R/Angstrom]', $
      TITLE='[S II] 6716/6731 Angstrom Doublet (Density Diagnostic)', $
      XRANGE=[6680, 6760], $
      FONT_SIZE=12, DIMENSIONS=[1000, 600], $
      NAME='Single Maxwellian')

    p4_double = PLOT(wave_bins1_double[mask_sii], spectrum1_double[mask_sii], $
      /OVERPLOT, COLOR=mp_colors[1], THICK=2, TRANSPARENCY=30, $
      NAME='Double Maxwellian')

    leg4 = LEGEND(TARGET=[p4, p4_double], $
      /RELATIVE, POSITION=[0.85, 0.85], $
      HORIZONTAL_ALIGNMENT='right', VERTICAL_ALIGNMENT='top', $
      /AUTO_TEXT_COLOR, FONT_SIZE=11)

    p4.XGRIDSTYLE = 1
    p4.YGRIDSTYLE = 1

    p4.SAVE, 'ipt_optical_sii_doublet_region.png', RESOLUTION=300
    PRINT, '  Saved: ipt_optical_sii_doublet_region.png'
  ENDIF

  ; ============================================================================
  ; PLOT 5: [O III] 4959/5007 REGION (NEBULAR LINES)
  ; ============================================================================
  ; The [O III] nebular lines are strong diagnostics for O++ abundance
  ; and electron temperature

  mask_oiii = WHERE(wave_bins1_single GE 4920 AND wave_bins1_single LE 5040, count_oiii)

  IF count_oiii GT 0 THEN BEGIN
    p5 = PLOT(wave_bins1_single[mask_oiii], spectrum1_single[mask_oiii], $
      COLOR=mp_colors[0], THICK=2, $
      XTITLE='Wavelength [Angstrom]', $
      YTITLE='Brightness [R/Angstrom]', $
      TITLE='[O III] 4959/5007 Angstrom (Nebular Lines)', $
      XRANGE=[4920, 5040], $
      FONT_SIZE=12, DIMENSIONS=[1000, 600], $
      NAME='Single Maxwellian')

    p5_double = PLOT(wave_bins1_double[mask_oiii], spectrum1_double[mask_oiii], $
      /OVERPLOT, COLOR=mp_colors[1], THICK=2, TRANSPARENCY=30, $
      NAME='Double Maxwellian')

    leg5 = LEGEND(TARGET=[p5, p5_double], $
      /RELATIVE, POSITION=[0.85, 0.85], $
      HORIZONTAL_ALIGNMENT='right', VERTICAL_ALIGNMENT='top', $
      /AUTO_TEXT_COLOR, FONT_SIZE=11)

    p5.XGRIDSTYLE = 1
    p5.YGRIDSTYLE = 1

    p5.SAVE, 'ipt_optical_oiii_nebular_region.png', RESOLUTION=300
    PRINT, '  Saved: ipt_optical_oiii_nebular_region.png'
  ENDIF

  ; ============================================================================
  ; PLOT 6: [S III] 6312 REGION (TEMPERATURE DIAGNOSTIC)
  ; ============================================================================

  mask_siii = WHERE(wave_bins1_single GE 6280 AND wave_bins1_single LE 6340, count_siii)

  IF count_siii GT 0 THEN BEGIN
    p6 = PLOT(wave_bins1_single[mask_siii], spectrum1_single[mask_siii], $
      COLOR=mp_colors[0], THICK=2, $
      XTITLE='Wavelength [Angstrom]', $
      YTITLE='Brightness [R/Angstrom]', $
      TITLE='[S III] 6312 Angstrom (Temperature Diagnostic)', $
      XRANGE=[6280, 6340], $
      FONT_SIZE=12, DIMENSIONS=[1000, 600], $
      NAME='Single Maxwellian')

    p6_double = PLOT(wave_bins1_double[mask_siii], spectrum1_double[mask_siii], $
      /OVERPLOT, COLOR=mp_colors[1], THICK=2, TRANSPARENCY=30, $
      NAME='Double Maxwellian')

    leg6 = LEGEND(TARGET=[p6, p6_double], $
      /RELATIVE, POSITION=[0.85, 0.85], $
      HORIZONTAL_ALIGNMENT='right', VERTICAL_ALIGNMENT='top', $
      /AUTO_TEXT_COLOR, FONT_SIZE=11)

    p6.XGRIDSTYLE = 1
    p6.YGRIDSTYLE = 1

    p6.SAVE, 'ipt_optical_siii_6312_region.png', RESOLUTION=300
    PRINT, '  Saved: ipt_optical_siii_6312_region.png'
  ENDIF

  ; ============================================================================
  ; KEY OPTICAL DIAGNOSTIC LINES ANALYSIS
  ; ============================================================================

  PRINT, ''
  PRINT, '======================================================================'
  PRINT, 'KEY OPTICAL DIAGNOSTIC LINES ANALYSIS'
  PRINT, '======================================================================'
  PRINT, 'Equatorial Line of Sight Results'
  PRINT, 'Instrumental FWHM: ', STRTRIM(STRING(fwhm, FORMAT='(F5.2)'),2), ' Angstrom'
  PRINT, ''

  ; Define key optical emission lines for analysis
  ; These are important diagnostics from ground-based observations
  line_names = ['[O II] 3726.0', '[O II] 3728.8', '[S II] 4068.6', '[S II] 4076.3', $
    '[O III] 4363.2', '[O III] 4958.9', '[O III] 5006.8', '[S III] 6312.1', $
    '[S II] 6716.4', '[S II] 6730.8', '[S III] 9068.6']
  line_wavelengths = [3726.03d, 3728.82d, 4068.6d, 4076.3d, 4363.21d, 4958.91d, $
    5006.84d, 6312.1d, 6716.4d, 6730.8d, 9068.6d]
  line_species = ['O+', 'O+', 'S+', 'S+', 'O++', 'O++', 'O++', 'S++', 'S+', 'S+', 'S++']

  ; Convert FWHM to Gaussian sigma for integration window
  sigma = fwhm / 2.35482d

  PRINT, STRING('Line', FORMAT='(A20)'), $
    STRING('Center [Ang]', FORMAT='(A12)'), $
    STRING('Single [R]', FORMAT='(A15)'), $
    STRING('Double [R]', FORMAT='(A15)'), $
    STRING('Enhancement', FORMAT='(A12)')
  PRINT, REPLICATE('-', 74)

  FOR i = 0, N_ELEMENTS(line_wavelengths)-1 DO BEGIN
    target_wav = line_wavelengths[i]

    ; Define integration region (±3 sigma from line center)
    wav_min = target_wav - 3.0d * sigma
    wav_max = target_wav + 3.0d * sigma

    ; Find indices within integration range
    idx = WHERE(wave_bins1_single GE wav_min AND wave_bins1_single LE wav_max, count)

    IF count GT 2 THEN BEGIN
      ; Integrate brightness over line profile using Simpson's rule
      int_single = IPT_SIMPSON_INTEGRATE(wave_bins1_single[idx], spectrum1_single[idx])
      int_double = IPT_SIMPSON_INTEGRATE(wave_bins1_double[idx], spectrum1_double[idx])

      ; Calculate enhancement factor from hot electrons
      IF int_single GT 0 THEN BEGIN
        enhancement = int_double / int_single
      ENDIF ELSE BEGIN
        enhancement = 0.0d
      ENDELSE

      PRINT, STRING(line_names[i], FORMAT='(A20)'), $
        STRING(target_wav, FORMAT='(F12.3)'), $
        STRING(int_single, FORMAT='(F15.1)'), $
        STRING(int_double, FORMAT='(F15.1)'), $
        STRING(enhancement, FORMAT='(F12.3)')
    ENDIF ELSE BEGIN
      PRINT, STRING(line_names[i], FORMAT='(A20)'), '  (outside wavelength range or insufficient data)'
    ENDELSE
  ENDFOR

  PRINT, REPLICATE('=', 74)

  ; ============================================================================
  ; OPTICAL LINE RATIO DIAGNOSTICS
  ; ============================================================================

  PRINT, ''
  PRINT, '======================================================================'
  PRINT, 'OPTICAL LINE RATIO DIAGNOSTICS'
  PRINT, '======================================================================'
  PRINT, ''
  PRINT, 'These ratios are commonly used for plasma diagnostics in ground-based'
  PRINT, 'observations of the Io Plasma Torus.'
  PRINT, ''

  ; [S II] 6716/6731 ratio (density sensitive)
  ; This doublet ratio varies from ~1.5 at low density to ~0.44 at high density
  idx_6716 = WHERE(wave_bins1_single GE 6716.4d - 3.0d*sigma AND $
    wave_bins1_single LE 6716.4d + 3.0d*sigma, count_6716)
  idx_6731 = WHERE(wave_bins1_single GE 6730.8d - 3.0d*sigma AND $
    wave_bins1_single LE 6730.8d + 3.0d*sigma, count_6731)

  IF count_6716 GT 2 AND count_6731 GT 2 THEN BEGIN
    int_6716_single = IPT_SIMPSON_INTEGRATE(wave_bins1_single[idx_6716], spectrum1_single[idx_6716])
    int_6731_single = IPT_SIMPSON_INTEGRATE(wave_bins1_single[idx_6731], spectrum1_single[idx_6731])
    IF int_6731_single GT 0 THEN ratio_sii_single = int_6716_single / int_6731_single

    int_6716_double = IPT_SIMPSON_INTEGRATE(wave_bins1_double[idx_6716], spectrum1_double[idx_6716])
    int_6731_double = IPT_SIMPSON_INTEGRATE(wave_bins1_double[idx_6731], spectrum1_double[idx_6731])
    IF int_6731_double GT 0 THEN ratio_sii_double = int_6716_double / int_6731_double

    PRINT, '[S II] 6716/6731 ratio (density diagnostic):'
    PRINT, '  Single Maxwellian: ', STRTRIM(STRING(ratio_sii_single, FORMAT='(F8.4)'),2)
    PRINT, '  Double Maxwellian: ', STRTRIM(STRING(ratio_sii_double, FORMAT='(F8.4)'),2)
    PRINT, '  (Low density limit ~1.5, high density limit ~0.44)'
    PRINT, ''
  ENDIF

  ; [O II] 3726/3729 ratio (density sensitive)
  idx_3726 = WHERE(wave_bins1_single GE 3726.03d - 3.0d*sigma AND $
    wave_bins1_single LE 3726.03d + 3.0d*sigma, count_3726)
  idx_3729 = WHERE(wave_bins1_single GE 3728.82d - 3.0d*sigma AND $
    wave_bins1_single LE 3728.82d + 3.0d*sigma, count_3729)

  IF count_3726 GT 2 AND count_3729 GT 2 THEN BEGIN
    int_3726_single = IPT_SIMPSON_INTEGRATE(wave_bins1_single[idx_3726], spectrum1_single[idx_3726])
    int_3729_single = IPT_SIMPSON_INTEGRATE(wave_bins1_single[idx_3729], spectrum1_single[idx_3729])
    IF int_3729_single GT 0 THEN ratio_oii_single = int_3726_single / int_3729_single

    int_3726_double = IPT_SIMPSON_INTEGRATE(wave_bins1_double[idx_3726], spectrum1_double[idx_3726])
    int_3729_double = IPT_SIMPSON_INTEGRATE(wave_bins1_double[idx_3729], spectrum1_double[idx_3729])
    IF int_3729_double GT 0 THEN ratio_oii_double = int_3726_double / int_3729_double

    PRINT, '[O II] 3726/3729 ratio (density diagnostic):'
    PRINT, '  Single Maxwellian: ', STRTRIM(STRING(ratio_oii_single, FORMAT='(F8.4)'),2)
    PRINT, '  Double Maxwellian: ', STRTRIM(STRING(ratio_oii_double, FORMAT='(F8.4)'),2)
    PRINT, ''
  ENDIF

  ; [O III] temperature diagnostic ratio
  ; 4363/(4959+5007) is sensitive to electron temperature
  idx_4363 = WHERE(wave_bins1_single GE 4363.21d - 3.0d*sigma AND $
    wave_bins1_single LE 4363.21d + 3.0d*sigma, count_4363)
  idx_4959 = WHERE(wave_bins1_single GE 4958.91d - 3.0d*sigma AND $
    wave_bins1_single LE 4958.91d + 3.0d*sigma, count_4959)
  idx_5007 = WHERE(wave_bins1_single GE 5006.84d - 3.0d*sigma AND $
    wave_bins1_single LE 5006.84d + 3.0d*sigma, count_5007)

  IF count_4363 GT 2 AND count_4959 GT 2 AND count_5007 GT 2 THEN BEGIN
    int_4363_single = IPT_SIMPSON_INTEGRATE(wave_bins1_single[idx_4363], spectrum1_single[idx_4363])
    int_4959_single = IPT_SIMPSON_INTEGRATE(wave_bins1_single[idx_4959], spectrum1_single[idx_4959])
    int_5007_single = IPT_SIMPSON_INTEGRATE(wave_bins1_single[idx_5007], spectrum1_single[idx_5007])
    IF (int_4959_single + int_5007_single) GT 0 THEN $
      ratio_oiii_single = int_4363_single / (int_4959_single + int_5007_single)

    int_4363_double = IPT_SIMPSON_INTEGRATE(wave_bins1_double[idx_4363], spectrum1_double[idx_4363])
    int_4959_double = IPT_SIMPSON_INTEGRATE(wave_bins1_double[idx_4959], spectrum1_double[idx_4959])
    int_5007_double = IPT_SIMPSON_INTEGRATE(wave_bins1_double[idx_5007], spectrum1_double[idx_5007])
    IF (int_4959_double + int_5007_double) GT 0 THEN $
      ratio_oiii_double = int_4363_double / (int_4959_double + int_5007_double)

    PRINT, '[O III] 4363/(4959+5007) ratio (temperature diagnostic):'
    PRINT, '  Single Maxwellian: ', STRTRIM(STRING(ratio_oiii_single, FORMAT='(F8.6)'),2)
    PRINT, '  Double Maxwellian: ', STRTRIM(STRING(ratio_oiii_double, FORMAT='(F8.6)'),2)
    PRINT, '  (Higher ratio indicates higher effective temperature)'
    PRINT, ''
  ENDIF

  ; [S III] temperature diagnostic
  ; 6312/9069 ratio is temperature sensitive
  idx_6312 = WHERE(wave_bins1_single GE 6312.1d - 3.0d*sigma AND $
    wave_bins1_single LE 6312.1d + 3.0d*sigma, count_6312)
  idx_9069 = WHERE(wave_bins1_single GE 9068.6d - 3.0d*sigma AND $
    wave_bins1_single LE 9068.6d + 3.0d*sigma, count_9069)

  IF count_6312 GT 2 AND count_9069 GT 2 THEN BEGIN
    int_6312_single = IPT_SIMPSON_INTEGRATE(wave_bins1_single[idx_6312], spectrum1_single[idx_6312])
    int_9069_single = IPT_SIMPSON_INTEGRATE(wave_bins1_single[idx_9069], spectrum1_single[idx_9069])
    IF int_9069_single GT 0 THEN ratio_siii_single = int_6312_single / int_9069_single

    int_6312_double = IPT_SIMPSON_INTEGRATE(wave_bins1_double[idx_6312], spectrum1_double[idx_6312])
    int_9069_double = IPT_SIMPSON_INTEGRATE(wave_bins1_double[idx_9069], spectrum1_double[idx_9069])
    IF int_9069_double GT 0 THEN ratio_siii_double = int_6312_double / int_9069_double

    PRINT, '[S III] 6312/9069 ratio (temperature diagnostic):'
    PRINT, '  Single Maxwellian: ', STRTRIM(STRING(ratio_siii_single, FORMAT='(F8.4)'),2)
    PRINT, '  Double Maxwellian: ', STRTRIM(STRING(ratio_siii_double, FORMAT='(F8.4)'),2)
  ENDIF

  PRINT, REPLICATE('=', 74)

  ; ============================================================================
  ; EXECUTION TIME AND SUMMARY
  ; ============================================================================

  end_time = SYSTIME(/SECONDS)
  elapsed_time = end_time - start_time

  PRINT, ''
  PRINT, '======================================================================'
  PRINT, 'EXECUTION SUMMARY'
  PRINT, '======================================================================'
  PRINT, 'Total execution time: ', STRTRIM(STRING(elapsed_time, FORMAT='(F10.2)'),2), ' seconds'
  PRINT, ''
  PRINT, 'Output files created:'
  PRINT, '  - ipt_optical_emission_integration_equatorial_comparison.png'
  PRINT, '  - ipt_optical_emission_integration_offequator_comparison.png'
  PRINT, '  - ipt_optical_emission_integration_four_panel.png'
  PRINT, '  - ipt_optical_sii_doublet_region.png'
  PRINT, '  - ipt_optical_oiii_nebular_region.png'
  PRINT, '  - ipt_optical_siii_6312_region.png'
  PRINT, ''
  PRINT, '======================================================================'
  PRINT, 'Variables available for inspection:'
  PRINT, '  plasma_model          - 3D plasma model structure'
  PRINT, '  tables_single/double  - Emission table structures'
  PRINT, '  wave_bins*_single     - Single Maxwellian wavelength grids [Angstrom]'
  PRINT, '  wave_bins*_double     - Double Maxwellian wavelength grids [Angstrom]'
  PRINT, '  spectrum*_single      - Single Maxwellian spectra [R/Angstrom]'
  PRINT, '  spectrum*_double      - Double Maxwellian spectra [R/Angstrom]'
  PRINT, '  lines*_single/double  - Emission line list structures'
  PRINT, '======================================================================'

  ; Stop execution to allow interactive inspection of results
  STOP

END
