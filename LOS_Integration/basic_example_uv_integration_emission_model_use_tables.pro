;+
; basic_example_uv_integration_emission_model_use_tables.pro
;
; Example Script for UV Emission Line-of-Sight Integration Through Io Plasma Torus
; ================================================================================
;
; This script demonstrates UV emission calculations through the Io Plasma Torus
; using proper line-of-sight integration with species-by-species interpolation
; from CHIANTI emission tables. The code handles both single and double Maxwellian
; electron distributions with realistic test geometries.
;
; PHYSICAL MODEL:
; - Ray tracing along line of sight through 3D plasma model
; - Trilinear interpolation of plasma parameters (Te, ne, n_ion) at each LOS point
; - Bilinear interpolation of single Maxwellian photon emission rates (2D tables)
; - Quadlinear interpolation of double Maxwellian photon emission rates (4D tables)
; - Per-ion photon emission rate calculation: emissivity / n_e [photons/s/ion]
; - Simpson's rule integration over line of sight for brightness in Rayleighs
; - Gaussian instrumental response convolution via analytic ERF formulation
; - Output brightness per wavelength: Rayleighs/Angstrom
;
; BRIGHTNESS CALCULATION METHOD:
; - At each LOS point: interpolate Te, ne, n_ion from 3D plasma model
; - Interpolate volume emissivity from CHIANTI tables [photons/s/cm^3]
; - Calculate per-ion emission rate: volume_emiss / n_e [photons/s/ion]
; - Brightness integrand: per_ion_rate * n_ion [photons/s/cm^3]
; - Integrate over LOS with Simpson's rule
; - Convert to Rayleighs: integral * R_J_cm * 1e-6
;
; TEST CASES:
; 1. Equatorial line of sight (single Maxwellian)
; 2. Off-equator line of sight (single Maxwellian)
; 3. Equatorial line of sight (double Maxwellian)
; 4. Off-equator line of sight (double Maxwellian)
;
; All lines of sight start at x=6 R_J and look through the plasma torus
; in the +y direction, sampling the peak emission region near Io's orbit.
;
; REQUIRED DATA FILES:
; - Plasma Model: ../3D_Torus_Model/jovian_plasma_interpolated_381x381x231.h5
; - Single Maxwellian: ../Emiss_tables/CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5
; - Double Maxwellian: ../Emiss_tables/CHIANTI_11.0.2_emiss_tables_double_maxwellian_24x10x24x12.h5
;
; OUTPUT:
; - UV emission spectra for single and double Maxwellian cases
; - Comparison plots showing hot electron enhancement
; - Analysis of key UV diagnostic lines
; - Publication-quality plots saved as PNG files
;
; WAVELENGTH RANGE:
; - UV: 550-2100 Angstroms (Europa-UVS, JUICE-UVS, HST/STIS)
;
; COORDINATE SYSTEMS AND UNITS:
; - Positions: Jupiter radii [R_J], System III right-handed Cartesian
; - Temperatures: electron volts [eV]
; - Densities: particles per cubic centimeter [cm^-3]
; - Wavelengths: Angstroms [Angstroms]
; - Brightnesses: Rayleighs [R], where 1 R = 10^6 photons s^-1 cm^-2 (4pi sr)^-1
;
; REFERENCES:
; - CHIANTI database: Dere et al. 1997; Del Zanna et al. 2020; Dufresne et al. 2024
; - IPT observations: Steffl et al. 2004a,b; Thomas et al. 2004; Bagenal & Delamere 2011
; - Emission modeling: Nerney et al. 2017, 2020, 2022, 2025a, 2025b
;
; AUTHOR: Edward (Eddie) G. Nerney
; INSTITUTION: Laboratory for Atmospheric and Space Physics, University of Colorado Boulder
; LICENSE: Open source for academic and research use - Magnetospheres of Outer Planets (MOP) community
; VERSION: 2.0
; DATE: November 2025
;
; CHIANTI ACKNOWLEDGMENT:
; CHIANTI is a collaborative project involving George Mason University, the University
; of Michigan (USA), University of Cambridge (UK), and NASA Goddard Space Flight Center (USA).
;-

; Load the IPT emission raytracer library
@IPT_emiss_MOP_community_code

PRO basic_example_uv_integration_emission_model_use_tables
  ;+
  ; NAME:
  ;   basic_example_uv_integration_emission_model_use_tables
  ;
  ; PURPOSE:
  ;   Main procedure demonstrating UV emission calculations through the Io Plasma
  ;   Torus using line-of-sight integration with proper species-by-species
  ;   interpolation from CHIANTI tables. Calculates per-ion photon emission rates
  ;   and integrates over LOS to obtain brightness in Rayleighs.
  ;-

  ; Record start time
  start_time = SYSTIME(/SECONDS)

  PRINT, '======================================================================'
  PRINT, 'JOVIAN UV EMISSION RAYTRACER - LINE-OF-SIGHT INTEGRATION'
  PRINT, 'Io Plasma Torus UV Emission Model'
  PRINT, 'MOP Community Code v2.0'
  PRINT, '======================================================================'
  PRINT, ''
  PRINT, 'This example demonstrates UV emission calculations through the'
  PRINT, 'Io Plasma Torus using proper line-of-sight integration with'
  PRINT, 'species-by-species interpolation from CHIANTI tables.'
  PRINT, ''
  PRINT, 'Method: Per-ion photon emission rates integrated over LOS'
  PRINT, '  1. Interpolate Te, ne, n_ion from 3D model along ray'
  PRINT, '  2. Interpolate volume emissivity from CHIANTI tables'
  PRINT, '  3. Calculate per-ion rate: emissivity / n_e'
  PRINT, '  4. Brightness integrand: per_ion_rate * n_ion'
  PRINT, '  5. Simpson integration -> Rayleighs'
  PRINT, ''

  ; ============================================================================
  ; LOAD DATA FILES
  ; ============================================================================

  PRINT, '======================================================================'
  PRINT, 'LOADING PLASMA MODEL AND EMISSION TABLES'
  PRINT, '======================================================================'
  PRINT, ''

  ; Load 3D plasma model
  plasma_filename = 'jovian_plasma_interpolated_381x381x231.h5'
  plasma_file = FILEPATH(plasma_filename, ROOT_DIR='..', SUBDIRECTORY='3D_Torus_Model')
  plasma_model = IPT_LOAD_PLASMA_MODEL(plasma_file)
  PRINT, ''

  ; Load single Maxwellian emission tables
  emission_filename_single = 'CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5'
  emission_file_single = FILEPATH(emission_filename_single, ROOT_DIR='..', SUBDIRECTORY='Emiss_tables')
  tables_single = IPT_LOAD_EMISSION_TABLES_SINGLE(emission_file_single)
  PRINT, ''

  ; Load double Maxwellian emission tables
  emission_filename_double = 'CHIANTI_11.0.2_emiss_tables_double_maxwellian_24x10x24x12.h5'
  emission_file_double = FILEPATH(emission_filename_double, ROOT_DIR='..', SUBDIRECTORY='Emiss_tables')
  tables_double = IPT_LOAD_EMISSION_TABLES_DOUBLE(emission_file_double)
  PRINT, ''

  ; ============================================================================
  ; INSTRUMENT AND INTEGRATION PARAMETERS
  ; ============================================================================
  wavelength_range = [550.0d, 2100.0d]  ; Angstroms (Europa-UVS, JUICE-UVS range)
  bin_width = 1.0d                       ; Angstroms (spectral bin width)
  fwhm = 6.0d                            ; Angstroms (instrumental FWHM)
  ds = 0.01d                             ; Integration step size in R_J

  PRINT, '======================================================================'
  PRINT, 'INTEGRATION PARAMETERS'
  PRINT, '======================================================================'
  PRINT, '  Wavelength range: ', STRTRIM(STRING(wavelength_range[0], FORMAT='(F6.0)'),2), $
    ' - ', STRTRIM(STRING(wavelength_range[1], FORMAT='(F6.0)'),2), ' Angstrom'
  PRINT, '  Spectral bin width: ', STRTRIM(STRING(bin_width, FORMAT='(F4.1)'),2), ' Angstrom'
  PRINT, '  Instrumental FWHM: ', STRTRIM(STRING(fwhm, FORMAT='(F4.1)'),2), ' Angstrom'
  PRINT, '  LOS step size: ', STRTRIM(STRING(ds, FORMAT='(F5.3)'),2), ' R_J'
  PRINT, ''

  ; ============================================================================
  ; TEST CASE 1: EQUATORIAL LINE OF SIGHT (SINGLE MAXWELLIAN)
  ; ============================================================================

  PRINT, ''
  PRINT, '======================================================================'
  PRINT, 'TEST CASE 1: EQUATORIAL LINE OF SIGHT (SINGLE MAXWELLIAN)'
  PRINT, '----------------------------------------------------------------------'

  slit_pos_vec1 = [6.0d, -20.0d, 0.0d]   ; Start at x=6 R_J, equator
  norm_vec1 = [0.0d, 1.0d, 0.0d]         ; Look in +y direction

  PRINT, 'Starting position: (', STRTRIM(STRING(slit_pos_vec1[0], FORMAT='(F5.1)'),2), ', ', $
    STRTRIM(STRING(slit_pos_vec1[1], FORMAT='(F6.1)'),2), ', ', $
    STRTRIM(STRING(slit_pos_vec1[2], FORMAT='(F5.1)'),2), ') R_J'
  PRINT, 'Direction vector: (', STRTRIM(STRING(norm_vec1[0], FORMAT='(F4.1)'),2), ', ', $
    STRTRIM(STRING(norm_vec1[1], FORMAT='(F4.1)'),2), ', ', $
    STRTRIM(STRING(norm_vec1[2], FORMAT='(F4.1)'),2), ')'
  PRINT, ''
  PRINT, 'Calculating UV spectrum for single Maxwellian...'

  IPT_CALCULATE_SPECTRUM_SINGLE, plasma_model, tables_single, $
    slit_pos_vec1, norm_vec1, wave_bins1_single, spectrum1_single, lines1_single, $
    WAVELENGTH_RANGE=wavelength_range, BIN_WIDTH=bin_width, FWHM=fwhm, DS=ds

  total_brightness1_single = IPT_SIMPSON_INTEGRATE(wave_bins1_single, spectrum1_single)
  peak_idx1_single = WHERE(spectrum1_single EQ MAX(spectrum1_single))

  PRINT, ''
  PRINT, 'Single Maxwellian Results:'
  PRINT, '  Total brightness: ', STRTRIM(STRING(total_brightness1_single, FORMAT='(F10.1)'),2), ' Rayleighs'
  PRINT, '  Peak brightness: ', STRTRIM(STRING(MAX(spectrum1_single), FORMAT='(F8.2)'),2), $
    ' R/Angstrom at ', STRTRIM(STRING(wave_bins1_single[peak_idx1_single[0]], FORMAT='(F7.1)'),2), ' Angstrom'

  ; ============================================================================
  ; TEST CASE 2: OFF-EQUATOR LINE OF SIGHT (SINGLE MAXWELLIAN)
  ; ============================================================================

  PRINT, ''
  PRINT, '======================================================================'
  PRINT, 'TEST CASE 2: OFF-EQUATOR LINE OF SIGHT (SINGLE MAXWELLIAN)'
  PRINT, '----------------------------------------------------------------------'

  slit_pos_vec2 = [6.0d, -20.0d, 0.5d]   ; 0.5 R_J above equator
  norm_vec2 = [0.0d, 1.0d, 0.0d]

  PRINT, 'Starting position: (', STRTRIM(STRING(slit_pos_vec2[0], FORMAT='(F5.1)'),2), ', ', $
    STRTRIM(STRING(slit_pos_vec2[1], FORMAT='(F6.1)'),2), ', ', $
    STRTRIM(STRING(slit_pos_vec2[2], FORMAT='(F5.1)'),2), ') R_J'
  PRINT, 'Direction vector: (', STRTRIM(STRING(norm_vec2[0], FORMAT='(F4.1)'),2), ', ', $
    STRTRIM(STRING(norm_vec2[1], FORMAT='(F4.1)'),2), ', ', $
    STRTRIM(STRING(norm_vec2[2], FORMAT='(F4.1)'),2), ')'
  PRINT, ''
  PRINT, 'Calculating UV spectrum for single Maxwellian...'

  IPT_CALCULATE_SPECTRUM_SINGLE, plasma_model, tables_single, $
    slit_pos_vec2, norm_vec2, wave_bins2_single, spectrum2_single, lines2_single, $
    WAVELENGTH_RANGE=wavelength_range, BIN_WIDTH=bin_width, FWHM=fwhm, DS=ds

  total_brightness2_single = IPT_SIMPSON_INTEGRATE(wave_bins2_single, spectrum2_single)

  PRINT, ''
  PRINT, 'Single Maxwellian Results:'
  PRINT, '  Total brightness: ', STRTRIM(STRING(total_brightness2_single, FORMAT='(F10.1)'),2), ' Rayleighs'

  ; ============================================================================
  ; TEST CASE 3: EQUATORIAL LINE OF SIGHT (DOUBLE MAXWELLIAN)
  ; ============================================================================

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
  PRINT, 'Calculating UV spectrum for double Maxwellian...'

  IPT_CALCULATE_SPECTRUM_DOUBLE, plasma_model, tables_single, tables_double, $
    slit_pos_vec1, norm_vec1, wave_bins1_double, spectrum1_double, lines1_double, $
    WAVELENGTH_RANGE=wavelength_range, BIN_WIDTH=bin_width, FWHM=fwhm, DS=ds

  total_brightness1_double = IPT_SIMPSON_INTEGRATE(wave_bins1_double, spectrum1_double)
  peak_idx1_double = WHERE(spectrum1_double EQ MAX(spectrum1_double))

  PRINT, ''
  PRINT, 'Double Maxwellian Results:'
  PRINT, '  Total brightness: ', STRTRIM(STRING(total_brightness1_double, FORMAT='(F10.1)'),2), ' Rayleighs'
  PRINT, '  Peak brightness: ', STRTRIM(STRING(MAX(spectrum1_double), FORMAT='(F8.2)'),2), $
    ' R/Angstrom at ', STRTRIM(STRING(wave_bins1_double[peak_idx1_double[0]], FORMAT='(F7.1)'),2), ' Angstrom'
  IF total_brightness1_single GT 0 THEN BEGIN
    PRINT, '  Enhancement factor: ', $
      STRTRIM(STRING(total_brightness1_double/total_brightness1_single, FORMAT='(F6.2)'),2)
  ENDIF

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
  PRINT, 'Calculating UV spectrum for double Maxwellian...'

  IPT_CALCULATE_SPECTRUM_DOUBLE, plasma_model, tables_single, tables_double, $
    slit_pos_vec2, norm_vec2, wave_bins2_double, spectrum2_double, lines2_double, $
    WAVELENGTH_RANGE=wavelength_range, BIN_WIDTH=bin_width, FWHM=fwhm, DS=ds

  total_brightness2_double = IPT_SIMPSON_INTEGRATE(wave_bins2_double, spectrum2_double)

  PRINT, ''
  PRINT, 'Double Maxwellian Results:'
  PRINT, '  Total brightness: ', STRTRIM(STRING(total_brightness2_double, FORMAT='(F10.1)'),2), ' Rayleighs'
  IF total_brightness2_single GT 0 THEN BEGIN
    PRINT, '  Enhancement factor: ', $
      STRTRIM(STRING(total_brightness2_double/total_brightness2_single, FORMAT='(F6.2)'),2)
  ENDIF

  ; ============================================================================
  ; PLOT RESULTS
  ; ============================================================================

  PRINT, ''
  PRINT, '======================================================================'
  PRINT, 'CREATING PLOTS'
  PRINT, '======================================================================'

  ; Define matplotlib default colors for consistency
  mp_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

  ; ============================================================================
  ; PLOT 1: EQUATORIAL SINGLE VS DOUBLE MAXWELLIAN
  ; ============================================================================

  p1 = PLOT(wave_bins1_single, spectrum1_single, $
    COLOR=mp_colors[0], THICK=2, $
    XTITLE='Wavelength [Angstrom]', $
    YTITLE='Brightness [R/Angstrom]', $
    TITLE='IPT UV Emission: Equatorial LOS (Proper LOS Integration)', $
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

  p1.SAVE, 'ipt_uv_emission_integration_equatorial_comparison.png', RESOLUTION=300
  PRINT, '  Saved: ipt_uv_emission_integration_equatorial_comparison.png'

  ; ============================================================================
  ; PLOT 2: OFF-EQUATOR SINGLE VS DOUBLE MAXWELLIAN
  ; ============================================================================

  p2 = PLOT(wave_bins2_single, spectrum2_single, $
    COLOR=mp_colors[0], THICK=2, $
    XTITLE='Wavelength [Angstrom]', $
    YTITLE='Brightness [R/Angstrom]', $
    TITLE='IPT UV Emission: Off-Equator LOS (z = 0.5 R_J)', $
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

  p2.SAVE, 'ipt_uv_emission_integration_offequator_comparison.png', RESOLUTION=300
  PRINT, '  Saved: ipt_uv_emission_integration_offequator_comparison.png'

  ; ============================================================================
  ; PLOT 3: FOUR-PANEL COMPARISON
  ; ============================================================================

  w = WINDOW(DIMENSIONS=[1400, 1400])

  ; Panel 1: Equatorial single
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

  ; Panel 2: Equatorial double
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

  ; Panel 3: Off-equator single
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

  ; Panel 4: Off-equator double
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

  w.SAVE, 'ipt_uv_emission_integration_four_panel.png', RESOLUTION=300
  PRINT, '  Saved: ipt_uv_emission_integration_four_panel.png'

  ; ============================================================================
  ; KEY UV LINE ANALYSIS
  ; ============================================================================

  PRINT, ''
  PRINT, '======================================================================'
  PRINT, 'KEY UV DIAGNOSTIC LINES ANALYSIS'
  PRINT, '======================================================================'
  PRINT, 'Equatorial Line of Sight Results'
  PRINT, 'Instrumental FWHM: ', STRTRIM(STRING(fwhm, FORMAT='(F5.2)'),2), ' Angstrom'
  PRINT, ''

  ; Define key UV lines for IPT diagnostics
  line_names = ['S III 680.4', 'S IV 744.9', 'S IV 750.2', 'O II 833.3', $
    'O III 833.7', 'S III 1012.5', 'S IV 1062.7', 'S III 1194.0']
  line_wavelengths = [680.4d, 744.9d, 750.2d, 833.3d, 833.7d, 1012.5d, 1062.7d, 1194.0d]

  ; Integration window: +/- 3 sigma from line center
  sigma = fwhm / 2.35482d

  PRINT, STRING('Line', FORMAT='(A20)'), $
    STRING('Center [Ang]', FORMAT='(A12)'), $
    STRING('Single [R]', FORMAT='(A15)'), $
    STRING('Double [R]', FORMAT='(A15)'), $
    STRING('Enhancement', FORMAT='(A12)')
  PRINT, REPLICATE('-', 74)

  FOR i = 0, N_ELEMENTS(line_wavelengths)-1 DO BEGIN
    target_wav = line_wavelengths[i]
    wav_min = target_wav - 3.0d * sigma
    wav_max = target_wav + 3.0d * sigma

    idx = WHERE(wave_bins1_single GE wav_min AND wave_bins1_single LE wav_max, count)

    IF count GT 2 THEN BEGIN
      ; Integrate brightness over wavelength window
      int_single = IPT_SIMPSON_INTEGRATE(wave_bins1_single[idx], spectrum1_single[idx])
      int_double = IPT_SIMPSON_INTEGRATE(wave_bins1_double[idx], spectrum1_double[idx])

      ; Calculate enhancement factor
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
    ENDIF
  ENDFOR

  PRINT, REPLICATE('=', 74)

  ; ============================================================================
  ; SUMMARY TABLE
  ; ============================================================================

  PRINT, ''
  PRINT, '======================================================================'
  PRINT, 'SUMMARY TABLE'
  PRINT, '======================================================================'
  PRINT, ''
  PRINT, STRING('Test Case', FORMAT='(A30)'), $
    STRING('Total [R]', FORMAT='(A15)'), $
    STRING('Enhancement', FORMAT='(A12)')
  PRINT, REPLICATE('-', 57)
  PRINT, STRING('Equatorial Single Maxwellian', FORMAT='(A30)'), $
    STRING(total_brightness1_single, FORMAT='(F15.1)'), $
    STRING(1.00, FORMAT='(F12.2)')
  PRINT, STRING('Equatorial Double Maxwellian', FORMAT='(A30)'), $
    STRING(total_brightness1_double, FORMAT='(F15.1)'), $
    STRING(total_brightness1_double/total_brightness1_single, FORMAT='(F12.2)')
  PRINT, STRING('Off-Equator Single Maxwellian', FORMAT='(A30)'), $
    STRING(total_brightness2_single, FORMAT='(F15.1)'), $
    STRING(1.00, FORMAT='(F12.2)')
  PRINT, STRING('Off-Equator Double Maxwellian', FORMAT='(A30)'), $
    STRING(total_brightness2_double, FORMAT='(F15.1)'), $
    STRING(total_brightness2_double/total_brightness2_single, FORMAT='(F12.2)')
  PRINT, REPLICATE('=', 57)

  ; ============================================================================
  ; EXECUTION TIME
  ; ============================================================================

  end_time = SYSTIME(/SECONDS)
  elapsed_time = end_time - start_time

  PRINT, ''
  PRINT, 'Total execution time: ', STRTRIM(STRING(elapsed_time, FORMAT='(F10.2)'),2), ' seconds'
  PRINT, ''
  PRINT, '======================================================================'
  PRINT, 'Variables available for inspection:'
  PRINT, '  plasma_model          - 3D plasma model structure'
  PRINT, '  tables_single/double  - Emission table structures'
  PRINT, '  wave_bins*            - Wavelength grids [Angstrom]'
  PRINT, '  spectrum*_single      - Single Maxwellian spectra [R/Angstrom]'
  PRINT, '  spectrum*_double      - Double Maxwellian spectra [R/Angstrom]'
  PRINT, '  lines*_single/double  - Emission line lists with brightness [R]'
  PRINT, '======================================================================'

  STOP

END