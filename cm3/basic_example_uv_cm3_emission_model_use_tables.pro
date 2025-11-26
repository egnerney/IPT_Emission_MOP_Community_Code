;+
; basic_example_uv_cm3_emission_model_use_tables.pro
;
; Example Script for UV Emission Modeling of Io Plasma Torus
; ===========================================================
;
; This script demonstrates the use of the IPT_cm3_emission_model library for
; calculating UV emission spectra from the Io Plasma Torus. It shows typical
; usage for both single and double Maxwellian electron distributions and includes
; proper parameter choices for realistic IPT conditions.
;
; The example uses emission tables pre-calculated from CHIANTI 11.0.2 atomic
; database and demonstrates:
; 1. Loading emission tables for single and double Maxwellian distributions
; 2. Calculating discrete emission line brightnesses
; 3. Convolving with instrument response functions
; 4. Analyzing key diagnostic emission lines
; 5. Creating publication-quality plots
;
; REQUIRED DATA FILES:
; - Single Maxwellian: CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5
; - Double Maxwellian: CHIANTI_11.0.2_emiss_tables_double_maxwellian_24x10x24x12.h5
; - Located in ../Emiss_tables/ directory relative to this script
;
; OUTPUT:
; - Emission spectra for single and double Maxwellian cases
; - Detailed analysis of key UV emission lines
; - Publication-quality plots saved as PNG files
;
; AUTHOR: Edward (Eddie) G. Nerney
; INSTITUTION: Laboratory for Atmospheric and Space Physics, University of Colorado Boulder
; VERSION: 1.0
; DATE: November 2025
;-

; Load the IPT emission model library
@IPT_cm3_emission_model

PRO basic_example_uv_cm3_emission_model_use_tables
;+
; NAME:
;   basic_example_uv_cm3_emission_model_use_tables
;
; PURPOSE:
;   Main procedure demonstrating UV emission calculations for the Io Plasma
;   Torus using pre-calculated CHIANTI emission tables. This version does not
;   require CHIANTI to be installed - only needs the .h5 files with tables.
;
; DESCRIPTION:
;   This procedure:
;   1. Loads pre-calculated emission tables from .h5 files
;   2. Sets up typical IPT plasma parameters
;   3. Calculates UV emission spectra using table interpolation
;   4. Demonstrates both single and double Maxwellian cases
;   5. Plots the resulting spectra with matplotlib-style colors
;
; OUTPUTS:
;   Creates plots showing simulated UV spectra in Rayleighs/Angstrom
;   Saves PNG files of the spectra
;
; REQUIRED FILES:
;   For single Maxwellian:
;     CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5
;   For double Maxwellian:
;     CHIANTI_11.0.2_emiss_tables_double_maxwellian_24x10x24x12.h5
;
; MODIFICATION HISTORY:
;   Written by E. G. Nerney, November 2025
;-

  ; Record start time for performance measurement
  start_time = SYSTIME(/SECONDS)
  
  PRINT, '=================================================================='
  PRINT, 'UV EMISSION MODEL FOR IO PLASMA TORUS'
  PRINT, 'Table-Based Implementation using CHIANTI 11.0.2'
  PRINT, '=================================================================='
  PRINT, ''
  
  ; ============================================================================
  ; WAVELENGTH GRID SETUP
  ; ============================================================================
  ; Define the wavelength range and resolution for the UV simulation
  ; Europa-UVS and JUICE-UVS wavelength range
  
  min_xwav = 550.0d        ; Minimum wavelength [Angstroms]
  max_xwav = 2100.0d       ; Maximum wavelength [Angstroms]
  xwav_bin_width = 1.0d    ; Spectral bin width [Angstroms]
  
  ; Calculate number of wavelength bins
  num_xwav_points = FLOOR((max_xwav - min_xwav)/xwav_bin_width) + 1
  
  ; Create wavelength grid (bin centers)
  xwav = xwav_bin_width * DINDGEN(num_xwav_points) + min_xwav
  
  ; ============================================================================
  ; INSTRUMENT PARAMETERS
  ; ============================================================================
  ; Full Width Half Maximum of instrument response function
  ; Best case scenario for IPT observations from JUICE/Europa-UVS
  fwhm = 6.0d  ; [Angstroms]
  
  ; ============================================================================
  ; PLASMA PARAMETERS
  ; ============================================================================
  ; Core/cold electron population parameters
  ; Typical IPT values from UV line ratio diagnostics
  Tec = 5.0d        ; Core electron temperature [eV]
  nec = 2200.0d     ; Core electron density [cm^-3] - peak torus density
  
  ; Hot electron component (suprathermal population)
  ; From Voyager in situ measurements and wave-particle heating models
  Teh = 270.0d      ; Hot electron temperature [eV]
  feh = 0.0025d     ; Fraction of hot electrons (0.25% typical)
  fec = 1.0d - feh  ; Fraction of cold electrons
  
  ; Calculate densities for double Maxwellian
  neh = nec * (1.0d/fec - 1.0d)  ; Hot electron density [cm^-3]
  ne_total = nec + neh           ; Total electron density [cm^-3]
  
  ; ============================================================================
  ; COLUMN DENSITIES
  ; ============================================================================
  ; Ion column densities [cm^-2]
  ; Based on Nerney et al. (2017), Steffl et al. (2004b), Thomas et al. (2004)
  ; These values correspond to ~10 R_J path length through the Io torus
  
  col_Sp = 1.2d13    ; S+ (S II) - singly ionized sulfur
  col_S2p = 4.2d13   ; S++ (S III) - dominant sulfur ion (~90% of S)
  col_S3p = 5.92d12  ; S+++ (S IV) - minor component
  col_S4p = 6.0d11   ; S++++ (S V) - highly ionized, trace amounts
  col_Op = 5.2d13    ; O+ (O II) - dominant ion overall (~80% of O)
  col_O2p = 5.92d12  ; O++ (O III) - secondary oxygen ion
  
  ; ============================================================================
  ; MATPLOTLIB DEFAULT COLORS
  ; ============================================================================
  ; Define colors matching matplotlib's default color cycle for consistency
  ; These are RGB triplets [R, G, B] with values 0-255
  mp_color_C0 = [31, 119, 180]    ; Blue
  mp_color_C1 = [255, 127, 14]    ; Orange
  mp_color_C2 = [44, 160, 44]     ; Green
  mp_color_C3 = [214, 39, 40]     ; Red
  mp_color_C4 = [148, 103, 189]   ; Purple
  
  ; ============================================================================
  ; LOAD SINGLE MAXWELLIAN TABLES
  ; ============================================================================
  PRINT, '=================================================================='
  PRINT, 'SINGLE MAXWELLIAN CALCULATION'
  PRINT, '=================================================================='
  
  ; Set file path for single Maxwellian HDF5 file
  emission_filename_single = 'CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5'
  emission_file_single = FILEPATH(emission_filename_single, $
                                  ROOT_DIR='..', SUBDIRECTORY='Emiss_tables')
  
  ; Load the tables
  tables_single = IPT_LOAD_SINGLE_MAXWELLIAN_TABLES(emission_file_single)
  
  ; Calculate emission line brightnesses
  yptsi_single = IPT_CALCULATE_EMISSION_SINGLE(tables_single, Tec, nec, $
    col_Sp, col_S2p, col_S3p, col_S4p, col_Op, col_O2p, $
    MIN_WAV=min_xwav, MAX_WAV=max_xwav, WAVELENGTHS=xwavi_single)
  
  ; Convolve with instrument response
  PRINT, 'Convolving with instrument response function...'
  ypts_single_maxwellian = IPT_SIMULATE_SPECTRUM_ERF(xwav, xwav_bin_width, $
    xwavi_single, yptsi_single, FWHM=fwhm)
  
  ; ============================================================================
  ; LOAD DOUBLE MAXWELLIAN TABLES
  ; ============================================================================
  PRINT, ''
  PRINT, '=================================================================='
  PRINT, 'DOUBLE MAXWELLIAN CALCULATION'
  PRINT, '=================================================================='
  
  ; Set file path for double Maxwellian HDF5 file
  emission_filename_double = 'CHIANTI_11.0.2_emiss_tables_double_maxwellian_24x10x24x12.h5'
  emission_file_double = FILEPATH(emission_filename_double, $
                                  ROOT_DIR='..', SUBDIRECTORY='Emiss_tables')
  
  ; Load the tables
  tables_double = IPT_LOAD_DOUBLE_MAXWELLIAN_TABLES(emission_file_double)
  
  ; Calculate emission line brightnesses
  yptsi_double = IPT_CALCULATE_EMISSION_DOUBLE(tables_double, Tec, Teh, $
    ne_total, feh, col_Sp, col_S2p, col_S3p, col_S4p, col_Op, col_O2p, $
    MIN_WAV=min_xwav, MAX_WAV=max_xwav, WAVELENGTHS=xwavi_double)
  
  ; Convolve with instrument response
  PRINT, 'Convolving with instrument response function...'
  ypts_double_maxwellian = IPT_SIMULATE_SPECTRUM_ERF(xwav, xwav_bin_width, $
    xwavi_double, yptsi_double, FWHM=fwhm)
  
  ; ============================================================================
  ; SUMMARY OUTPUT
  ; ============================================================================
  PRINT, ''
  PRINT, '=================================================================='
  PRINT, 'SIMULATION COMPLETE'
  PRINT, '=================================================================='
  PRINT, 'Plasma parameters:'
  PRINT, '  Single Maxwellian: Te = ', STRTRIM(STRING(Tec, FORMAT='(F5.1)'),2), $
         ' eV, ne = ', STRTRIM(STRING(nec, FORMAT='(F7.1)'),2), ' cm^-3'
  PRINT, '  Double Maxwellian: Tec = ', STRTRIM(STRING(Tec, FORMAT='(F5.1)'),2), $
         ' eV, Teh = ', STRTRIM(STRING(Teh, FORMAT='(F6.1)'),2), ' eV'
  PRINT, '                     feh = ', STRTRIM(STRING(feh, FORMAT='(F7.4)'),2), $
         ', ne_total = ', STRTRIM(STRING(ne_total, FORMAT='(F8.1)'),2), ' cm^-3'
  PRINT, ''
  PRINT, 'Column densities [cm^-2]:'
  PRINT, '  S+    : ', STRTRIM(STRING(col_Sp, FORMAT='(E10.2)'),2)
  PRINT, '  S++   : ', STRTRIM(STRING(col_S2p, FORMAT='(E10.2)'),2)
  PRINT, '  S+++  : ', STRTRIM(STRING(col_S3p, FORMAT='(E10.2)'),2)
  PRINT, '  S++++ : ', STRTRIM(STRING(col_S4p, FORMAT='(E10.2)'),2)
  PRINT, '  O+    : ', STRTRIM(STRING(col_Op, FORMAT='(E10.2)'),2)
  PRINT, '  O++   : ', STRTRIM(STRING(col_O2p, FORMAT='(E10.2)'),2)
  PRINT, ''
  PRINT, 'Number of emission lines in wavelength range:'
  PRINT, '  Single Maxwellian: ', STRTRIM(STRING(N_ELEMENTS(xwavi_single)),2)
  PRINT, '  Double Maxwellian: ', STRTRIM(STRING(N_ELEMENTS(xwavi_double)),2)
  PRINT, ''
  PRINT, 'Output spectrum statistics:'
  PRINT, '  Single Maxwellian peak: ', STRTRIM(STRING(MAX(ypts_single_maxwellian), $
         FORMAT='(F10.2)'),2), ' R/Angstrom'
  PRINT, '  Double Maxwellian peak: ', STRTRIM(STRING(MAX(ypts_double_maxwellian), $
         FORMAT='(F10.2)'),2), ' R/Angstrom'
  
  ; Calculate integrated brightness using Simpson's rule
  int_single = IPT_SIMPSON_INTEGRATE(xwav, ypts_single_maxwellian)
  int_double = IPT_SIMPSON_INTEGRATE(xwav, ypts_double_maxwellian)
  PRINT, '  Integrated brightness (single): ', STRTRIM(STRING(int_single, $
         FORMAT='(F12.0)'),2), ' R'
  PRINT, '  Integrated brightness (double): ', STRTRIM(STRING(int_double, $
         FORMAT='(F12.0)'),2), ' R'
  
  ; ============================================================================
  ; PLOTTING
  ; ============================================================================
  PRINT, ''
  PRINT, 'Creating plots...'
  
  ; Plot 1: Single Maxwellian spectrum
  p1 = PLOT(xwav, ypts_single_maxwellian, $
    COLOR=mp_color_C0, $
    THICK=2, $
    XTITLE='Wavelength [Angstrom]', $
    YTITLE='Brightness [R/Angstrom]', $
    TITLE='IPT UV Emission: Single Maxwellian (Te=' + $
          STRTRIM(STRING(Tec, FORMAT='(F4.1)'),2) + ' eV, ne=' + $
          STRTRIM(STRING(nec, FORMAT='(F6.0)'),2) + ' cm$^{-3}$)', $
    XRANGE=[min_xwav, max_xwav], $
    YRANGE=[0, MAX(ypts_single_maxwellian)*1.1], $
    FONT_SIZE=12, $
    DIMENSIONS=[1200, 500])
  
  ; Add grid
  p1.XGRIDSTYLE = 1
  p1.YGRIDSTYLE = 1
  p1.XMINOR = 0
  p1.YMINOR = 0
  
  p1.SAVE, 'ipt_uv_emission_single_maxwellian.png', RESOLUTION=300
  PRINT, '  Saved: ipt_uv_emission_single_maxwellian.png'
  
  ; Plot 2: Double Maxwellian spectrum
  p2 = PLOT(xwav, ypts_double_maxwellian, $
    COLOR=mp_color_C1, $
    THICK=2, $
    XTITLE='Wavelength [Angstrom]', $
    YTITLE='Brightness [R/Angstrom]', $
    TITLE='IPT UV Emission: Double Maxwellian (Tec=' + $
          STRTRIM(STRING(Tec, FORMAT='(F4.1)'),2) + ' eV, Teh=' + $
          STRTRIM(STRING(Teh, FORMAT='(F5.0)'),2) + ' eV, feh=' + $
          STRTRIM(STRING(feh, FORMAT='(F6.4)'),2) + ')', $
    XRANGE=[min_xwav, max_xwav], $
    YRANGE=[0, MAX(ypts_double_maxwellian)*1.1], $
    FONT_SIZE=12, $
    DIMENSIONS=[1200, 500])
  
  p2.XGRIDSTYLE = 1
  p2.YGRIDSTYLE = 1
  
  p2.SAVE, 'ipt_uv_emission_double_maxwellian.png', RESOLUTION=300
  PRINT, '  Saved: ipt_uv_emission_double_maxwellian.png'
  
  ; Plot 3: Comparison plot with both spectra
  p3 = PLOT(xwav, ypts_single_maxwellian, $
    COLOR=mp_color_C0, $
    THICK=2, $
    XTITLE='Wavelength [Angstrom]', $
    YTITLE='Brightness [R/Angstrom]', $
    TITLE='IPT UV Emission Comparison: Single vs Double Maxwellian', $
    XRANGE=[min_xwav, max_xwav], $
    YRANGE=[0, MAX([MAX(ypts_single_maxwellian), MAX(ypts_double_maxwellian)])*1.1], $
    FONT_SIZE=12, $
    DIMENSIONS=[1200, 600], $
    NAME='Single Maxwellian')
  
  p3_double = PLOT(xwav, ypts_double_maxwellian, $
    /OVERPLOT, $
    COLOR=mp_color_C1, $
    THICK=2, $
    TRANSPARENCY=30, $
    NAME='Double Maxwellian')
  
  ; Add legend
  leg = LEGEND(TARGET=[p3, p3_double], $
    /RELATIVE, $
    POSITION=[0.85, 0.85], $
    HORIZONTAL_ALIGNMENT='right', $
    VERTICAL_ALIGNMENT='top', $
    /AUTO_TEXT_COLOR, $
    FONT_SIZE=11)
  
  p3.XGRIDSTYLE = 1
  p3.YGRIDSTYLE = 1
  
  p3.SAVE, 'ipt_uv_emission_comparison.png', RESOLUTION=300
  PRINT, '  Saved: ipt_uv_emission_comparison.png'
  
  ; Plot 4: Two-panel comparison (like Python version)
  ; Create window with two plots
  w = WINDOW(DIMENSIONS=[1200, 1000])
  
  ; Top panel: Single Maxwellian
  p4a = PLOT(xwav, ypts_single_maxwellian, $
    COLOR=mp_color_C0, $
    THICK=2, $
    XTITLE='Wavelength [Angstrom]', $
    YTITLE='Brightness [R/Angstrom]', $
    TITLE='IPT UV Emission: Single Maxwellian (Te=' + $
          STRTRIM(STRING(Tec, FORMAT='(F4.1)'),2) + ' eV, ne=' + $
          STRTRIM(STRING(nec, FORMAT='(F6.0)'),2) + ' cm$^{-3}$)', $
    XRANGE=[min_xwav, max_xwav], $
    YRANGE=[0, MAX(ypts_single_maxwellian)*1.1], $
    FONT_SIZE=11, $
    LAYOUT=[1,2,1], $
    /CURRENT)
  
  p4a.XGRIDSTYLE = 1
  p4a.YGRIDSTYLE = 1
  
  ; Bottom panel: Double Maxwellian
  p4b = PLOT(xwav, ypts_double_maxwellian, $
    COLOR=mp_color_C1, $
    THICK=2, $
    XTITLE='Wavelength [Angstrom]', $
    YTITLE='Brightness [R/Angstrom]', $
    TITLE='IPT UV Emission: Double Maxwellian (Tec=' + $
          STRTRIM(STRING(Tec, FORMAT='(F4.1)'),2) + ' eV, Teh=' + $
          STRTRIM(STRING(Teh, FORMAT='(F5.0)'),2) + ' eV, feh=' + $
          STRTRIM(STRING(feh, FORMAT='(F6.4)'),2) + ')', $
    XRANGE=[min_xwav, max_xwav], $
    YRANGE=[0, MAX(ypts_double_maxwellian)*1.1], $
    FONT_SIZE=11, $
    LAYOUT=[1,2,2], $
    /CURRENT)
  
  p4b.XGRIDSTYLE = 1
  p4b.YGRIDSTYLE = 1
  
  w.SAVE, 'ipt_uv_emission_single_vs_double.png', RESOLUTION=300
  PRINT, '  Saved: ipt_uv_emission_single_vs_double.png'
  
  ; ============================================================================
  ; ANALYZE KEY UV EMISSION LINES
  ; ============================================================================
  PRINT, ''
  PRINT, '=================================================================='
  PRINT, 'DETAILED UV EMISSION LINE ANALYSIS'
  PRINT, '=================================================================='
  PRINT, 'Analyzing key diagnostic lines for IPT UV spectroscopy'
  PRINT, 'Instrumental FWHM: ', STRTRIM(STRING(fwhm, FORMAT='(F5.2)'),2), ' Angstrom'
  
  ; Define key UV emission lines for analysis
  ; These are important diagnostics for IPT plasma parameters
  line_names = ['S III 680.4', 'S IV 744.9', 'S IV 750.2', 'O II 833.3', $
                'O III 833.7', 'S III 1012.5', 'S IV 1062.7', 'S III 1194.0']
  line_wavelengths = [680.4d, 744.9d, 750.2d, 833.3d, 833.7d, 1012.5d, 1062.7d, 1194.0d]
  line_species = ['S++', 'S+++', 'S+++', 'O+', 'O++', 'S++', 'S+++', 'S++']
  
  ; Analyze each line
  sigma = fwhm / 2.35482d  ; Convert FWHM to sigma
  
  PRINT, ''
  PRINT, STRING('Line', FORMAT='(A20)'), $
         STRING('Peak [Ang]', FORMAT='(A12)'), $
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
    idx = WHERE(xwav GE wav_min AND xwav LE wav_max, count)
    
    IF count GT 2 THEN BEGIN
      ; Integrate using Simpson's rule
      int_single_line = IPT_SIMPSON_INTEGRATE(xwav[idx], ypts_single_maxwellian[idx])
      int_double_line = IPT_SIMPSON_INTEGRATE(xwav[idx], ypts_double_maxwellian[idx])
      
      ; Calculate enhancement factor
      enhancement = int_double_line / int_single_line
      
      ; Find peak wavelength
      peak_idx = WHERE(ypts_double_maxwellian[idx] EQ MAX(ypts_double_maxwellian[idx]))
      peak_wav = xwav[idx[peak_idx[0]]]
      
      PRINT, STRING(line_names[i] + ' Ang', FORMAT='(A20)'), $
             STRING(peak_wav, FORMAT='(F12.3)'), $
             STRING(int_single_line, FORMAT='(F15.1)'), $
             STRING(int_double_line, FORMAT='(F15.1)'), $
             STRING(enhancement, FORMAT='(F12.3)')
    ENDIF
  ENDFOR
  
  PRINT, REPLICATE('=', 74)
  
  ; ============================================================================
  ; EXECUTION TIME
  ; ============================================================================
  end_time = SYSTIME(/SECONDS)
  elapsed_time = end_time - start_time
  
  PRINT, ''
  PRINT, 'Execution time: ', STRTRIM(STRING(elapsed_time, FORMAT='(F8.4)'),2), ' seconds'
  PRINT, ''
  PRINT, '=================================================================='
  PRINT, 'Variables available for inspection:'
  PRINT, '  xwav                    - wavelength grid [Angstrom]'
  PRINT, '  ypts_single_maxwellian  - single Maxwellian spectrum [R/Ang]'
  PRINT, '  ypts_double_maxwellian  - double Maxwellian spectrum [R/Ang]'
  PRINT, '  xwavi_single/double     - emission line wavelengths [Ang]'
  PRINT, '  yptsi_single/double     - emission line brightnesses [R]'
  PRINT, '  tables_single/double    - loaded emission table structures'
  PRINT, '=================================================================='
  
  ; Stop to allow interactive examination of results
  STOP
  
END