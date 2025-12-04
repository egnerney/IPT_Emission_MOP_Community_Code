;+
; basic_example_optical_cm3_emission_model_use_tables.pro
;
; Example Script for Optical Emission Modeling of Io Plasma Torus
; ================================================================
;
; This script demonstrates the use of the IPT_cm3_emission_model library for
; calculating optical emission spectra from the Io Plasma Torus. It shows typical
; usage for both single and double Maxwellian electron distributions with parameters
; appropriate for ground-based optical spectroscopy.
;
; The example uses emission tables pre-calculated from CHIANTI 11.0.2 atomic
; database and demonstrates:
; 1. Loading emission tables for single and double Maxwellian distributions
; 2. Calculating discrete emission line brightnesses in optical wavelengths
; 3. Convolving with ground-based telescope instrument response functions
; 4. Analyzing key diagnostic optical emission lines
; 5. Creating publication-quality plots for optical spectroscopy
;
; REQUIRED DATA FILES:
; - Single Maxwellian: CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5
; - Double Maxwellian: CHIANTI_11.0.2_emiss_tables_double_maxwellian_16x8x12x10.h5
; - Located in ../../Emiss_tables/ directory relative to this script
;
; OUTPUT:
; - Optical emission spectra for single and double Maxwellian cases
; - Detailed analysis of key optical emission lines ([S II], [O II], [O III], [S III])
; - Publication-quality plots saved as PNG files
;
; WAVELENGTH RANGE:
; - Optical: 3000-10000 Angstroms (suitable for ground-based telescopes)
;
; AUTHOR: Edward (Eddie) G. Nerney
; INSTITUTION: Laboratory for Atmospheric and Space Physics, University of Colorado Boulder
; VERSION: 1.0
; DATE: November 2025
;-

; Load the IPT emission model library
@IPT_cm3_emission_model

PRO basic_example_optical_cm3_emission_model_use_tables
;+
; NAME:
;   basic_example_optical_cm3_emission_model_use_tables
;
; PURPOSE:
;   Main procedure demonstrating optical emission calculations for the Io Plasma
;   Torus using pre-calculated CHIANTI emission tables. This version is designed
;   for ground-based telescope observations and analysis.
;
; DESCRIPTION:
;   This procedure:
;   1. Loads pre-calculated emission tables from .h5 files
;   2. Sets up typical IPT plasma parameters
;   3. Calculates optical emission spectra (3000-10000 Angstrom)
;   4. Demonstrates both single and double Maxwellian cases
;   5. Analyzes key optical diagnostic lines ([S II], [O II], [O III])
;   6. Plots the resulting spectra with matplotlib-style colors
;
; OUTPUTS:
;   Creates plots showing simulated optical spectra in Rayleighs/Angstrom
;   Saves PNG files of the spectra
;   Prints summary of key optical line diagnostics
;
; REQUIRED FILES:
;   For single Maxwellian:
;     CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5
;   For double Maxwellian:
;     CHIANTI_11.0.2_emiss_tables_double_maxwellian_16x8x12x10.h5
;
; OPTICAL LINE DIAGNOSTICS:
;   Key lines analyzed include:
;   - [O II] 3726/3729 Angstrom doublet (density diagnostic)
;   - [S II] 4069/4076 Angstrom (auroral lines)
;   - [O III] 4363 Angstrom (temperature diagnostic)
;   - [O III] 4959/5007 Angstrom (nebular lines)
;   - [S III] 6312 Angstrom (temperature diagnostic)
;   - [S II] 6716/6731 Angstrom (density diagnostic doublet)
;   - [S III] 9069 Angstrom (near-IR forbidden line)
;
; MODIFICATION HISTORY:
;   Written by E. G. Nerney, November 2025
;-

  ; Record start time for performance measurement
  start_time = SYSTIME(/SECONDS)
  
  PRINT, '=================================================================='
  PRINT, 'OPTICAL EMISSION MODEL FOR IO PLASMA TORUS'
  PRINT, 'Ground-Based Telescope Spectroscopy'
  PRINT, 'Table-Based Implementation using CHIANTI 11.0.2'
  PRINT, '=================================================================='
  PRINT, ''
  
  ; ============================================================================
  ; WAVELENGTH GRID SETUP - OPTICAL RANGE
  ; ============================================================================
  ; Define the wavelength range and resolution for optical simulation
  ; Ground-based optical telescope wavelength range
  
  min_xwav = 3000.0d       ; Minimum wavelength [Angstroms] - blue optical
  max_xwav = 10000.0d      ; Maximum wavelength [Angstroms] - near-infrared limit
  xwav_bin_width = 2.0d    ; Spectral bin width [Angstroms] - typical for R~3000
  
  ; Calculate number of wavelength bins
  num_xwav_points = FLOOR((max_xwav - min_xwav)/xwav_bin_width) + 1
  
  ; Create wavelength grid (bin centers)
  xwav = xwav_bin_width * DINDGEN(num_xwav_points) + min_xwav
  
  ; ============================================================================
  ; INSTRUMENT PARAMETERS
  ; ============================================================================
  ; Full Width Half Maximum of instrument response function
  ; Typical for ground-based spectrograph with R~3000
  fwhm = 3.0d  ; [Angstroms]
  
  ; ============================================================================
  ; PLASMA PARAMETERS
  ; ============================================================================
  ; Core/cold electron population parameters
  ; Typical IPT values from UV/optical line ratio diagnostics
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
  PRINT, 'SINGLE MAXWELLIAN CALCULATION - OPTICAL WAVELENGTHS'
  PRINT, '=================================================================='
  
  ; Set file path for single Maxwellian HDF5 file
  emission_filename_single = 'CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5'
  emission_file_single = FILEPATH(emission_filename_single, $
                                  ROOT_DIR='..', SUBDIRECTORY=['..', 'Emiss_tables'])
  
  ; Load the tables
  tables_single = IPT_LOAD_SINGLE_MAXWELLIAN_TABLES(emission_file_single)
  
  ; Calculate emission line brightnesses in optical range
  yptsi_single = IPT_CALCULATE_EMISSION_SINGLE(tables_single, Tec, nec, $
    col_Sp, col_S2p, col_S3p, col_S4p, col_Op, col_O2p, $
    MIN_WAV=min_xwav, MAX_WAV=max_xwav, WAVELENGTHS=xwavi_single)
  
  ; Convolve with instrument response
  PRINT, 'Convolving with ground-based telescope instrument response...'
  ypts_single_maxwellian = IPT_SIMULATE_SPECTRUM_ERF(xwav, xwav_bin_width, $
    xwavi_single, yptsi_single, FWHM=fwhm)
  
  ; ============================================================================
  ; LOAD DOUBLE MAXWELLIAN TABLES
  ; ============================================================================
  PRINT, ''
  PRINT, '=================================================================='
  PRINT, 'DOUBLE MAXWELLIAN CALCULATION - OPTICAL WAVELENGTHS'
  PRINT, '=================================================================='
  
  ; Set file path for double Maxwellian HDF5 file
  emission_filename_double = 'CHIANTI_11.0.2_emiss_tables_double_maxwellian_16x8x12x10.h5'
  emission_file_double = FILEPATH(emission_filename_double, $
                                  ROOT_DIR='..', SUBDIRECTORY=['..', 'Emiss_tables'])
  
  ; Load the tables
  tables_double = IPT_LOAD_DOUBLE_MAXWELLIAN_TABLES(emission_file_double)
  
  ; Calculate emission line brightnesses in optical range
  yptsi_double = IPT_CALCULATE_EMISSION_DOUBLE(tables_double, Tec, Teh, $
    ne_total, feh, col_Sp, col_S2p, col_S3p, col_S4p, col_Op, col_O2p, $
    MIN_WAV=min_xwav, MAX_WAV=max_xwav, WAVELENGTHS=xwavi_double)
  
  ; Convolve with instrument response
  PRINT, 'Convolving with ground-based telescope instrument response...'
  ypts_double_maxwellian = IPT_SIMULATE_SPECTRUM_ERF(xwav, xwav_bin_width, $
    xwavi_double, yptsi_double, FWHM=fwhm)
  
  ; ============================================================================
  ; SUMMARY OUTPUT
  ; ============================================================================
  PRINT, ''
  PRINT, '=================================================================='
  PRINT, 'OPTICAL SIMULATION COMPLETE'
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
  PRINT, 'Optical wavelength range:'
  PRINT, '  ', STRTRIM(STRING(min_xwav, FORMAT='(F8.0)'),2), ' - ', $
         STRTRIM(STRING(max_xwav, FORMAT='(F8.0)'),2), ' Angstrom'
  PRINT, '  Spectral resolution (FWHM): ', STRTRIM(STRING(fwhm, FORMAT='(F5.1)'),2), $
         ' Angstrom'
  PRINT, '  Bin width: ', STRTRIM(STRING(xwav_bin_width, FORMAT='(F5.1)'),2), ' Angstrom'
  PRINT, ''
  PRINT, 'Number of emission lines in optical wavelength range:'
  PRINT, '  Single Maxwellian: ', STRTRIM(STRING(N_ELEMENTS(xwavi_single)),2)
  PRINT, '  Double Maxwellian: ', STRTRIM(STRING(N_ELEMENTS(xwavi_double)),2)
  PRINT, ''
  PRINT, 'Output optical spectrum statistics:'
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
  ; PLOTTING - OPTICAL SPECTRA
  ; ============================================================================
  PRINT, ''
  PRINT, 'Creating optical spectrum plots...'
  
  ; Plot 1: Single Maxwellian optical spectrum
  p1 = PLOT(xwav, ypts_single_maxwellian, $
    COLOR=mp_color_C0, $
    THICK=2, $
    XTITLE='Wavelength [Angstrom]', $
    YTITLE='Brightness [R/Angstrom]', $
    TITLE='IPT Optical Emission: Single Maxwellian (Te=' + $
          STRTRIM(STRING(Tec, FORMAT='(F4.1)'),2) + ' eV, ne=' + $
          STRTRIM(STRING(nec, FORMAT='(F6.0)'),2) + ' cm$^{-3}$)', $
    XRANGE=[min_xwav, max_xwav], $
    YRANGE=[0, MAX(ypts_single_maxwellian)*1.1], $
    FONT_SIZE=12, $
    DIMENSIONS=[1200, 500])
  
  p1.XGRIDSTYLE = 1
  p1.YGRIDSTYLE = 1
  
  p1.SAVE, 'ipt_optical_emission_single_maxwellian.png', RESOLUTION=300
  PRINT, '  Saved: ipt_optical_emission_single_maxwellian.png'
  
  ; Plot 2: Double Maxwellian optical spectrum
  p2 = PLOT(xwav, ypts_double_maxwellian, $
    COLOR=mp_color_C1, $
    THICK=2, $
    XTITLE='Wavelength [Angstrom]', $
    YTITLE='Brightness [R/Angstrom]', $
    TITLE='IPT Optical Emission: Double Maxwellian (Tec=' + $
          STRTRIM(STRING(Tec, FORMAT='(F4.1)'),2) + ' eV, Teh=' + $
          STRTRIM(STRING(Teh, FORMAT='(F5.0)'),2) + ' eV, feh=' + $
          STRTRIM(STRING(feh, FORMAT='(F6.4)'),2) + ')', $
    XRANGE=[min_xwav, max_xwav], $
    YRANGE=[0, MAX(ypts_double_maxwellian)*1.1], $
    FONT_SIZE=12, $
    DIMENSIONS=[1200, 500])
  
  p2.XGRIDSTYLE = 1
  p2.YGRIDSTYLE = 1
  
  p2.SAVE, 'ipt_optical_emission_double_maxwellian.png', RESOLUTION=300
  PRINT, '  Saved: ipt_optical_emission_double_maxwellian.png'
  
  ; Plot 3: Comparison plot with both spectra
  p3 = PLOT(xwav, ypts_single_maxwellian, $
    COLOR=mp_color_C0, $
    THICK=2, $
    XTITLE='Wavelength [Angstrom]', $
    YTITLE='Brightness [R/Angstrom]', $
    TITLE='IPT Optical Emission Comparison: Single vs Double Maxwellian', $
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
  
  p3.SAVE, 'ipt_optical_emission_comparison.png', RESOLUTION=300
  PRINT, '  Saved: ipt_optical_emission_comparison.png'
  
  ; Plot 4: Two-panel comparison
  w = WINDOW(DIMENSIONS=[1200, 1000])
  
  ; Top panel: Single Maxwellian
  p4a = PLOT(xwav, ypts_single_maxwellian, $
    COLOR=mp_color_C0, $
    THICK=2, $
    XTITLE='Wavelength [Angstrom]', $
    YTITLE='Brightness [R/Angstrom]', $
    TITLE='IPT Optical Emission: Single Maxwellian (Te=' + $
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
    TITLE='IPT Optical Emission: Double Maxwellian (Tec=' + $
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
  
  w.SAVE, 'ipt_optical_emission_single_vs_double.png', RESOLUTION=300
  PRINT, '  Saved: ipt_optical_emission_single_vs_double.png'
  
  ; Plot 5: Zoom on [S II] 6716/6731 doublet region
  ; This doublet is a key density diagnostic for ground-based observations
  mask_sii = WHERE(xwav GE 6700 AND xwav LE 6750)
  
  p5 = PLOT(xwav[mask_sii], ypts_single_maxwellian[mask_sii], $
    COLOR=mp_color_C0, $
    THICK=2, $
    XTITLE='Wavelength [Angstrom]', $
    YTITLE='Brightness [R/Angstrom]', $
    TITLE='[S II] 6716/6731 Angstrom Doublet Region (Density Diagnostic)', $
    XRANGE=[6700, 6750], $
    FONT_SIZE=12, $
    DIMENSIONS=[1000, 600], $
    NAME='Single Maxwellian')
  
  p5_double = PLOT(xwav[mask_sii], ypts_double_maxwellian[mask_sii], $
    /OVERPLOT, $
    COLOR=mp_color_C1, $
    THICK=2, $
    TRANSPARENCY=30, $
    NAME='Double Maxwellian')
  
  leg5 = LEGEND(TARGET=[p5, p5_double], $
    /RELATIVE, $
    POSITION=[0.85, 0.85], $
    HORIZONTAL_ALIGNMENT='right', $
    VERTICAL_ALIGNMENT='top', $
    /AUTO_TEXT_COLOR, $
    FONT_SIZE=11)
  
  p5.XGRIDSTYLE = 1
  p5.YGRIDSTYLE = 1
  
  p5.SAVE, 'ipt_optical_sii_doublet_region.png', RESOLUTION=300
  PRINT, '  Saved: ipt_optical_sii_doublet_region.png'
  
  ; Plot 6: Zoom on [S III] 6312 and surrounding region
  mask_siii = WHERE(xwav GE 6280 AND xwav LE 6330)
  
  p6 = PLOT(xwav[mask_siii], ypts_single_maxwellian[mask_siii], $
    COLOR=mp_color_C0, $
    THICK=2, $
    XTITLE='Wavelength [Angstrom]', $
    YTITLE='Brightness [R/Angstrom]', $
    TITLE='[S III] 6312 Angstrom Region (Temperature Diagnostic)', $
    XRANGE=[6280, 6330], $
    FONT_SIZE=12, $
    DIMENSIONS=[1000, 600], $
    NAME='Single Maxwellian')
  
  p6_double = PLOT(xwav[mask_siii], ypts_double_maxwellian[mask_siii], $
    /OVERPLOT, $
    COLOR=mp_color_C1, $
    THICK=2, $
    TRANSPARENCY=30, $
    NAME='Double Maxwellian')
  
  leg6 = LEGEND(TARGET=[p6, p6_double], $
    /RELATIVE, $
    POSITION=[0.85, 0.85], $
    HORIZONTAL_ALIGNMENT='right', $
    VERTICAL_ALIGNMENT='top', $
    /AUTO_TEXT_COLOR, $
    FONT_SIZE=11)
  
  p6.XGRIDSTYLE = 1
  p6.YGRIDSTYLE = 1
  
  p6.SAVE, 'ipt_optical_siii_6312_region.png', RESOLUTION=300
  PRINT, '  Saved: ipt_optical_siii_6312_region.png'
  
  ; ============================================================================
  ; ANALYZE KEY OPTICAL EMISSION LINES
  ; ============================================================================
  PRINT, ''
  PRINT, '=================================================================='
  PRINT, 'DETAILED OPTICAL EMISSION LINE ANALYSIS'
  PRINT, '=================================================================='
  PRINT, 'Analyzing key diagnostic lines for IPT optical spectroscopy'
  PRINT, 'Instrumental FWHM: ', STRTRIM(STRING(fwhm, FORMAT='(F5.2)'),2), ' Angstrom'
  
  ; Define key optical emission lines for analysis
  ; These are important diagnostics from ground-based observations
  line_names = ['[O II] 3726.0', '[O II] 3728.8', '[S II] 4068.6', '[S II] 4076.3', $
                '[O III] 4363.2', '[O III] 4958.9', '[O III] 5006.8', '[S III] 6312.1', $
                '[S II] 6716.4', '[S II] 6730.8', '[S III] 9068.6']
  line_wavelengths = [3726.03d, 3728.82d, 4068.6d, 4076.3d, 4363.21d, 4958.91d, $
                      5006.84d, 6312.1d, 6716.4d, 6730.8d, 9068.6d]
  line_species = ['O+', 'O+', 'S+', 'S+', 'O++', 'O++', 'O++', 'S++', 'S+', 'S+', 'S++']
  
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
      IF int_single_line GT 0 THEN BEGIN
        enhancement = int_double_line / int_single_line
      ENDIF ELSE BEGIN
        enhancement = 0.0d
      ENDELSE
      
      ; Find peak wavelength
      peak_idx = WHERE(ypts_double_maxwellian[idx] EQ MAX(ypts_double_maxwellian[idx]))
      peak_wav = xwav[idx[peak_idx[0]]]
      
      PRINT, STRING(line_names[i], FORMAT='(A20)'), $
             STRING(peak_wav, FORMAT='(F12.3)'), $
             STRING(int_single_line, FORMAT='(F15.1)'), $
             STRING(int_double_line, FORMAT='(F15.1)'), $
             STRING(enhancement, FORMAT='(F12.3)')
    ENDIF ELSE BEGIN
      PRINT, STRING(line_names[i], FORMAT='(A20)'), $
             '  (outside wavelength range or insufficient data)'
    ENDELSE
  ENDFOR
  
  PRINT, REPLICATE('=', 74)
  
  ; ============================================================================
  ; DENSITY DIAGNOSTIC RATIOS
  ; ============================================================================
  PRINT, ''
  PRINT, '=================================================================='
  PRINT, 'OPTICAL LINE RATIO DIAGNOSTICS'
  PRINT, '=================================================================='
  
  ; Calculate [S II] 6716/6731 ratio (density sensitive)
  idx_6716 = WHERE(xwav GE 6716.4d - 3.0d*sigma AND xwav LE 6716.4d + 3.0d*sigma, count_6716)
  idx_6731 = WHERE(xwav GE 6730.8d - 3.0d*sigma AND xwav LE 6730.8d + 3.0d*sigma, count_6731)
  
  IF count_6716 GT 2 AND count_6731 GT 2 THEN BEGIN
    int_6716_single = IPT_SIMPSON_INTEGRATE(xwav[idx_6716], ypts_single_maxwellian[idx_6716])
    int_6731_single = IPT_SIMPSON_INTEGRATE(xwav[idx_6731], ypts_single_maxwellian[idx_6731])
    ratio_sii_single = int_6716_single / int_6731_single
    
    int_6716_double = IPT_SIMPSON_INTEGRATE(xwav[idx_6716], ypts_double_maxwellian[idx_6716])
    int_6731_double = IPT_SIMPSON_INTEGRATE(xwav[idx_6731], ypts_double_maxwellian[idx_6731])
    ratio_sii_double = int_6716_double / int_6731_double
    
    PRINT, '[S II] 6716/6731 ratio (density diagnostic):'
    PRINT, '  Single Maxwellian: ', STRTRIM(STRING(ratio_sii_single, FORMAT='(F8.4)'),2)
    PRINT, '  Double Maxwellian: ', STRTRIM(STRING(ratio_sii_double, FORMAT='(F8.4)'),2)
    PRINT, '  (Low density limit ~1.5, high density limit ~0.44)'
  ENDIF
  
  ; Calculate [O II] 3726/3729 ratio (also density sensitive)
  idx_3726 = WHERE(xwav GE 3726.03d - 3.0d*sigma AND xwav LE 3726.03d + 3.0d*sigma, count_3726)
  idx_3729 = WHERE(xwav GE 3728.82d - 3.0d*sigma AND xwav LE 3728.82d + 3.0d*sigma, count_3729)
  
  IF count_3726 GT 2 AND count_3729 GT 2 THEN BEGIN
    int_3726_single = IPT_SIMPSON_INTEGRATE(xwav[idx_3726], ypts_single_maxwellian[idx_3726])
    int_3729_single = IPT_SIMPSON_INTEGRATE(xwav[idx_3729], ypts_single_maxwellian[idx_3729])
    ratio_oii_single = int_3726_single / int_3729_single
    
    int_3726_double = IPT_SIMPSON_INTEGRATE(xwav[idx_3726], ypts_double_maxwellian[idx_3726])
    int_3729_double = IPT_SIMPSON_INTEGRATE(xwav[idx_3729], ypts_double_maxwellian[idx_3729])
    ratio_oii_double = int_3726_double / int_3729_double
    
    PRINT, ''
    PRINT, '[O II] 3726/3729 ratio (density diagnostic):'
    PRINT, '  Single Maxwellian: ', STRTRIM(STRING(ratio_oii_single, FORMAT='(F8.4)'),2)
    PRINT, '  Double Maxwellian: ', STRTRIM(STRING(ratio_oii_double, FORMAT='(F8.4)'),2)
  ENDIF
  
  ; Calculate [O III] temperature diagnostic ratio
  idx_4363 = WHERE(xwav GE 4363.21d - 3.0d*sigma AND xwav LE 4363.21d + 3.0d*sigma, count_4363)
  idx_4959 = WHERE(xwav GE 4958.91d - 3.0d*sigma AND xwav LE 4958.91d + 3.0d*sigma, count_4959)
  idx_5007 = WHERE(xwav GE 5006.84d - 3.0d*sigma AND xwav LE 5006.84d + 3.0d*sigma, count_5007)
  
  IF count_4363 GT 2 AND count_4959 GT 2 AND count_5007 GT 2 THEN BEGIN
    int_4363_single = IPT_SIMPSON_INTEGRATE(xwav[idx_4363], ypts_single_maxwellian[idx_4363])
    int_4959_single = IPT_SIMPSON_INTEGRATE(xwav[idx_4959], ypts_single_maxwellian[idx_4959])
    int_5007_single = IPT_SIMPSON_INTEGRATE(xwav[idx_5007], ypts_single_maxwellian[idx_5007])
    ratio_oiii_single = int_4363_single / (int_4959_single + int_5007_single)
    
    int_4363_double = IPT_SIMPSON_INTEGRATE(xwav[idx_4363], ypts_double_maxwellian[idx_4363])
    int_4959_double = IPT_SIMPSON_INTEGRATE(xwav[idx_4959], ypts_double_maxwellian[idx_4959])
    int_5007_double = IPT_SIMPSON_INTEGRATE(xwav[idx_5007], ypts_double_maxwellian[idx_5007])
    ratio_oiii_double = int_4363_double / (int_4959_double + int_5007_double)
    
    PRINT, ''
    PRINT, '[O III] 4363/(4959+5007) ratio (temperature diagnostic):'
    PRINT, '  Single Maxwellian: ', STRTRIM(STRING(ratio_oiii_single, FORMAT='(F8.6)'),2)
    PRINT, '  Double Maxwellian: ', STRTRIM(STRING(ratio_oiii_double, FORMAT='(F8.6)'),2)
  ENDIF
  
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