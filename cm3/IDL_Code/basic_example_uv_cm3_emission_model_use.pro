;+
; basic_example_uv_cm3_emission_model_use.pro
;
; Example Script for UV Emission Modeling of Io Plasma Torus
; Direct CHIANTI Implementation
; ===========================================================
;
; This script demonstrates the use of direct CHIANTI database calls for
; calculating UV emission spectra from the Io Plasma Torus. It provides
; the same inputs and outputs as the table-interpolation version but
; calculates emissivities directly from the CHIANTI atomic database.
;
; The example uses CHIANTI 11.0.2 atomic database and demonstrates:
; 1. Direct calculation of emission line emissivities via emiss_calc
; 2. Calculating discrete emission line brightnesses
; 3. Convolving with instrument response functions
; 4. Analyzing key diagnostic emission lines
; 5. Creating publication-quality plots
;
; REQUIRED SOFTWARE:
; - IDL 8.0 or later
; - CHIANTI 11.0.2 atomic database (properly installed and configured)
; - SSW (SolarSoftWare) or standalone CHIANTI installation
;
; OUTPUT:
; - Emission spectra for single and double Maxwellian cases
; - Detailed analysis of key UV emission lines
; - Publication-quality plots saved as PNG files
;
; PHYSICS:
; - Volume emission rates calculated directly from CHIANTI (photons/s/cm^3)
; - Electron impact excitation with proper atomic physics
; - Optically thin plasma approximation (valid for IPT)
; - Line-of-sight integration via column densities
;
; COORDINATE SYSTEMS AND UNITS:
; - Temperatures: electron volts [eV]
; - Densities: particles per cubic centimeter [cm^-3]
; - Column densities: particles per square centimeter [cm^-2]
; - Wavelengths: Angstroms [Angstroms]
; - Emissivities: photons per second per cubic centimeter [photons s^-1 cm^-3]
; - Brightnesses: Rayleighs [R], where 1 R = 10^6 photons s^-1 cm^-2 (4pi sr)^-1
;
; REFERENCES:
; - CHIANTI database: Dere et al. 1997; Del Zanna et al. 2020; Dufresne et al. 2024
; - IPT observations: Steffl et al. 2004a,b; Thomas et al. 2004; Bagenal & Delamere 2011
; - Emission modeling: Nerney et al. 2017, 2020, 2022, 2025a, 2025b
; - Electron distributions: Meyer-Vernet & Moncuquet 1989; Moncuquet et al. 2002
;
; AUTHOR: Edward (Eddie) G. Nerney
; INSTITUTION: Laboratory for Atmospheric and Space Physics, University of Colorado Boulder
; LICENSE: Open source for academic and research use
; VERSION: 1.0
; DATE: November 2025
;
; CHIANTI ACKNOWLEDGMENT:
; CHIANTI is a collaborative project involving George Mason University, the University
; of Michigan (USA), University of Cambridge (UK), and NASA Goddard Space Flight Center (USA).
;-

;==============================================================================
FUNCTION IPT_CHIANTI_CALCULATE_EMISSION_SINGLE, temperature, density, $
  col_Sp, col_S2p, col_S3p, col_S4p, col_Op, col_O2p, $
  MIN_WAV=min_wav, MAX_WAV=max_wav, WAVELENGTHS=wavelengths
  ;+
  ; NAME:
  ;   IPT_CHIANTI_CALCULATE_EMISSION_SINGLE
  ;
  ; PURPOSE:
  ;   Calculate IPT emission line brightnesses for a single Maxwellian electron
  ;   distribution by directly calling the CHIANTI atomic database. This function
  ;   provides the same outputs as IPT_CALCULATE_EMISSION_SINGLE but computes
  ;   emissivities on-the-fly rather than interpolating from pre-computed tables.
  ;
  ; CATEGORY:
  ;   Io Plasma Torus, Emission Modeling, Atomic Physics
  ;
  ; CALLING SEQUENCE:
  ;   brightnesses = IPT_CHIANTI_CALCULATE_EMISSION_SINGLE(temperature, density, $
  ;                    col_Sp, col_S2p, col_S3p, col_S4p, col_Op, col_O2p, $
  ;                    MIN_WAV=min_wav, MAX_WAV=max_wav, WAVELENGTHS=wavelengths)
  ;
  ; INPUTS:
  ;   temperature - Electron temperature [eV] (scalar)
  ;   density     - Electron density [cm^-3] (scalar)
  ;   col_Sp      - S+ (S II) column density [cm^-2]
  ;   col_S2p     - S++ (S III) column density [cm^-2]
  ;   col_S3p     - S+++ (S IV) column density [cm^-2]
  ;   col_S4p     - S++++ (S V) column density [cm^-2]
  ;   col_Op      - O+ (O II) column density [cm^-2]
  ;   col_O2p     - O++ (O III) column density [cm^-2]
  ;
  ; KEYWORDS:
  ;   MIN_WAV     - Minimum wavelength [Angstroms], default=550
  ;   MAX_WAV     - Maximum wavelength [Angstroms], default=2100
  ;   WAVELENGTHS - OUTPUT: Array of emission line wavelengths [Angstroms]
  ;
  ; OUTPUTS:
  ;   Returns 1D array of line brightnesses in Rayleighs for all emission lines
  ;   within the specified wavelength range, sorted by wavelength.
  ;
  ; CHIANTI PARAMETERS USED:
  ;   /no_de     - Returns emissivities in photon units (drops hc/lambda factor)
  ;   radt=1.d   - Negligible background radiation (cold space environment)
  ;   /quiet     - Suppress informational messages
  ;   /NOPROT    - Exclude proton collision rates (appropriate for IPT)
  ;   /NOIONREC  - Use fixed ion fractions (ionization equilibrium assumed separate)
  ;   /NO_RREC   - Exclude radiative recombination contributions
  ;
  ; PHYSICS:
  ;   The brightness of an emission line in Rayleighs is:
  ;     I = epsilon(Te, ne) * N_ion * 10^-6  [Rayleighs]
  ;
  ;   where:
  ;   - epsilon(Te, ne) is the volume emissivity [photons s^-1 cm^-3]
  ;   - N_ion is the ion column density [cm^-2]
  ;   - 10^-6 converts to Rayleighs (1 R = 10^6 photons s^-1 cm^-2 (4pi sr)^-1)
  ;
  ; MODIFICATION HISTORY:
  ;   Written by E. G. Nerney, November 2025
  ;-
  COMPILE_OPT IDL2

  ; Set default wavelength range for Europa-UVS/JUICE-UVS instruments
  IF ~KEYWORD_SET(max_wav) THEN max_wav = 2100.0d
  IF ~KEYWORD_SET(min_wav) THEN min_wav = 550.0d

  PRINT, 'Calculating single Maxwellian emission using direct CHIANTI calls...'
  PRINT, '  Te = ', STRTRIM(STRING(temperature, FORMAT='(F6.2)'),2), ' eV'
  PRINT, '  ne = ', STRTRIM(STRING(density, FORMAT='(F8.1)'),2), ' cm^-3'

  ; ============================================================================
  ; UNIT CONVERSION
  ; ============================================================================
  ; CHIANTI requires temperature in log10(Kelvin) and density in log10(cm^-3)
  ; Conversion factor from eV to Kelvin (NIST CODATA 2018)
  eV_to_K = 11604.51812d

  ; Convert temperature to log10(Kelvin) for CHIANTI
  log10_temp = ALOG10(temperature * eV_to_K)

  ; Convert density to log10(cm^-3) for CHIANTI
  log10_dens = ALOG10(density)

  ; ============================================================================
  ; CALCULATE EMISSIVITIES FOR SULFUR IONS
  ; ============================================================================
  ; Atomic number for Sulfur = 16
  ; Ion stages: 1=neutral, 2=S+, 3=S++, 4=S+++, 5=S++++
  ;
  ; CHIANTI emiss_calc parameters:
  ;   /no_de     - Photon units (essential for emission measure analyses)
  ;   radt=1.d   - Background radiation temperature = 1 K (negligible)
  ;   /quiet     - Suppress diagnostic messages
  ;   /NOPROT    - No proton collisions (minor effect in IPT)
  ;   /NOIONREC  - No ionization/recombination (ion fractions fixed)
  ;   /NO_RREC   - No radiative recombination (continuum contribution)

  PRINT, '  Calling CHIANTI for sulfur ions...'

  ; S+ (S II) - singly ionized sulfur
  s2em = emiss_calc(16, 2, temp=log10_temp, dens=log10_dens, $
    /no_de, radt=1.d, /quiet, /NOPROT, /NOIONREC, /NO_RREC)

  ; S++ (S III) - doubly ionized sulfur (dominant sulfur ion in IPT)
  s3em = emiss_calc(16, 3, temp=log10_temp, dens=log10_dens, $
    /no_de, radt=1.d, /quiet, /NOPROT, /NOIONREC, /NO_RREC)

  ; S+++ (S IV) - triply ionized sulfur
  s4em = emiss_calc(16, 4, temp=log10_temp, dens=log10_dens, $
    /no_de, radt=1.d, /quiet, /NOPROT, /NOIONREC, /NO_RREC)

  ; S++++ (S V) - quadruply ionized sulfur
  s5em = emiss_calc(16, 5, temp=log10_temp, dens=log10_dens, $
    /no_de, radt=1.d, /quiet, /NOPROT, /NOIONREC, /NO_RREC)

  ; ============================================================================
  ; CALCULATE EMISSIVITIES FOR OXYGEN IONS
  ; ============================================================================
  ; Atomic number for Oxygen = 8
  ; Ion stages: 1=neutral, 2=O+, 3=O++

  PRINT, '  Calling CHIANTI for oxygen ions...'

  ; O+ (O II) - singly ionized oxygen (dominant ion overall in IPT)
  o2em = emiss_calc(8, 2, temp=log10_temp, dens=log10_dens, $
    /no_de, radt=1.d, /quiet, /NOPROT, /NOIONREC, /NO_RREC)

  ; O++ (O III) - doubly ionized oxygen
  o3em = emiss_calc(8, 3, temp=log10_temp, dens=log10_dens, $
    /no_de, radt=1.d, /quiet, /NOPROT, /NOIONREC, /NO_RREC)

  ; ============================================================================
  ; EXTRACT WAVELENGTHS AND EMISSIVITIES
  ; ============================================================================
  ; The emiss_calc output structure contains:
  ;   .lambda - Array of emission line wavelengths [Angstroms]
  ;   .em     - Array of emissivities [photons s^-1 cm^-3]
  ;             For single Te/ne: shape is (1, 1, n_lines)

  ; Extract wavelengths for each species
  xwavi_sp = s2em.lambda
  xwavi_s2p = s3em.lambda
  xwavi_s3p = s4em.lambda
  xwavi_s4p = s5em.lambda
  xwavi_op = o2em.lambda
  xwavi_o2p = o3em.lambda

  ; Extract emissivities and reshape from (1,1,n_lines) to (n_lines)
  emiss_sp = REFORM(s2em.em)
  emiss_s2p = REFORM(s3em.em)
  emiss_s3p = REFORM(s4em.em)
  emiss_s4p = REFORM(s5em.em)
  emiss_op = REFORM(o2em.em)
  emiss_o2p = REFORM(o3em.em)

  ; ============================================================================
  ; CONVERT TO OBSERVABLE BRIGHTNESSES
  ; ============================================================================
  ; Multiply volume emission rate by column density and convert to Rayleighs
  ; Emissivities from CHIANTI are in photons/s/cm^3
  ; 1 Rayleigh = 10^6 photons/s/cm^2/4pi sr (hence the 1d-6 factor)

  sp_brightness = (1d-6) * emiss_sp * col_Sp      ; S II lines
  s2p_brightness = (1d-6) * emiss_s2p * col_S2p   ; S III lines
  s3p_brightness = (1d-6) * emiss_s3p * col_S3p   ; S IV lines
  s4p_brightness = (1d-6) * emiss_s4p * col_S4p   ; S V lines
  op_brightness = (1d-6) * emiss_op * col_Op      ; O II lines
  o2p_brightness = (1d-6) * emiss_o2p * col_O2p   ; O III lines

  ; ============================================================================
  ; COMBINE WAVELENGTHS AND BRIGHTNESSES
  ; ============================================================================
  ; Concatenate all species into single arrays
  xwavi_all = [xwavi_sp, xwavi_s2p, xwavi_s3p, xwavi_s4p, xwavi_op, xwavi_o2p]
  yptsi_all = [sp_brightness, s2p_brightness, s3p_brightness, $
    s4p_brightness, op_brightness, o2p_brightness]

  ; Sort by wavelength for proper spectral ordering
  wsort = SORT(xwavi_all)
  xwavi_all = xwavi_all[wsort]
  yptsi_all = yptsi_all[wsort]

  ; ============================================================================
  ; FILTER TO REQUESTED WAVELENGTH RANGE
  ; ============================================================================
  avgwav = (min_wav + max_wav)/2.0d
  wrange = WHERE(ABS(xwavi_all - avgwav) LE (avgwav - min_wav), nlines)

  IF nlines GT 0 THEN BEGIN
    wavelengths = xwavi_all[wrange]
    brightnesses = yptsi_all[wrange]
    PRINT, '  Number of emission lines in range: ', STRTRIM(STRING(nlines),2)
  ENDIF ELSE BEGIN
    ; Return empty arrays if no lines in range
    wavelengths = [0d]
    brightnesses = [0d]
    PRINT, '  Warning: No emission lines found in wavelength range'
  ENDELSE

  RETURN, brightnesses
END

;==============================================================================
FUNCTION IPT_CHIANTI_CALCULATE_EMISSION_DOUBLE, core_temp, hot_temp, $
  total_density, hot_fraction, $
  col_Sp, col_S2p, col_S3p, col_S4p, col_Op, col_O2p, $
  MIN_WAV=min_wav, MAX_WAV=max_wav, WAVELENGTHS=wavelengths
  ;+
  ; NAME:
  ;   IPT_CHIANTI_CALCULATE_EMISSION_DOUBLE
  ;
  ; PURPOSE:
  ;   Calculate IPT emission line brightnesses for a double Maxwellian electron
  ;   distribution by directly calling the CHIANTI atomic database. This function
  ;   properly accounts for the nonlinear enhancement of high-excitation transitions
  ;   by suprathermal electrons using CHIANTI's sum_mwl_coeff capability.
  ;
  ; CATEGORY:
  ;   Io Plasma Torus, Emission Modeling, Atomic Physics
  ;
  ; CALLING SEQUENCE:
  ;   brightnesses = IPT_CHIANTI_CALCULATE_EMISSION_DOUBLE(core_temp, hot_temp, $
  ;                    total_density, hot_fraction, col_Sp, col_S2p, col_S3p, $
  ;                    col_S4p, col_Op, col_O2p, MIN_WAV=min_wav, MAX_WAV=max_wav, $
  ;                    WAVELENGTHS=wavelengths)
  ;
  ; INPUTS:
  ;   core_temp      - Core (cold) electron temperature [eV] (scalar)
  ;   hot_temp       - Hot electron temperature [eV] (scalar)
  ;   total_density  - Total electron density [cm^-3] (scalar)
  ;   hot_fraction   - Fraction of hot electrons (0 < feh < 1) (scalar)
  ;   col_Sp         - S+ column density [cm^-2]
  ;   col_S2p        - S++ column density [cm^-2]
  ;   col_S3p        - S+++ column density [cm^-2]
  ;   col_S4p        - S++++ column density [cm^-2]
  ;   col_Op         - O+ column density [cm^-2]
  ;   col_O2p        - O++ column density [cm^-2]
  ;
  ; KEYWORDS:
  ;   MIN_WAV     - Minimum wavelength [Angstroms], default=550
  ;   MAX_WAV     - Maximum wavelength [Angstroms], default=2100
  ;   WAVELENGTHS - OUTPUT: Emission line wavelengths [Angstroms]
  ;
  ; OUTPUTS:
  ;   Returns 1D array of line brightnesses in Rayleighs.
  ;
  ; PHYSICAL MODEL:
  ;   The double Maxwellian electron distribution function is:
  ;     f(v) = (1 - feh) * f_Maxwell(v, Tec) + feh * f_Maxwell(v, Teh)
  ;
  ;   CHIANTI's sum_mwl_coeff keyword properly weights the emissivity
  ;   contributions from each temperature component:
  ;     epsilon_total = (1-feh) * epsilon(Tec, ne) + feh * epsilon(Teh, ne)
  ;
  ;   This arises from:
  ;   - Wave-particle interactions (electron cyclotron waves)
  ;   - Pickup ion acceleration and subsequent electron heating
  ;   - Magnetic reconnection events
  ;   - Centrifugally driven instabilities
  ;
  ;   Effects on emission:
  ;   - High-excitation lines (e.g., S IV, S V) are strongly enhanced
  ;   - Enhancement is nonlinear (not simply proportional to feh)
  ;   - Ratio of high to low excitation lines increases
  ;   - Overall brightness can increase significantly
  ;
  ; TYPICAL IPT PARAMETERS:
  ;   - Core temp: 3-8 eV (from UV line ratios)
  ;   - Hot temp: 100-500 eV (from in situ measurements)
  ;   - Hot fraction: 0.001-0.01 (0.1-1%)
  ;   - Total density: 1000-3000 cm^-3
  ;
  ; REFERENCES:
  ;   - Meyer-Vernet & Moncuquet (1989): Theoretical framework
  ;   - Moncuquet et al. (2002): Voyager observations
  ;   - Nerney et al. (2020, 2022): UV spectroscopic constraints
  ;
  ; MODIFICATION HISTORY:
  ;   Written by E. G. Nerney, November 2025
  ;-
  COMPILE_OPT IDL2

  ; Set default wavelength range for Europa-UVS/JUICE-UVS instruments
  IF ~KEYWORD_SET(max_wav) THEN max_wav = 2100.0d
  IF ~KEYWORD_SET(min_wav) THEN min_wav = 550.0d

  PRINT, 'Calculating double Maxwellian emission using direct CHIANTI calls...'
  PRINT, '  Tec = ', STRTRIM(STRING(core_temp, FORMAT='(F6.2)'),2), ' eV'
  PRINT, '  Teh = ', STRTRIM(STRING(hot_temp, FORMAT='(F6.1)'),2), ' eV'
  PRINT, '  ne_total = ', STRTRIM(STRING(total_density, FORMAT='(F8.1)'),2), ' cm^-3'
  PRINT, '  feh = ', STRTRIM(STRING(hot_fraction, FORMAT='(F8.6)'),2)

  ; ============================================================================
  ; UNIT CONVERSION
  ; ============================================================================
  ; CHIANTI requires temperature in log10(Kelvin) and density in log10(cm^-3)
  ; Conversion factor from eV to Kelvin (NIST CODATA 2018)
  eV_to_K = 11604.51812d

  ; Convert temperatures to log10(Kelvin) for CHIANTI
  log10_tec = ALOG10(core_temp * eV_to_K)
  log10_teh = ALOG10(hot_temp * eV_to_K)

  ; Convert density to log10(cm^-3) for CHIANTI
  log10_dens = ALOG10(total_density)

  ; ============================================================================
  ; SETUP DOUBLE MAXWELLIAN PARAMETERS
  ; ============================================================================
  ; CHIANTI's sum_mwl_coeff keyword allows weighted sum of multiple temperatures
  ; The array of temperatures and corresponding fractional weights
  log10_temps = [log10_tec, log10_teh]  ; Array of 2 temperatures
  fec = 1.0d - hot_fraction              ; Cold electron fraction
  fracs = [fec, hot_fraction]            ; Weights that sum to 1.0

  ; ============================================================================
  ; CALCULATE EMISSIVITIES FOR SULFUR IONS
  ; ============================================================================
  ; Atomic number for Sulfur = 16
  ; Ion stages: 1=neutral, 2=S+, 3=S++, 4=S+++, 5=S++++
  ;
  ; CHIANTI emiss_calc parameters for double Maxwellian:
  ;   temp           - Array of log10(T) values for multi-thermal calculation
  ;   sum_mwl_coeff  - Fractional weights for each temperature component
  ;   dens           - Electron density (single value used for all components)
  ;   /no_de         - Photon units
  ;   radt=1.d       - Negligible background radiation
  ;   /quiet         - Suppress diagnostic messages
  ;   /NOPROT        - No proton collisions
  ;   /NOIONREC      - No ionization/recombination
  ;   /NO_RREC       - No radiative recombination

  PRINT, '  Calling CHIANTI for sulfur ions (double Maxwellian)...'

  ; S+ (S II) - singly ionized sulfur
  s2em = emiss_calc(16, 2, temp=log10_temps, sum_mwl_coeff=fracs, $
    dens=log10_dens, /no_de, radt=1.d, /quiet, $
    /NOPROT, /NOIONREC, /NO_RREC)

  ; S++ (S III) - doubly ionized sulfur (dominant sulfur ion in IPT)
  s3em = emiss_calc(16, 3, temp=log10_temps, sum_mwl_coeff=fracs, $
    dens=log10_dens, /no_de, radt=1.d, /quiet, $
    /NOPROT, /NOIONREC, /NO_RREC)

  ; S+++ (S IV) - triply ionized sulfur (enhanced by hot electrons)
  s4em = emiss_calc(16, 4, temp=log10_temps, sum_mwl_coeff=fracs, $
    dens=log10_dens, /no_de, radt=1.d, /quiet, $
    /NOPROT, /NOIONREC, /NO_RREC)

  ; S++++ (S V) - quadruply ionized sulfur (strongly enhanced by hot electrons)
  s5em = emiss_calc(16, 5, temp=log10_temps, sum_mwl_coeff=fracs, $
    dens=log10_dens, /no_de, radt=1.d, /quiet, $
    /NOPROT, /NOIONREC, /NO_RREC)

  ; ============================================================================
  ; CALCULATE EMISSIVITIES FOR OXYGEN IONS
  ; ============================================================================
  ; Atomic number for Oxygen = 8

  PRINT, '  Calling CHIANTI for oxygen ions (double Maxwellian)...'

  ; O+ (O II) - singly ionized oxygen (dominant ion overall in IPT)
  o2em = emiss_calc(8, 2, temp=log10_temps, sum_mwl_coeff=fracs, $
    dens=log10_dens, /no_de, radt=1.d, /quiet, $
    /NOPROT, /NOIONREC, /NO_RREC)

  ; O++ (O III) - doubly ionized oxygen
  o3em = emiss_calc(8, 3, temp=log10_temps, sum_mwl_coeff=fracs, $
    dens=log10_dens, /no_de, radt=1.d, /quiet, $
    /NOPROT, /NOIONREC, /NO_RREC)

  ; ============================================================================
  ; EXTRACT WAVELENGTHS AND EMISSIVITIES
  ; ============================================================================
  ; For sum_mwl_coeff calculation, emiss_calc returns the weighted sum directly
  ; The output .em has shape that needs to be reformed to 1D

  ; Extract wavelengths for each species
  xwavi_sp = s2em.lambda
  xwavi_s2p = s3em.lambda
  xwavi_s3p = s4em.lambda
  xwavi_s4p = s5em.lambda
  xwavi_op = o2em.lambda
  xwavi_o2p = o3em.lambda

  ; Extract emissivities and reshape to 1D arrays
  emiss_sp = REFORM(s2em.em)
  emiss_s2p = REFORM(s3em.em)
  emiss_s3p = REFORM(s4em.em)
  emiss_s4p = REFORM(s5em.em)
  emiss_op = REFORM(o2em.em)
  emiss_o2p = REFORM(o3em.em)

  ; ============================================================================
  ; CONVERT TO OBSERVABLE BRIGHTNESSES
  ; ============================================================================
  ; Multiply volume emission rate by column density and convert to Rayleighs

  sp_brightness = (1d-6) * emiss_sp * col_Sp      ; S II lines
  s2p_brightness = (1d-6) * emiss_s2p * col_S2p   ; S III lines
  s3p_brightness = (1d-6) * emiss_s3p * col_S3p   ; S IV lines
  s4p_brightness = (1d-6) * emiss_s4p * col_S4p   ; S V lines
  op_brightness = (1d-6) * emiss_op * col_Op      ; O II lines
  o2p_brightness = (1d-6) * emiss_o2p * col_O2p   ; O III lines

  ; ============================================================================
  ; COMBINE AND SORT
  ; ============================================================================
  xwavi_all = [xwavi_sp, xwavi_s2p, xwavi_s3p, xwavi_s4p, xwavi_op, xwavi_o2p]
  yptsi_all = [sp_brightness, s2p_brightness, s3p_brightness, $
    s4p_brightness, op_brightness, o2p_brightness]

  ; Sort by wavelength for proper spectral ordering
  wsort = SORT(xwavi_all)
  xwavi_all = xwavi_all[wsort]
  yptsi_all = yptsi_all[wsort]

  ; ============================================================================
  ; FILTER TO REQUESTED WAVELENGTH RANGE
  ; ============================================================================
  avgwav = (min_wav + max_wav)/2.0d
  wrange = WHERE(ABS(xwavi_all - avgwav) LE (avgwav - min_wav), nlines)

  IF nlines GT 0 THEN BEGIN
    wavelengths = xwavi_all[wrange]
    brightnesses = yptsi_all[wrange]
    PRINT, '  Number of emission lines in range: ', STRTRIM(STRING(nlines),2)
  ENDIF ELSE BEGIN
    wavelengths = [0d]
    brightnesses = [0d]
    PRINT, '  Warning: No emission lines found in wavelength range'
  ENDELSE

  RETURN, brightnesses
END

;==============================================================================
FUNCTION IPT_SIMULATE_SPECTRUM_ERF, wavelength_grid, bin_width, $
  line_wavelengths, line_brightnesses, FWHM=fwhm
  ;+
  ; NAME:
  ;   IPT_SIMULATE_SPECTRUM_ERF
  ;
  ; PURPOSE:
  ;   Convolve discrete emission lines with a Gaussian instrument response
  ;   function to create a realistic observed spectrum. Uses analytical error
  ;   function (ERF) formulation for exact integration of Gaussian profiles
  ;   over wavelength bins.
  ;
  ; CATEGORY:
  ;   Io Plasma Torus, Spectral Modeling, Instrument Response
  ;
  ; CALLING SEQUENCE:
  ;   spectrum = IPT_SIMULATE_SPECTRUM_ERF(wavelength_grid, bin_width, $
  ;                line_wavelengths, line_brightnesses, FWHM=fwhm)
  ;
  ; INPUTS:
  ;   wavelength_grid   - Output wavelength grid (bin centers) [Angstroms]
  ;   bin_width         - Width of wavelength bins [Angstroms] (scalar)
  ;   line_wavelengths  - Wavelengths of discrete emission lines [Angstroms]
  ;   line_brightnesses - Integrated brightnesses of emission lines [Rayleighs]
  ;
  ; KEYWORDS:
  ;   FWHM - Full Width at Half Maximum of Gaussian instrument response [Angstroms]
  ;          Default = 6.0 (appropriate for Europa-UVS/JUICE-UVS)
  ;
  ; OUTPUTS:
  ;   Returns 1D array of spectral brightness per unit wavelength at each
  ;   wavelength grid point. Units: [Rayleighs/Angstrom]
  ;
  ; METHOD:
  ;   For each emission line, calculates the exact integral of a Gaussian line
  ;   profile over each wavelength bin using the error function:
  ;
  ;   Flux_bin = (Brightness_line / 2) * [erf(z_upper) - erf(z_lower)]
  ;
  ;   where:
  ;     z_upper = (lambda_upper - lambda_line) / (sigma * sqrt(2))
  ;     z_lower = (lambda_lower - lambda_line) / (sigma * sqrt(2))
  ;     sigma = FWHM / (2 * sqrt(2 * ln(2)))
  ;
  ; MODIFICATION HISTORY:
  ;   Written by E. G. Nerney, November 2025
  ;-
  COMPILE_OPT IDL2

  ; Set default FWHM if not specified (appropriate for Europa-UVS/JUICE-UVS)
  IF ~KEYWORD_SET(fwhm) THEN fwhm = 6.0d

  ; Convert FWHM to the error function scaling parameter
  ; For a Gaussian: FWHM = 2*sqrt(2*ln(2)) * sigma
  ; We need: 1/(sigma*sqrt(2)) for the ERF argument
  rootc = 2.0d * SQRT(ALOG(2.0d)) / fwhm

  ; Get dimensions
  n_bins = N_ELEMENTS(wavelength_grid)
  n_lines = N_ELEMENTS(line_brightnesses)

  ; Find lines with positive brightness
  good_idx = WHERE(line_brightnesses GT 0.0d, n_good)

  ; If no lines have positive brightness, return zero spectrum
  IF n_good EQ 0 THEN RETURN, DBLARR(n_bins)

  ; Extract only lines with positive brightness for efficiency
  active_wavelengths = line_wavelengths[good_idx]
  active_brightnesses = line_brightnesses[good_idx]

  ; Create 2D arrays for fully vectorized computation
  ; Dimensions: (n_good_lines, n_bins)
  line_wl_2d = REBIN(active_wavelengths, n_good, n_bins)
  grid_2d = REBIN(REFORM(wavelength_grid, 1, n_bins), n_good, n_bins)
  bright_2d = REBIN(active_brightnesses, n_good, n_bins)

  ; Calculate wavelength difference from each line to each bin center
  delta_lambda = grid_2d - line_wl_2d

  ; Calculate ERF arguments for upper and lower bin edges
  half_bin = bin_width / 2.0d
  z_upper = (delta_lambda + half_bin) * rootc
  z_lower = (delta_lambda - half_bin) * rootc

  ; Compute the integrated Gaussian contribution over each bin
  contributions = bright_2d * 0.5d * (ERF(z_upper) - ERF(z_lower))

  ; Sum contributions from all emission lines for each wavelength bin
  spectrum = TOTAL(contributions, 1)

  ; Convert to brightness per unit wavelength [Rayleighs/Angstrom]
  spectrum /= bin_width

  RETURN, spectrum
END

;==============================================================================
FUNCTION IPT_SIMPSON_INTEGRATE, x, y
  ;+
  ; NAME:
  ;   IPT_SIMPSON_INTEGRATE
  ;
  ; PURPOSE:
  ;   Perform numerical integration using composite Simpson's rule for
  ;   irregularly spaced data points. Provides accurate integration for
  ;   smooth functions (O(h^4) accuracy).
  ;
  ; CATEGORY:
  ;   Numerical Methods, Integration
  ;
  ; CALLING SEQUENCE:
  ;   integral = IPT_SIMPSON_INTEGRATE(x, y)
  ;
  ; INPUTS:
  ;   x - Array of x values (independent variable), monotonically increasing
  ;   y - Array of y values (function values at each x position)
  ;
  ; OUTPUTS:
  ;   Returns the approximate integral (scalar) of y(x) dx.
  ;
  ; MODIFICATION HISTORY:
  ;   Written by E. G. Nerney, November 2025
  ;-
  COMPILE_OPT IDL2

  ; Get number of data points
  n = N_ELEMENTS(x)

  ; Handle edge cases
  IF n LT 2 THEN RETURN, 0.0d
  IF n EQ 2 THEN RETURN, 0.5d * (y[0] + y[1]) * (x[1] - x[0])

  ; Number of complete Simpson's rule applications (pairs of intervals)
  n_pairs = (n - 1) / 2

  ; Extract point triplets for vectorized computation
  indices = LINDGEN(n_pairs) * 2L

  ; Extract x and y values for all triplets
  x0 = x[indices]
  x1 = x[indices + 1L]
  x2 = x[indices + 2L]

  y0 = y[indices]
  y1 = y[indices + 1L]
  y2 = y[indices + 2L]

  ; Calculate interval widths
  h1 = x1 - x0
  h2 = x2 - x1
  h = h1 + h2

  ; Apply adaptive Simpson's rule formula for irregular spacing
  contributions = (h / 6.0d) * ((2.0d - h2/h1) * y0 + $
    (h * h / (h1 * h2)) * y1 + $
    (2.0d - h1/h2) * y2)

  ; Sum all Simpson's rule contributions
  integral = TOTAL(contributions, /DOUBLE)

  ; Handle remaining interval if number of intervals is odd
  IF (n MOD 2) EQ 0 THEN BEGIN
    integral += 0.5d * (y[n-2] + y[n-1]) * (x[n-1] - x[n-2])
  ENDIF

  RETURN, integral
END

;==============================================================================
PRO basic_example_uv_cm3_emission_model_use
  ;+
  ; NAME:
  ;   basic_example_uv_cm3_emission_model_use
  ;
  ; PURPOSE:
  ;   Main procedure demonstrating UV emission calculations for the Io Plasma
  ;   Torus using direct CHIANTI atomic database calls. This version requires
  ;   CHIANTI to be installed and properly configured.
  ;
  ; DESCRIPTION:
  ;   This procedure:
  ;   1. Sets up typical IPT plasma parameters
  ;   2. Directly calls CHIANTI to calculate UV emission emissivities
  ;   3. Demonstrates both single and double Maxwellian electron distributions
  ;   4. Convolves with instrument response function
  ;   5. Plots the resulting spectra
  ;
  ; OUTPUTS:
  ;   Creates plots showing simulated UV spectra in Rayleighs/Angstrom
  ;   Saves PNG files of the spectra
  ;
  ; REQUIRED SOFTWARE:
  ;   - CHIANTI 11.0.2 atomic database
  ;   - IDL 8.0 or later
  ;   - SSW or standalone CHIANTI installation
  ;
  ; MODIFICATION HISTORY:
  ;   Written by E. G. Nerney, November 2025
  ;-

  ; Record start time for performance measurement
  start_time = SYSTIME(/SECONDS)

  PRINT, '=================================================================='
  PRINT, 'UV EMISSION MODEL FOR IO PLASMA TORUS'
  PRINT, 'Direct CHIANTI Implementation using CHIANTI 11.0.2'
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
  ; SINGLE MAXWELLIAN CALCULATION
  ; ============================================================================
  PRINT, '=================================================================='
  PRINT, 'SINGLE MAXWELLIAN CALCULATION'
  PRINT, '=================================================================='

  ; Calculate emission line brightnesses using direct CHIANTI calls
  yptsi_single = IPT_CHIANTI_CALCULATE_EMISSION_SINGLE(Tec, nec, $
    col_Sp, col_S2p, col_S3p, col_S4p, col_Op, col_O2p, $
    MIN_WAV=min_xwav, MAX_WAV=max_xwav, WAVELENGTHS=xwavi_single)

  ; Convolve with instrument response
  PRINT, 'Convolving with instrument response function...'
  ypts_single_maxwellian = IPT_SIMULATE_SPECTRUM_ERF(xwav, xwav_bin_width, $
    xwavi_single, yptsi_single, FWHM=fwhm)

  ; ============================================================================
  ; DOUBLE MAXWELLIAN CALCULATION
  ; ============================================================================
  PRINT, ''
  PRINT, '=================================================================='
  PRINT, 'DOUBLE MAXWELLIAN CALCULATION'
  PRINT, '=================================================================='

  ; Calculate emission line brightnesses using direct CHIANTI calls
  yptsi_double = IPT_CHIANTI_CALCULATE_EMISSION_DOUBLE(Tec, Teh, $
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

  ; Plot 4: Two-panel comparison
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

    ; Define integration region (+/- 3 sigma from line center)
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
  PRINT, '=================================================================='

  ; Stop to allow interactive examination of results
  STOP

END