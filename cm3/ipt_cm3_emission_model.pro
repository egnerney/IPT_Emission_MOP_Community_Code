;+
; IPT_cm3_emission_model.pro
;
; Io Plasma Torus (IPT) UV/Optical Emission Model Library
; ========================================================
;
; This library provides comprehensive tools for calculating emission spectra from 
; the Io Plasma Torus in Jupiter's magnetosphere using pre-calculated CHIANTI 
; atomic database emission tables. The code supports both single and double 
; Maxwellian electron distributions and includes realistic instrument response 
; modeling for UV and optical wavelength ranges.
;
; PHYSICAL MODEL:
; - Pre-calculated volume emission rates from CHIANTI 11.0.2 (photons/s)
; - Temperature and density dependent emissivities for single Maxwellian distributions
; - Four-dimensional interpolation for double Maxwellian distributions  
; - Optically thin plasma approximation (valid for IPT)
; - Electron impact excitation with proper atomic physics
; - Line-of-sight integration via column densities
;
; APPLICABLE WAVELENGTH RANGES:
; - UV instruments: 550-2100 Angstroms (Europa-UVS, JUICE-UVS, HST/STIS)
; - Optical instruments: 3000-10000 Angstroms (ground-based telescopes)
; - Supports any wavelength range where CHIANTI data is available
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
FUNCTION IPT_INTERPOLATE_EMISSIVITY_2D_SPECIES, temp_eV, dens_cm3, $
  temp_arr, dens_arr, emiss_table_species
  ;+
  ; NAME:
  ;   IPT_INTERPOLATE_EMISSIVITY_2D_SPECIES
  ;
  ; PURPOSE:
  ;   Interpolate 2D emissivity rates for a single species from single
  ;   Maxwellian emission tables at a specific temperature and density.
  ;   Performs bilinear interpolation in log10(Te)-log10(ne) space.
  ;
  ; CATEGORY:
  ;   Io Plasma Torus, Emission Modeling, Atomic Physics
  ;
  ; CALLING SEQUENCE:
  ;   emiss_interp = IPT_INTERPOLATE_EMISSIVITY_2D_SPECIES(temp_eV, dens_cm3, $
  ;                    temp_arr, dens_arr, emiss_table_species)
  ;
  ; INPUTS:
  ;   temp_eV             - Electron temperature [eV] to interpolate at (scalar)
  ;   dens_cm3            - Electron density [cm^-3] to interpolate at (scalar)
  ;   temp_arr            - Temperature grid from tables [eV], shape (n_temp,)
  ;                         Must be monotonically increasing
  ;   dens_arr            - Density grid from tables [cm^-3], shape (n_dens,)
  ;                         Must be monotonically increasing
  ;   emiss_table_species - Emissivity table for one species
  ;                         Shape: (n_temp, n_dens, n_lines_species)
  ;                         Units: [photons s^-1 cm^-3] (volume emission rates)
  ;
  ; OUTPUTS:
  ;   Returns 1D array of interpolated emissivity rates for all emission
  ;   lines of this species. Shape: (n_lines_species,)
  ;   Units: [photons s^-1 cm^-3]
  ;
  ; METHOD:
  ;   Bilinear interpolation in log10(Te)-log10(ne) space. The table grids are
  ;   uniformly spaced in log10 space, so we:
  ;     1. Convert Te and ne to log10 space
  ;     2. Find bracketing indices in the log-spaced grids
  ;     3. Calculate fractional positions between grid points
  ;     4. Extract the 4 corner values for each emission line
  ;     5. Combine corners with bilinear weights
  ;
  ;   For each emission line, bilinear interpolation uses the formula:
  ;     f(x,y) = (1-wx)*(1-wy)*f00 + (1-wx)*wy*f01 + wx*(1-wy)*f10 + wx*wy*f11
  ;
  ;   where f00, f01, f10, f11 are the four corner values and wx, wy are the
  ;   fractional distances in the x (temperature) and y (density) directions.
  ;
  ; PHYSICS NOTES:
  ;   Emissivity scaling with plasma parameters:
  ;   - Collisional excitation: epsilon ~ n_e^2
  ;   - Boltzmann factor: epsilon ~ exp(-Delta_E / kT)
  ;   - Collision rate: epsilon ~ T^(-1/2)
  ;
  ;   The strong nonlinear temperature dependence necessitates interpolation in
  ;   log space to accurately capture power-law behavior between grid points.
  ;
  ;   Volume emissivity units are [photons s^-1 cm^-3]. Multiply by ion column
  ;   density [cm^-2] and convert to Rayleighs (factor of 10^-6) to get
  ;   observable line brightness.
  ;
  ; EXAMPLE:
  ;   ; Interpolate S++ emission at Te=5.0 eV, ne=2200 cm^-3
  ;   emiss_s2p = IPT_INTERPOLATE_EMISSIVITY_2D_SPECIES(5.0d, 2200.0d, $
  ;                 temp_arr, dens_arr, yptsi_struct.S2P)
  ;
  ;   ; Convert to brightness by multiplying by S++ column density
  ;   brightness_rayleighs = emiss_s2p * n_s2p * 1.0d-6
  ;
  ; PERFORMANCE:
  ;   Fully vectorized - interpolates all emission lines simultaneously.
  ;   Approximately 10-50x faster than looping over emission lines.
  ;   Typical execution time: ~0.05 ms for tables with 100-500 emission lines.
  ;
  ; MODIFICATION HISTORY:
  ;   Written by E. G. Nerney, November 2025
  ;-
  COMPILE_OPT IDL2, HIDDEN

  ; Get dimensions of the emission table
  dims = SIZE(emiss_table_species, /DIMENSIONS)
  n_temp = dims[0]
  n_dens = dims[1]
  n_lines = dims[2]

  ; Convert input parameters to log10 space for interpolation
  ; Tables are uniformly sampled in log10(Te) and log10(ne)
  log10_temp = ALOG10(temp_eV)
  log10_dens = ALOG10(dens_cm3)
  log10_temp_arr = ALOG10(temp_arr)
  log10_dens_arr = ALOG10(dens_arr)

  ; Find the bracketing indices in the log-space grid using VALUE_LOCATE
  ; VALUE_LOCATE returns the index i such that array[i] <= value < array[i+1]
  i_temp = VALUE_LOCATE(log10_temp_arr, log10_temp)
  i_dens = VALUE_LOCATE(log10_dens_arr, log10_dens)

  ; Clamp indices to valid interpolation range [0, n-2]
  ; This handles edge cases where the query point is at or beyond grid boundaries
  i_temp = 0 > i_temp < (n_temp - 2)
  i_dens = 0 > i_dens < (n_dens - 2)

  ; Calculate fractional position between grid points
  ; frac = 0 at lower grid point, frac = 1 at upper grid point
  frac_temp = (log10_temp - log10_temp_arr[i_temp]) / $
    (log10_temp_arr[i_temp + 1] - log10_temp_arr[i_temp])
  frac_dens = (log10_dens - log10_dens_arr[i_dens]) / $
    (log10_dens_arr[i_dens + 1] - log10_dens_arr[i_dens])

  ; Clamp fractional positions to [0, 1] for safety
  ; Prevents extrapolation artifacts at grid boundaries
  frac_temp = 0.0d > frac_temp < 1.0d
  frac_dens = 0.0d > frac_dens < 1.0d

  ; Extract the 4 corner values of the 2D cell for all emission lines
  ; Using REFORM ensures we get pure 1D arrays of shape (n_lines,)
  ; Corners are labeled by their position in the (temp, dens) grid:
  ;   v00 = lower temp, lower dens
  ;   v01 = lower temp, upper dens
  ;   v10 = upper temp, lower dens
  ;   v11 = upper temp, upper dens
  v00 = REFORM(emiss_table_species[i_temp,   i_dens,   *], n_lines)
  v01 = REFORM(emiss_table_species[i_temp,   i_dens+1, *], n_lines)
  v10 = REFORM(emiss_table_species[i_temp+1, i_dens,   *], n_lines)
  v11 = REFORM(emiss_table_species[i_temp+1, i_dens+1, *], n_lines)

  ; Perform vectorized bilinear interpolation for all emission lines
  ; Standard bilinear interpolation formula:
  ;   f(x,y) = (1-x)*(1-y)*f00 + (1-x)*y*f01 + x*(1-y)*f10 + x*y*f11
  ; where x = frac_temp and y = frac_dens
  emiss_interp = v00 * (1.0d - frac_temp) * (1.0d - frac_dens) + $
    v01 * (1.0d - frac_temp) * frac_dens + $
    v10 * frac_temp * (1.0d - frac_dens) + $
    v11 * frac_temp * frac_dens

  ; Set any negative values to zero
  ; Negative values can arise from numerical artifacts near table boundaries
  ; or in regions with very low emissivity approaching machine precision
  neg_idx = WHERE(emiss_interp LT 0.0d, n_neg)
  IF n_neg GT 0 THEN emiss_interp[neg_idx] = 0.0d

  RETURN, emiss_interp
END

;==============================================================================
FUNCTION IPT_INTERPOLATE_EMISSIVITY_4D_SPECIES, Tec_eV, Teh_eV, ne_cm3, feh, $
  tec_arr, teh_arr, ne_arr, feh_arr, emiss_table_species
  ;+
  ; NAME:
  ;   IPT_INTERPOLATE_EMISSIVITY_4D_SPECIES
  ;
  ; PURPOSE:
  ;   Interpolate 4D emissivity rates for a single species from double
  ;   Maxwellian emission tables. Performs quadralinear (4D linear)
  ;   interpolation in log10 parameter space.
  ;
  ; CATEGORY:
  ;   Io Plasma Torus, Emission Modeling, Atomic Physics
  ;
  ; CALLING SEQUENCE:
  ;   emiss_interp = IPT_INTERPOLATE_EMISSIVITY_4D_SPECIES(Tec_eV, Teh_eV, $
  ;                    ne_cm3, feh, tec_arr, teh_arr, ne_arr, feh_arr, $
  ;                    emiss_table_species)
  ;
  ; INPUTS:
  ;   Tec_eV              - Core electron temperature [eV] (scalar)
  ;   Teh_eV              - Hot electron temperature [eV] (scalar)
  ;   ne_cm3              - Total electron density [cm^-3] (scalar)
  ;   feh                 - Hot electron fraction, 0 < feh < 1 (scalar)
  ;   tec_arr             - Core temperature grid [eV], shape (n_tec,)
  ;   teh_arr             - Hot temperature grid [eV], shape (n_teh,)
  ;   ne_arr              - Density grid [cm^-3], shape (n_ne,)
  ;   feh_arr             - Hot fraction grid, shape (n_feh,)
  ;   emiss_table_species - Emissivity table for one species
  ;                         Shape: (n_tec, n_teh, n_ne, n_feh, n_lines_species)
  ;
  ; OUTPUTS:
  ;   Returns 1D array of interpolated volume emission rates for all emission
  ;   lines of this species. Shape: (n_lines_species,)
  ;   Units: [photons cm^-3 s^-1]
  ;
  ; METHOD:
  ;   Fully vectorized quadralinear interpolation using pre-computed corner
  ;   weights applied simultaneously to all emission lines via matrix operations.
  ;
  ; MODIFICATION HISTORY:
  ;   Written by E. G. Nerney, November 2025
  ;   Optimized for vectorized matrix operations, November 2025
  ;-
  COMPILE_OPT IDL2, HIDDEN

  ; Get dimensions of the emission table
  dims = SIZE(emiss_table_species, /DIMENSIONS)
  n_tec = dims[0]
  n_teh = dims[1]
  n_ne = dims[2]
  n_feh = dims[3]
  n_lines = dims[4]

  ; Convert all input parameters to log10 space for interpolation
  log10_tec = ALOG10(Tec_eV)
  log10_teh = ALOG10(Teh_eV)
  log10_ne = ALOG10(ne_cm3)
  log10_feh = ALOG10(feh)

  log10_tec_arr = ALOG10(tec_arr)
  log10_teh_arr = ALOG10(teh_arr)
  log10_ne_arr = ALOG10(ne_arr)
  log10_feh_arr = ALOG10(feh_arr)

  ; Find bracketing indices and clamp to valid range
  itec = 0 > VALUE_LOCATE(log10_tec_arr, log10_tec) < (n_tec - 2)
  iteh = 0 > VALUE_LOCATE(log10_teh_arr, log10_teh) < (n_teh - 2)
  ine = 0 > VALUE_LOCATE(log10_ne_arr, log10_ne) < (n_ne - 2)
  ifeh = 0 > VALUE_LOCATE(log10_feh_arr, log10_feh) < (n_feh - 2)

  ; Calculate and clamp interpolation weights
  wtec = 0.0d > ((log10_tec - log10_tec_arr[itec]) / $
    (log10_tec_arr[itec + 1] - log10_tec_arr[itec])) < 1.0d
  wteh = 0.0d > ((log10_teh - log10_teh_arr[iteh]) / $
    (log10_teh_arr[iteh + 1] - log10_teh_arr[iteh])) < 1.0d
  wne = 0.0d > ((log10_ne - log10_ne_arr[ine]) / $
    (log10_ne_arr[ine + 1] - log10_ne_arr[ine])) < 1.0d
  wfeh = 0.0d > ((log10_feh - log10_feh_arr[ifeh]) / $
    (log10_feh_arr[ifeh + 1] - log10_feh_arr[ifeh])) < 1.0d

  ; Pre-compute complementary weights
  w0tec = 1.0d - wtec
  w0teh = 1.0d - wteh
  w0ne = 1.0d - wne
  w0feh = 1.0d - wfeh

  ; Compute all 16 corner weights as a vector
  corner_weights = [ $
    w0tec * w0teh * w0ne * w0feh, $  ; 0000
    wtec  * w0teh * w0ne * w0feh, $  ; 1000
    w0tec * wteh  * w0ne * w0feh, $  ; 0100
    wtec  * wteh  * w0ne * w0feh, $  ; 1100
    w0tec * w0teh * wne  * w0feh, $  ; 0010
    wtec  * w0teh * wne  * w0feh, $  ; 1010
    w0tec * wteh  * wne  * w0feh, $  ; 0110
    wtec  * wteh  * wne  * w0feh, $  ; 1110
    w0tec * w0teh * w0ne * wfeh,  $  ; 0001
    wtec  * w0teh * w0ne * wfeh,  $  ; 1001
    w0tec * wteh  * w0ne * wfeh,  $  ; 0101
    wtec  * wteh  * w0ne * wfeh,  $  ; 1101
    w0tec * w0teh * wne  * wfeh,  $  ; 0011
    wtec  * w0teh * wne  * wfeh,  $  ; 1011
    w0tec * wteh  * wne  * wfeh,  $  ; 0111
    wtec  * wteh  * wne  * wfeh   $  ; 1111
    ]

  ; Build corner values array: shape (16, n_lines)
  ; Using REFORM to ensure proper 1D shape for each corner
  corners = DBLARR(16, n_lines)
  corners[0,*]  = REFORM(emiss_table_species[itec,   iteh,   ine,   ifeh,   *], n_lines)
  corners[1,*]  = REFORM(emiss_table_species[itec+1, iteh,   ine,   ifeh,   *], n_lines)
  corners[2,*]  = REFORM(emiss_table_species[itec,   iteh+1, ine,   ifeh,   *], n_lines)
  corners[3,*]  = REFORM(emiss_table_species[itec+1, iteh+1, ine,   ifeh,   *], n_lines)
  corners[4,*]  = REFORM(emiss_table_species[itec,   iteh,   ine+1, ifeh,   *], n_lines)
  corners[5,*]  = REFORM(emiss_table_species[itec+1, iteh,   ine+1, ifeh,   *], n_lines)
  corners[6,*]  = REFORM(emiss_table_species[itec,   iteh+1, ine+1, ifeh,   *], n_lines)
  corners[7,*]  = REFORM(emiss_table_species[itec+1, iteh+1, ine+1, ifeh,   *], n_lines)
  corners[8,*]  = REFORM(emiss_table_species[itec,   iteh,   ine,   ifeh+1, *], n_lines)
  corners[9,*]  = REFORM(emiss_table_species[itec+1, iteh,   ine,   ifeh+1, *], n_lines)
  corners[10,*] = REFORM(emiss_table_species[itec,   iteh+1, ine,   ifeh+1, *], n_lines)
  corners[11,*] = REFORM(emiss_table_species[itec+1, iteh+1, ine,   ifeh+1, *], n_lines)
  corners[12,*] = REFORM(emiss_table_species[itec,   iteh,   ine+1, ifeh+1, *], n_lines)
  corners[13,*] = REFORM(emiss_table_species[itec+1, iteh,   ine+1, ifeh+1, *], n_lines)
  corners[14,*] = REFORM(emiss_table_species[itec,   iteh+1, ine+1, ifeh+1, *], n_lines)
  corners[15,*] = REFORM(emiss_table_species[itec+1, iteh+1, ine+1, ifeh+1, *], n_lines)

  ; Matrix multiply: weights (1x16) # corners (16 x n_lines) -> (1 x n_lines)
  ; REFORM to get pure 1D array output
  emiss_interp = REFORM(corner_weights # corners, n_lines)

  ; Clamp negative values to zero
  neg_idx = WHERE(emiss_interp LT 0.0d, n_neg)
  IF n_neg GT 0 THEN emiss_interp[neg_idx] = 0.0d

  RETURN, emiss_interp
END
;==============================================================================
FUNCTION IPT_LOAD_SINGLE_MAXWELLIAN_TABLES, filename
;+
; NAME:
;   IPT_LOAD_SINGLE_MAXWELLIAN_TABLES
;
; PURPOSE:
;   Load single Maxwellian emission tables from HDF5 file and organize
;   by species for efficient interpolation.
;
; CATEGORY:
;   Io Plasma Torus, Data I/O, Emission Modeling
;
; CALLING SEQUENCE:
;   tables = IPT_LOAD_SINGLE_MAXWELLIAN_TABLES(filename)
;
; INPUTS:
;   filename - Full path to HDF5 file containing single Maxwellian
;              emission tables (e.g., 'CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5')
;
; OUTPUTS:
;   Returns a structure containing:
;     .temp_arr   - Temperature grid [eV], shape (n_temp,)
;     .dens_arr   - Density grid [cm^-3], shape (n_dens,)
;     .xwavi      - Structure of wavelength arrays by species
;                   e.g., .xwavi.SP contains S+ wavelengths
;     .yptsi      - Structure of emissivity tables by species
;                   e.g., .yptsi.SP has shape (n_temp, n_dens, n_lines_sp)
;
; HDF5 FILE STRUCTURE:
;   The HDF5 file contains:
;   - T : temperature grid [eV], shape (n_T,)
;   - n : density grid [cm^-3], shape (n_n,)
;   - emissivity : volume emission rates [photons cm^-3 s^-1], 
;                  shape (n_T, n_n, n_lines)
;   - wavelength : emission line wavelengths [Angstroms], shape (n_lines,)
;   - species : species labels, shape (n_lines,)
;
; SPECIES NOTATION:
;   - 'S'    : S I (neutral sulfur) - negligible in IPT
;   - 'SP'   : S II (S+) - singly ionized sulfur
;   - 'S2P'  : S III (S++) - doubly ionized sulfur (dominant S ion)
;   - 'S3P'  : S IV (S+++) - triply ionized sulfur
;   - 'S4P'  : S V (S++++) - quadruply ionized sulfur
;   - 'O'    : O I (neutral oxygen) - negligible in IPT
;   - 'OP'   : O II (O+) - singly ionized oxygen (dominant ion overall)
;   - 'O2P'  : O III (O++) - doubly ionized oxygen
;   - 'NAP'  : Na II (Na+) - trace species from Io
;
; NOTES:
;   Due to IDL's column-major ordering vs Python's row-major ordering,
;   the HDF5 arrays are transposed after reading to restore the correct
;   dimension ordering for IDL.
;
; EXAMPLE:
;   tables = IPT_LOAD_SINGLE_MAXWELLIAN_TABLES('path/to/tables.h5')
;   PRINT, 'Temperature range:', MIN(tables.temp_arr), '-', MAX(tables.temp_arr), ' eV'
;
; MODIFICATION HISTORY:
;   Written by E. G. Nerney, November 2025
;-

  PRINT, 'Loading single Maxwellian emission tables...'
  PRINT, '  File: ' + filename
  
  ; Open HDF5 file
  file_id = H5F_OPEN(filename)
  
  ; Read temperature and density grids
  temp_arr = H5D_READ(H5D_OPEN(file_id, 'T'))
  dens_arr = H5D_READ(H5D_OPEN(file_id, 'n'))
  
  ; Read wavelength and species arrays
  wavelength_all = H5D_READ(H5D_OPEN(file_id, 'wavelength'))
  species_all = H5D_READ(H5D_OPEN(file_id, 'species'))
  
  ; Read volume emission rate table
  ; Python stored as (n_T, n_n, n_lines) in row-major
  ; IDL reads as (n_lines, n_n, n_T) due to column-major ordering
  ; Need to TRANSPOSE to get back to (n_T, n_n, n_lines)
  vol_emiss_rate_all = H5D_READ(H5D_OPEN(file_id, 'emiss'))
  vol_emiss_rate_all = TRANSPOSE(vol_emiss_rate_all, [2, 1, 0])
  
  ; Close HDF5 file
  H5F_CLOSE, file_id
  
  PRINT, '  HDF5 file read successfully'
  PRINT, '  Temperature range: ', STRTRIM(STRING(MIN(temp_arr), FORMAT='(F8.3)'),2), $
         ' - ', STRTRIM(STRING(MAX(temp_arr), FORMAT='(F8.1)'),2), ' eV'
  PRINT, '  Density range: ', STRTRIM(STRING(MIN(dens_arr), FORMAT='(F8.1)'),2), $
         ' - ', STRTRIM(STRING(MAX(dens_arr), FORMAT='(F10.1)'),2), ' cm^-3'
  PRINT, '  Grid size: ', STRTRIM(STRING(N_ELEMENTS(temp_arr)),2), ' x ', $
         STRTRIM(STRING(N_ELEMENTS(dens_arr)),2)
  PRINT, '  Total emission lines: ', STRTRIM(STRING(N_ELEMENTS(wavelength_all)),2)
  
  ; Reconstruct YPTSI and XWAVI structures by species
  ; This allows species-by-species interpolation which is essential for
  ; maintaining proper line-species associations
  species_names = ['S', 'SP', 'S2P', 'S3P', 'S4P', 'O', 'OP', 'O2P', 'NAP']
  
  ; Initialize structures with dummy values
  xwavi = CREATE_STRUCT('DUMMY', 0.0d)
  yptsi = CREATE_STRUCT('DUMMY', DBLARR(1,1,1))
  first_valid = 1
  
  PRINT, '  Organizing by species:'
  FOR i = 0, N_ELEMENTS(species_names)-1 DO BEGIN
    species_name = species_names[i]
    
    ; Find indices for this species
    idx = WHERE(species_all EQ species_name, count)
    
    IF count GT 0 THEN BEGIN
      ; Extract wavelengths for this species
      waves = wavelength_all[idx]
      
      ; Extract volume emission rates for this species
      ; Shape: (n_T, n_n, n_lines_for_species)
      emiss = vol_emiss_rate_all[*, *, idx]
      
      ; Add to structures
      field_name = STRUPCASE(species_name)
      IF first_valid EQ 1 THEN BEGIN
        xwavi = CREATE_STRUCT(field_name, waves)
        yptsi = CREATE_STRUCT(field_name, emiss)
        first_valid = 0
      ENDIF ELSE BEGIN
        xwavi = CREATE_STRUCT(xwavi, field_name, waves)
        yptsi = CREATE_STRUCT(yptsi, field_name, emiss)
      ENDELSE
      
      PRINT, '    ', species_name, ': ', STRTRIM(STRING(count),2), ' lines'
    ENDIF
  ENDFOR
  
  ; Create output structure
  tables = {temp_arr: temp_arr, $
            dens_arr: dens_arr, $
            xwavi: xwavi, $
            yptsi: yptsi}
  
  PRINT, 'Single Maxwellian tables loaded successfully'
  PRINT, ''
  
  RETURN, tables
END

;==============================================================================
FUNCTION IPT_LOAD_DOUBLE_MAXWELLIAN_TABLES, filename
;+
; NAME:
;   IPT_LOAD_DOUBLE_MAXWELLIAN_TABLES
;
; PURPOSE:
;   Load double Maxwellian emission tables from HDF5 file and organize
;   by species for efficient 4D interpolation.
;
; CATEGORY:
;   Io Plasma Torus, Data I/O, Emission Modeling
;
; CALLING SEQUENCE:
;   tables = IPT_LOAD_DOUBLE_MAXWELLIAN_TABLES(filename)
;
; INPUTS:
;   filename - Full path to HDF5 file containing double Maxwellian
;              emission tables
;
; OUTPUTS:
;   Returns a structure containing:
;     .tec_arr   - Core temperature grid [eV], shape (n_tec,)
;     .teh_arr   - Hot temperature grid [eV], shape (n_teh,)
;     .ne_arr    - Total density grid [cm^-3], shape (n_ne,)
;     .feh_arr   - Hot electron fraction grid, shape (n_feh,)
;     .xwavi     - Structure of wavelength arrays by species
;     .yptsi     - Structure of emissivity tables by species
;                  e.g., .yptsi.SP has shape (n_tec, n_teh, n_ne, n_feh, n_lines_sp)
;
; HDF5 FILE STRUCTURE:
;   - T_cold : core temperature grid [eV], shape (n_Tc,)
;   - T_hot : hot temperature grid [eV], shape (n_Th,)
;   - n : total density grid [cm^-3], shape (n_n,)
;   - feh : hot electron fraction grid, shape (n_feh,)
;   - emissivity : volume emission rates [photons cm^-3 s^-1],
;                  shape (n_Tc, n_Th, n_n, n_feh, n_lines)
;   - wavelength : emission line wavelengths [Angstroms], shape (n_lines,)
;   - species : species labels, shape (n_lines,)
;
; PHYSICAL MODEL:
;   The double Maxwellian electron distribution function is:
;     f(v) = (1 - feh) * f_Maxwell(v, Tec) + feh * f_Maxwell(v, Teh)
;   
;   where:
;   - Tec is the core electron temperature
;   - Teh is the hot electron temperature (suprathermal tail)
;   - feh is the hot electron fraction
;
;   Typical IPT Parameters:
;   - Core temp: 3-8 eV (from UV line ratios)
;   - Hot temp: 100-500 eV (from in situ measurements)
;   - Hot fraction: 0.001-0.01 (0.1-1%)
;   - Total density: 1000-3000 cm^-3
;
; EXAMPLE:
;   tables = IPT_LOAD_DOUBLE_MAXWELLIAN_TABLES('path/to/tables.h5')
;   PRINT, 'Hot fraction range:', MIN(tables.feh_arr), '-', MAX(tables.feh_arr)
;
; MODIFICATION HISTORY:
;   Written by E. G. Nerney, November 2025
;-

  PRINT, 'Loading double Maxwellian emission tables...'
  PRINT, '  File: ' + filename
  
  ; Open HDF5 file
  file_id = H5F_OPEN(filename)
  
  ; Read temperature, density, and fraction grids
  tec_arr = H5D_READ(H5D_OPEN(file_id, 'T_cold'))
  teh_arr = H5D_READ(H5D_OPEN(file_id, 'T_hot'))
  ne_arr = H5D_READ(H5D_OPEN(file_id, 'n'))
  feh_arr = H5D_READ(H5D_OPEN(file_id, 'feh'))
  
  ; Read wavelength and species arrays
  wavelength_all = H5D_READ(H5D_OPEN(file_id, 'wavelength'))
  species_all = H5D_READ(H5D_OPEN(file_id, 'species'))
  
  ; Read volume emission rate table
  ; Python stored as (n_Tc, n_Th, n_n, n_feh, n_lines) in row-major
  ; IDL reads as (n_lines, n_feh, n_n, n_Th, n_Tc) due to column-major ordering
  ; Need to TRANSPOSE to get back to (n_Tc, n_Th, n_n, n_feh, n_lines)
  vol_emiss_rate_all = H5D_READ(H5D_OPEN(file_id, 'emiss'))
  vol_emiss_rate_all = TRANSPOSE(vol_emiss_rate_all, [4, 3, 2, 1, 0])
  
  ; Close HDF5 file
  H5F_CLOSE, file_id
  
  PRINT, '  HDF5 file read successfully'
  PRINT, '  Core temperature range: ', STRTRIM(STRING(MIN(tec_arr), FORMAT='(F8.3)'),2), $
         ' - ', STRTRIM(STRING(MAX(tec_arr), FORMAT='(F8.1)'),2), ' eV'
  PRINT, '  Hot temperature range: ', STRTRIM(STRING(MIN(teh_arr), FORMAT='(F8.1)'),2), $
         ' - ', STRTRIM(STRING(MAX(teh_arr), FORMAT='(F8.1)'),2), ' eV'
  PRINT, '  Density range: ', STRTRIM(STRING(MIN(ne_arr), FORMAT='(F8.1)'),2), $
         ' - ', STRTRIM(STRING(MAX(ne_arr), FORMAT='(F10.1)'),2), ' cm^-3'
  PRINT, '  Hot fraction range: ', STRTRIM(STRING(MIN(feh_arr), FORMAT='(E10.4)'),2), $
         ' - ', STRTRIM(STRING(MAX(feh_arr), FORMAT='(F8.4)'),2)
  PRINT, '  Grid size: ', STRTRIM(STRING(N_ELEMENTS(tec_arr)),2), ' x ', $
         STRTRIM(STRING(N_ELEMENTS(teh_arr)),2), ' x ', $
         STRTRIM(STRING(N_ELEMENTS(ne_arr)),2), ' x ', $
         STRTRIM(STRING(N_ELEMENTS(feh_arr)),2)
  PRINT, '  Total emission lines: ', STRTRIM(STRING(N_ELEMENTS(wavelength_all)),2)
  
  ; Reconstruct YPTSI and XWAVI structures by species
  species_names = ['S', 'SP', 'S2P', 'S3P', 'S4P', 'O', 'OP', 'O2P', 'NAP']
  
  ; Initialize structures
  xwavi = CREATE_STRUCT('DUMMY', 0.0d)
  yptsi = CREATE_STRUCT('DUMMY', DBLARR(1,1,1,1,1))
  first_valid = 1
  
  PRINT, '  Organizing by species:'
  FOR i = 0, N_ELEMENTS(species_names)-1 DO BEGIN
    species_name = species_names[i]
    
    ; Find indices for this species
    idx = WHERE(species_all EQ species_name, count)
    
    IF count GT 0 THEN BEGIN
      ; Extract wavelengths for this species
      waves = wavelength_all[idx]
      
      ; Extract volume emission rates for this species
      ; Shape: (n_Tc, n_Th, n_n, n_feh, n_lines_for_species)
      emiss = vol_emiss_rate_all[*, *, *, *, idx]
      
      ; Add to structures
      field_name = STRUPCASE(species_name)
      IF first_valid EQ 1 THEN BEGIN
        xwavi = CREATE_STRUCT(field_name, waves)
        yptsi = CREATE_STRUCT(field_name, emiss)
        first_valid = 0
      ENDIF ELSE BEGIN
        xwavi = CREATE_STRUCT(xwavi, field_name, waves)
        yptsi = CREATE_STRUCT(yptsi, field_name, emiss)
      ENDELSE
      
      PRINT, '    ', species_name, ': ', STRTRIM(STRING(count),2), ' lines'
    ENDIF
  ENDFOR
  
  ; Create output structure
  tables = {tec_arr: tec_arr, $
            teh_arr: teh_arr, $
            ne_arr: ne_arr, $
            feh_arr: feh_arr, $
            xwavi: xwavi, $
            yptsi: yptsi}
  
  PRINT, 'Double Maxwellian tables loaded successfully'
  PRINT, ''
  
  RETURN, tables
END

;==============================================================================
FUNCTION IPT_CALCULATE_EMISSION_SINGLE, tables, temperature, density, $
  col_Sp, col_S2p, col_S3p, col_S4p, col_Op, col_O2p, $
  MIN_WAV=min_wav, MAX_WAV=max_wav, WAVELENGTHS=wavelengths
;+
; NAME:
;   IPT_CALCULATE_EMISSION_SINGLE
;
; PURPOSE:
;   Calculate IPT emission line brightnesses for a single Maxwellian electron
;   distribution using pre-calculated emission tables. Performs species-by-species
;   interpolation to maintain proper line-species associations.
;
; CATEGORY:
;   Io Plasma Torus, Emission Modeling, Spectral Analysis
;
; CALLING SEQUENCE:
;   brightnesses = IPT_CALCULATE_EMISSION_SINGLE(tables, temperature, density, $
;                    col_Sp, col_S2p, col_S3p, col_S4p, col_Op, col_O2p, $
;                    MIN_WAV=min_wav, MAX_WAV=max_wav, WAVELENGTHS=wavelengths)
;
; INPUTS:
;   tables      - Structure from IPT_LOAD_SINGLE_MAXWELLIAN_TABLES
;   temperature - Electron temperature [eV]
;   density     - Electron density [cm^-3]
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
; PHYSICS:
;   The brightness of an emission line in Rayleighs is:
;     I = epsilon(Te, ne) * N_ion * 10^-6  [Rayleighs]
;   
;   where:
;   - epsilon(Te, ne) is the volume emissivity [photons s^-1 cm^-3]
;   - N_ion is the ion column density [cm^-2]
;   - 10^-6 converts to Rayleighs (1 R = 10^6 photons s^-1 cm^-2 (4pi sr)^-1)
;   
;   This assumes:
;   - Optically thin plasma (no self-absorption)
;   - Uniform temperature and density along line of sight
;   - Ion column density approximates line-of-sight integral
;
; TYPICAL IPT COLUMN DENSITIES:
;   For ~6 R_J path length through the torus:
;   - O+: 5 x 10^13 cm^-2 (dominant ion)
;   - S++: 4 x 10^13 cm^-2 (dominant sulfur ion)
;   - S+: 1 x 10^13 cm^-2
;   - O++, S+++: ~6 x 10^12 cm^-2
;   - S++++: ~6 x 10^11 cm^-2
;
; EXAMPLE:
;   brightnesses = IPT_CALCULATE_EMISSION_SINGLE(tables, 5.0d, 2200.0d, $
;                    1.2d13, 4.2d13, 5.92d12, 6.0d11, 5.2d13, 5.92d12, $
;                    MIN_WAV=550.0d, MAX_WAV=2100.0d, WAVELENGTHS=xwavi)
;
; MODIFICATION HISTORY:
;   Written by E. G. Nerney, November 2025
;-

  ; Set default wavelength range
  IF ~KEYWORD_SET(max_wav) THEN max_wav = 2100.0d
  IF ~KEYWORD_SET(min_wav) THEN min_wav = 550.0d
  
  PRINT, 'Calculating single Maxwellian emission using tables...'
  PRINT, '  Te = ', STRTRIM(STRING(temperature, FORMAT='(F6.2)'),2), ' eV'
  PRINT, '  ne = ', STRTRIM(STRING(density, FORMAT='(F8.1)'),2), ' cm^-3'
  
  ; ============================================================================
  ; INTERPOLATE EMISSIVITIES FOR EACH ION SPECIES
  ; ============================================================================
  ; Each species is interpolated separately to maintain proper line associations
  
  ; S+ (S II)
  emiss_sp = IPT_INTERPOLATE_EMISSIVITY_2D_SPECIES(temperature, density, $
    tables.temp_arr, tables.dens_arr, tables.yptsi.SP)
  
  ; S++ (S III) - dominant sulfur ion in IPT
  emiss_s2p = IPT_INTERPOLATE_EMISSIVITY_2D_SPECIES(temperature, density, $
    tables.temp_arr, tables.dens_arr, tables.yptsi.S2P)
  
  ; S+++ (S IV)
  emiss_s3p = IPT_INTERPOLATE_EMISSIVITY_2D_SPECIES(temperature, density, $
    tables.temp_arr, tables.dens_arr, tables.yptsi.S3P)
  
  ; S++++ (S V)
  emiss_s4p = IPT_INTERPOLATE_EMISSIVITY_2D_SPECIES(temperature, density, $
    tables.temp_arr, tables.dens_arr, tables.yptsi.S4P)
  
  ; O+ (O II) - dominant ion overall in IPT
  emiss_op = IPT_INTERPOLATE_EMISSIVITY_2D_SPECIES(temperature, density, $
    tables.temp_arr, tables.dens_arr, tables.yptsi.OP)
  
  ; O++ (O III)
  emiss_o2p = IPT_INTERPOLATE_EMISSIVITY_2D_SPECIES(temperature, density, $
    tables.temp_arr, tables.dens_arr, tables.yptsi.O2P)
  
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
  xwavi_all = [tables.xwavi.SP, tables.xwavi.S2P, tables.xwavi.S3P, $
               tables.xwavi.S4P, tables.xwavi.OP, tables.xwavi.O2P]
  
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
FUNCTION IPT_CALCULATE_EMISSION_DOUBLE, tables, core_temp, hot_temp, $
  total_density, hot_fraction, $
  col_Sp, col_S2p, col_S3p, col_S4p, col_Op, col_O2p, $
  MIN_WAV=min_wav, MAX_WAV=max_wav, WAVELENGTHS=wavelengths
;+
; NAME:
;   IPT_CALCULATE_EMISSION_DOUBLE
;
; PURPOSE:
;   Calculate IPT emission line brightnesses for a double Maxwellian electron
;   distribution using pre-calculated 4D emission tables. This properly accounts
;   for the nonlinear enhancement of high-excitation transitions by suprathermal
;   electrons.
;
; CATEGORY:
;   Io Plasma Torus, Emission Modeling, Spectral Analysis
;
; CALLING SEQUENCE:
;   brightnesses = IPT_CALCULATE_EMISSION_DOUBLE(tables, core_temp, hot_temp, $
;                    total_density, hot_fraction, col_Sp, col_S2p, col_S3p, $
;                    col_S4p, col_Op, col_O2p, MIN_WAV=min_wav, MAX_WAV=max_wav, $
;                    WAVELENGTHS=wavelengths)
;
; INPUTS:
;   tables         - Structure from IPT_LOAD_DOUBLE_MAXWELLIAN_TABLES
;   core_temp      - Core (cold) electron temperature [eV]
;   hot_temp       - Hot electron temperature [eV]
;   total_density  - Total electron density [cm^-3]
;   hot_fraction   - Fraction of hot electrons (0 < feh < 1)
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
;   The double Maxwellian distribution represents:
;     f(v) = (1 - feh) * f_Maxwell(v, Tec) + feh * f_Maxwell(v, Teh)
;   
;   This arises from:
;   - Wave-particle interactions (electron cyclotron waves)
;   - Pickup ion acceleration and subsequent electron heating
;   - Magnetic reconnection events
;   - Centrifugally driven instabilities
;   
;   Effects on emission:
;   - High-excitation lines (e.g., S IV, S V) are strongly enhanced
;   - Enhancement is nonlinear (not proportional to feh)
;   - Ratio of high to low excitation lines increases
;   - Overall brightness can increase significantly
;
; TYPICAL PARAMETERS:
;   For Io Plasma Torus:
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
; EXAMPLE:
;   brightnesses = IPT_CALCULATE_EMISSION_DOUBLE(tables, 5.0d, 270.0d, $
;                    2205.5d, 0.0025d, 1.2d13, 4.2d13, 5.92d12, 6.0d11, $
;                    5.2d13, 5.92d12, WAVELENGTHS=xwavi)
;
; MODIFICATION HISTORY:
;   Written by E. G. Nerney, November 2025
;-

  ; Set default wavelength range
  IF ~KEYWORD_SET(max_wav) THEN max_wav = 2100.0d
  IF ~KEYWORD_SET(min_wav) THEN min_wav = 550.0d
  
  PRINT, 'Calculating double Maxwellian emission using 4D tables...'
  PRINT, '  Tec = ', STRTRIM(STRING(core_temp, FORMAT='(F6.2)'),2), ' eV'
  PRINT, '  Teh = ', STRTRIM(STRING(hot_temp, FORMAT='(F6.1)'),2), ' eV'
  PRINT, '  ne_total = ', STRTRIM(STRING(total_density, FORMAT='(F8.1)'),2), ' cm^-3'
  PRINT, '  feh = ', STRTRIM(STRING(hot_fraction, FORMAT='(F8.6)'),2)
  
  ; ============================================================================
  ; INTERPOLATE 4D EMISSIVITIES FOR EACH ION SPECIES
  ; ============================================================================
  ; Each species is interpolated separately in 4D space
  
  ; S+ (S II)
  emiss_sp = IPT_INTERPOLATE_EMISSIVITY_4D_SPECIES(core_temp, hot_temp, $
    total_density, hot_fraction, tables.tec_arr, tables.teh_arr, $
    tables.ne_arr, tables.feh_arr, tables.yptsi.SP)
  
  ; S++ (S III)
  emiss_s2p = IPT_INTERPOLATE_EMISSIVITY_4D_SPECIES(core_temp, hot_temp, $
    total_density, hot_fraction, tables.tec_arr, tables.teh_arr, $
    tables.ne_arr, tables.feh_arr, tables.yptsi.S2P)
  
  ; S+++ (S IV) - enhanced by hot electrons
  emiss_s3p = IPT_INTERPOLATE_EMISSIVITY_4D_SPECIES(core_temp, hot_temp, $
    total_density, hot_fraction, tables.tec_arr, tables.teh_arr, $
    tables.ne_arr, tables.feh_arr, tables.yptsi.S3P)
  
  ; S++++ (S V) - strongly enhanced by hot electrons
  emiss_s4p = IPT_INTERPOLATE_EMISSIVITY_4D_SPECIES(core_temp, hot_temp, $
    total_density, hot_fraction, tables.tec_arr, tables.teh_arr, $
    tables.ne_arr, tables.feh_arr, tables.yptsi.S4P)
  
  ; O+ (O II)
  emiss_op = IPT_INTERPOLATE_EMISSIVITY_4D_SPECIES(core_temp, hot_temp, $
    total_density, hot_fraction, tables.tec_arr, tables.teh_arr, $
    tables.ne_arr, tables.feh_arr, tables.yptsi.OP)
  
  ; O++ (O III)
  emiss_o2p = IPT_INTERPOLATE_EMISSIVITY_4D_SPECIES(core_temp, hot_temp, $
    total_density, hot_fraction, tables.tec_arr, tables.teh_arr, $
    tables.ne_arr, tables.feh_arr, tables.yptsi.O2P)
  
  ; ============================================================================
  ; CONVERT TO OBSERVABLE BRIGHTNESSES
  ; ============================================================================
  sp_brightness = (1d-6) * emiss_sp * col_Sp
  s2p_brightness = (1d-6) * emiss_s2p * col_S2p
  s3p_brightness = (1d-6) * emiss_s3p * col_S3p
  s4p_brightness = (1d-6) * emiss_s4p * col_S4p
  op_brightness = (1d-6) * emiss_op * col_Op
  o2p_brightness = (1d-6) * emiss_o2p * col_O2p
  
  ; ============================================================================
  ; COMBINE AND SORT
  ; ============================================================================
  xwavi_all = [tables.xwavi.SP, tables.xwavi.S2P, tables.xwavi.S3P, $
               tables.xwavi.S4P, tables.xwavi.OP, tables.xwavi.O2P]
  
  yptsi_all = [sp_brightness, s2p_brightness, s3p_brightness, $
               s4p_brightness, op_brightness, o2p_brightness]
  
  ; Sort by wavelength
  wsort = SORT(xwavi_all)
  xwavi_all = xwavi_all[wsort]
  yptsi_all = yptsi_all[wsort]
  
  ; Filter to wavelength range
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
;                       Shape: (n_wavelengths,)
;   bin_width         - Width of wavelength bins [Angstroms] (scalar)
;   line_wavelengths  - Wavelengths of discrete emission lines [Angstroms]
;                       Shape: (n_lines,)
;   line_brightnesses - Integrated brightnesses of emission lines [Rayleighs]
;                       Shape: (n_lines,)
;
; KEYWORDS:
;   FWHM - Full Width at Half Maximum of Gaussian instrument response [Angstroms]
;          Default = 6.0 (appropriate for Europa-UVS/JUICE-UVS)
;          Scalar value applied to all emission lines
;
; OUTPUTS:
;   Returns 1D array of spectral brightness per unit wavelength at each 
;   wavelength grid point. Shape: (n_wavelengths,)
;   Units: [Rayleighs/Angstrom]
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
;     lambda_upper = bin_center + bin_width/2 (upper bin edge)
;     lambda_lower = bin_center - bin_width/2 (lower bin edge)
;     sigma = FWHM / (2 * sqrt(2 * ln(2)))  [Gaussian width parameter]
;   
;   The error function relationship arises from:
;     integral[exp(-x^2)] = sqrt(pi) * erf(x)
;   
;   This method exactly integrates the Gaussian profile rather than sampling 
;   at the bin center, providing superior accuracy especially when:
;   - Bins are wide compared to line width (FWHM)
;   - Lines fall near bin edges
;   - High precision is required
;
; GAUSSIAN LINE PROFILE:
;   The normalized Gaussian profile is:
;     G(lambda) = 1/(sigma*sqrt(2*pi)) * exp(-(lambda-lambda_0)^2 / (2*sigma^2))
;   
;   Key parameters:
;   - FWHM = 2*sqrt(2*ln(2)) * sigma ≈ 2.355 * sigma
;   - sigma = FWHM / 2.355
;   - Total integrated flux = 1 (for normalized profile)
;   
;   Relationship to velocity width:
;   - Delta_v = c * FWHM / lambda_0
;   - For lambda = 1200 Angstrom, FWHM = 6 Angstrom -> Delta_v ≈ 1500 km/s
;
; INSTRUMENT RESPONSE:
;   Typical FWHM values for UV/optical spectrographs:
;   
;   Space-based UV instruments:
;   - HST/STIS G140M: ~0.6 Angstrom (~40 km/s at 1200 Angstrom)
;   - HST/STIS G230M: ~1.0 Angstrom
;   - HST/STIS low-res: ~6 Angstrom
;   - Europa-UVS: ~6 Angstrom (FWHM varies with wavelength)
;   - JUICE-UVS: ~6 Angstrom
;   - MAVEN/IUVS: ~0.6-1.2 Angstrom
;   
;   Ground-based optical instruments:
;   - High-resolution echelle (R~50000): ~0.1-0.2 Angstrom
;   - Medium-resolution (R~5000): ~1-2 Angstrom
;   - Low-resolution (R~1000): ~5-10 Angstrom
;   
;   The FWHM typically varies with wavelength due to optical aberrations,
;   detector pixel size, and grating dispersion.
;
; PHYSICAL INTERPRETATION:
;   The observed spectrum is the convolution of intrinsic line emission with
;   the instrument response:
;     I_obs(lambda) = integral[I_intrinsic(lambda') * R(lambda - lambda') dlambda']
;   
;   For discrete lines and Gaussian response:
;     I_obs(lambda) = sum_i[Brightness_i * Gaussian(lambda - lambda_i)]
;   
;   Bin integration accounts for finite detector pixel width, preventing
;   aliasing and providing photometrically accurate results.
;
; EXAMPLE:
;   ; Create wavelength grid from 800-1200 Angstroms with 1 Angstrom bins
;   wavelength = 800.0d + DINDGEN(401)
;   
;   ; Define emission lines (O III, S III, S IV multiplets)
;   line_wl = [832.9d, 906.9d, 1031.9d, 1062.7d, 1073.0d, 1099.1d]
;   line_br = [50.0d, 30.0d, 100.0d, 80.0d, 40.0d, 60.0d]  ; Rayleighs
;   
;   ; Simulate spectrum with 6 Angstrom resolution
;   spectrum = IPT_SIMULATE_SPECTRUM_ERF(wavelength, 1.0d, line_wl, line_br, $
;                                         FWHM=6.0d)
;   
;   ; Plot result
;   p = PLOT(wavelength, spectrum, XTITLE='Wavelength [$\Angstrom$]', $
;            YTITLE='Brightness [R/$\Angstrom$]', THICK=2)
;
; PERFORMANCE:
;   Fully vectorized implementation processes all emission lines and wavelength
;   bins simultaneously. Typical execution times:
;   - 100 lines, 1000 wavelength points: ~1 ms
;   - 1000 lines, 5000 wavelength points: ~50 ms
;   
;   Memory usage: O(n_lines * n_wavelengths) for temporary arrays
;
; NOTES:
;   - Lines with zero or negative brightness are automatically skipped
;   - Output units are brightness per unit wavelength [R/Angstrom], suitable
;     for comparison with spectrograph data
;   - For photon counting instruments, multiply by bin_width to get counts per bin
;   - This function assumes the same FWHM for all lines; wavelength-dependent
;     FWHM can be implemented by calling this function multiple times for
;     different wavelength regions
;
; MODIFICATION HISTORY:
;   Written by E. G. Nerney, November 2025
;-
  COMPILE_OPT IDL2, HIDDEN
  
  ; Set default FWHM if not specified
  ; Default value appropriate for Europa-UVS and JUICE-UVS instruments
  IF ~KEYWORD_SET(fwhm) THEN fwhm = 6.0d
  
  ; Convert FWHM to the error function scaling parameter
  ; For a Gaussian: FWHM = 2*sqrt(2*ln(2)) * sigma
  ; We need: 1/(sigma*sqrt(2)) for the ERF argument
  ; Therefore: 1/(sigma*sqrt(2)) = 2*sqrt(ln(2))/FWHM
  ; This parameter converts wavelength differences to ERF arguments
  rootc = 2.0d * SQRT(ALOG(2.0d)) / fwhm
  
  ; Get dimensions
  n_bins = N_ELEMENTS(wavelength_grid)      ; Number of output wavelength bins
  n_lines = N_ELEMENTS(line_brightnesses)   ; Number of discrete emission lines
  
  ; Find lines with positive brightness (skip zero or negative lines)
  good_idx = WHERE(line_brightnesses GT 0.0d, n_good)
  
  ; If no lines have positive brightness, return zero spectrum
  IF n_good EQ 0 THEN RETURN, DBLARR(n_bins)
  
  ; Extract only lines with positive brightness for efficiency
  active_wavelengths = line_wavelengths[good_idx]
  active_brightnesses = line_brightnesses[good_idx]
  
  ; Create 2D arrays for fully vectorized computation
  ; Dimensions: (n_good_lines, n_bins)
  ; REBIN replicates arrays efficiently without explicit loops
  
  ; Replicate emission line wavelengths across all bins
  ; Shape: (n_good, n_bins)
  line_wl_2d = REBIN(active_wavelengths, n_good, n_bins)
  
  ; Replicate wavelength grid across all lines
  ; Shape: (n_good, n_bins)
  grid_2d = REBIN(REFORM(wavelength_grid, 1, n_bins), n_good, n_bins)
  
  ; Replicate emission line brightnesses across all bins
  ; Shape: (n_good, n_bins)
  bright_2d = REBIN(active_brightnesses, n_good, n_bins)
  
  ; Calculate wavelength difference from each line to each bin center
  ; Shape: (n_good, n_bins)
  delta_lambda = grid_2d - line_wl_2d
  
  ; Calculate ERF arguments for upper and lower bin edges
  ; Upper edge: bin_center + bin_width/2
  ; Lower edge: bin_center - bin_width/2
  ; The factor rootc = 1/(sigma*sqrt(2)) converts to ERF argument
  half_bin = bin_width / 2.0d
  z_upper = (delta_lambda + half_bin) * rootc
  z_lower = (delta_lambda - half_bin) * rootc
  
  ; Compute the integrated Gaussian contribution over each bin
  ; Factor of 0.5 comes from the normalization of the error function
  ; erf(infinity) = 1, so integral from -inf to x is (1 + erf(x))/2
  ; For integral from a to b: (erf(b) - erf(a))/2
  ; Shape: (n_good, n_bins)
  contributions = bright_2d * 0.5d * (ERF(z_upper) - ERF(z_lower))
  
  ; Sum contributions from all emission lines for each wavelength bin
  ; TOTAL with dimension 1 sums over the first dimension (lines)
  ; Output shape: (n_bins,)
  spectrum = TOTAL(contributions, 1)
  
  ; Convert from integrated brightness per bin to brightness per unit wavelength
  ; This normalization gives units of [Rayleighs/Angstrom]
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
  ;   smooth functions, significantly more accurate than trapezoidal rule.
  ;
  ; CATEGORY:
  ;   Numerical Methods, Integration
  ;
  ; CALLING SEQUENCE:
  ;   integral = IPT_SIMPSON_INTEGRATE(x, y)
  ;
  ; INPUTS:
  ;   x - Array of x values (independent variable), shape (n,)
  ;       Must be monotonically increasing
  ;   y - Array of y values (dependent variable), shape (n,)
  ;       Function values at each x position
  ;
  ; OUTPUTS:
  ;   Returns the approximate integral (scalar) of y(x) dx over the range
  ;   [x[0], x[n-1]].
  ;
  ; METHOD:
  ;   Uses composite Simpson's rule with adaptive formulation for irregularly
  ;   spaced data. Simpson's rule approximates the integrand as a parabola
  ;   over each pair of adjacent intervals, providing O(h^4) accuracy for
  ;   smooth functions (where h is the grid spacing).
  ;
  ;   For each pair of consecutive intervals [x[i], x[i+1], x[i+2]]:
  ;     h1 = x[i+1] - x[i]
  ;     h2 = x[i+2] - x[i+1]
  ;     h = h1 + h2
  ;
  ;     Integral ≈ (h/6) * [α*y[i] + β*y[i+1] + γ*y[i+2]]
  ;
  ;     where:
  ;       α = (2*h1 - h2)/h1
  ;       β = (h1 + h2)^2/(h1*h2)
  ;       γ = (2*h2 - h1)/h2
  ;
  ;   For regularly spaced points (h1 = h2 = h), this reduces to the classical
  ;   Simpson's 1/3 rule:
  ;     Integral ≈ (h/3) * [y[i] + 4*y[i+1] + y[i+2]]
  ;
  ;   If the number of intervals (n-1) is odd, the final interval is integrated
  ;   using the trapezoidal rule to maintain consistency.
  ;
  ; ALGORITHM DETAILS:
  ;   Composite Simpson's rule divides the integration range into pairs of
  ;   adjacent intervals. The method:
  ;
  ;   1. Fits a quadratic (parabola) through each triplet of points
  ;   2. Integrates the quadratic exactly
  ;   3. Sums contributions from all interval pairs
  ;
  ;   Accuracy comparison for smooth functions:
  ;   - Trapezoidal rule: O(h^2) error
  ;   - Simpson's rule: O(h^4) error  [100x better for h=0.1]
  ;   - Romberg integration: O(h^6) error
  ;
  ;   For non-smooth functions (discontinuities, sharp peaks), Simpson's rule
  ;   may not provide significant improvement over trapezoidal rule.
  ;
  ; EDGE CASES:
  ;   - n < 2: Returns 0.0 (undefined integral)
  ;   - n = 2: Uses trapezoidal rule (only option with 2 points)
  ;   - n = 3: Uses single Simpson's rule application (exact for quadratics)
  ;   - n even: Uses Simpson's rule for (n-2) points, trapezoidal for last interval
  ;   - n odd: Uses Simpson's rule for all (n-1) intervals
  ;
  ; EXAMPLES:
  ;   ; Example 1: Integrate a smooth spectrum over a wavelength range
  ;   idx = WHERE(wavelength GE 680.0d AND wavelength LE 690.0d, count)
  ;   IF count GT 1 THEN BEGIN
  ;     brightness = IPT_SIMPSON_INTEGRATE(wavelength[idx], spectrum[idx])
  ;     PRINT, 'Total brightness: ', brightness, ' Rayleighs'
  ;   ENDIF
  ;
  ;   ; Example 2: Compute column density from density profile
  ;   ; n(z) = density as function of height z
  ;   column_density = IPT_SIMPSON_INTEGRATE(z_grid, density_profile)
  ;
  ;   ; Example 3: Compare with analytical result (parabola is exact)
  ;   x_test = [0.0d, 1.0d, 2.0d, 3.0d, 4.0d]
  ;   y_test = x_test^2  ; y = x^2
  ;   result = IPT_SIMPSON_INTEGRATE(x_test, y_test)
  ;   analytical = 4.0d^3 / 3.0d  ; integral of x^2 from 0 to 4
  ;   PRINT, 'Numerical: ', result, '  Analytical: ', analytical
  ;   ; Result: Both give 21.333... (exact to machine precision)
  ;
  ; PERFORMANCE:
  ;   Fully vectorized implementation processes all interval pairs simultaneously.
  ;   Typical execution times:
  ;   - 100 points: ~0.01 ms
  ;   - 10,000 points: ~0.5 ms
  ;   - 1,000,000 points: ~50 ms
  ;
  ;   Approximately 5-10x faster than loop-based implementation.
  ;
  ; NOTES:
  ;   - Input arrays must have the same length
  ;   - x values should be monotonically increasing for physical interpretation
  ;   - For closed-loop integrals or periodic functions, consider adding the
  ;     first point again at the end to close the loop
  ;   - Simpson's rule is most accurate when the integrand is smooth (continuous
  ;     derivatives up to 4th order)
  ;
  ; SEE ALSO:
  ;   IDL built-in INT_TABULATED (uses trapezoidal rule by default)
  ;
  ; MODIFICATION HISTORY:
  ;   Written by E. G. Nerney, November 2025
  ;-
  COMPILE_OPT IDL2, HIDDEN

  ; Get number of data points
  n = N_ELEMENTS(x)

  ; Handle edge cases
  IF n LT 2 THEN RETURN, 0.0d

  ; For two points, use trapezoidal rule (only option)
  IF n EQ 2 THEN RETURN, 0.5d * (y[0] + y[1]) * (x[1] - x[0])

  ; Determine number of complete Simpson's rule applications (pairs of intervals)
  ; Each application uses 3 points (2 intervals), with adjacent applications sharing
  ; the middle point. Therefore, n points give (n-1)/2 complete pairs if n is odd,
  ; or (n-2)/2 complete pairs with one leftover interval if n is even.
  n_pairs = (n - 1) / 2

  ; Extract point triplets for vectorized computation
  ; For n=7: indices 0,1,2,3,4,5,6
  ; Triplets: (0,1,2), (2,3,4), (4,5,6)
  ; Start indices: 0, 2, 4 -> 0:4:2 or 0:n-3:2
  indices = LINDGEN(n_pairs) * 2L

  ; Extract x and y values for all triplets simultaneously
  ; Each array has length n_pairs
  x0 = x[indices]          ; Left points
  x1 = x[indices + 1L]     ; Middle points
  x2 = x[indices + 2L]     ; Right points

  y0 = y[indices]
  y1 = y[indices + 1L]
  y2 = y[indices + 2L]

  ; Calculate interval widths for all pairs
  h1 = x1 - x0  ; Width of first interval in each pair
  h2 = x2 - x1  ; Width of second interval in each pair
  h = h1 + h2   ; Total width of the pair

  ; Apply adaptive Simpson's rule formula for irregular spacing
  ; This reduces to classical Simpson's 1/3 rule when h1 = h2
  ;
  ; Coefficient breakdown:
  ;   y0: (h/6) * (2 - h2/h1) = (h/6) * (2*h1 - h2)/h1
  ;   y1: (h/6) * h^2/(h1*h2)
  ;   y2: (h/6) * (2 - h1/h2) = (h/6) * (2*h2 - h1)/h2
  ;
  ; For h1 = h2 = h: coefficients become h/3, 4h/3, h/3 (classical formula)
  contributions = (h / 6.0d) * ((2.0d - h2/h1) * y0 + $
    (h * h / (h1 * h2)) * y1 + $
    (2.0d - h1/h2) * y2)

  ; Sum all Simpson's rule contributions
  integral = TOTAL(contributions, /DOUBLE)

  ; Handle remaining interval if number of intervals is odd
  ; n even => (n-1) intervals, and (n-1) odd => one leftover interval
  ; Example: n=6 has 5 intervals; Simpson's handles 4, leaving 1
  IF (n MOD 2) EQ 0 THEN BEGIN
    ; Add trapezoidal rule for the final interval [x[n-2], x[n-1]]
    integral += 0.5d * (y[n-2] + y[n-1]) * (x[n-1] - x[n-2])
  ENDIF

  RETURN, integral
END
;==============================================================================
PRO IPT_cm3_emission_model
;+
; NAME:
;   IPT_cm3_emission_model
;
; PURPOSE:
;   Dummy procedure to load the module. This allows the use of 
;   @IPT_cm3_emission_model syntax to compile all functions in this file.
;
; CATEGORY:
;   Io Plasma Torus, Module Loader
;
; CALLING SEQUENCE:
;   @IPT_cm3_emission_model
;
; DESCRIPTION:
;   This procedure prints a banner indicating that the IPT emission model
;   library has been loaded and lists available functions.
;
; MODIFICATION HISTORY:
;   Written by E. G. Nerney, November 2025
;-

  PRINT, ''
  PRINT, '=========================================================================='
  PRINT, 'IPT UV/Optical Emission Model Library Loaded'
  PRINT, 'MOP Community Code v1.0'
  PRINT, '=========================================================================='
  PRINT, ''
  PRINT, 'Available functions:'
  PRINT, '  IPT_LOAD_SINGLE_MAXWELLIAN_TABLES    - Load single Maxwellian tables'
  PRINT, '  IPT_LOAD_DOUBLE_MAXWELLIAN_TABLES    - Load double Maxwellian tables'
  PRINT, '  IPT_CALCULATE_EMISSION_SINGLE        - Calculate single Maxwellian emission'
  PRINT, '  IPT_CALCULATE_EMISSION_DOUBLE        - Calculate double Maxwellian emission'
  PRINT, '  IPT_SIMULATE_SPECTRUM_ERF            - Convolve with instrument response'
  PRINT, '  IPT_INTERPOLATE_EMISSIVITY_2D_SPECIES - 2D interpolation (single species)'
  PRINT, '  IPT_INTERPOLATE_EMISSIVITY_4D_SPECIES - 4D interpolation (single species)'
  PRINT, '  IPT_SIMPSON_INTEGRATE                - Numerical integration'
  PRINT, ''
  PRINT, 'For UV observations: use wavelength range 550-2100 Angstroms'
  PRINT, 'For optical observations: use wavelength range 3000-10000 Angstroms'
  PRINT, ''
  PRINT, 'CHIANTI atomic database v11.0.2'
  PRINT, 'Dere et al. 1997; Del Zanna et al. 2020; Dufresne et al. 2024'
  PRINT, ''
  PRINT, 'Author: E. G. Nerney, LASP, University of Colorado Boulder'
  PRINT, '=========================================================================='
  PRINT, ''

END