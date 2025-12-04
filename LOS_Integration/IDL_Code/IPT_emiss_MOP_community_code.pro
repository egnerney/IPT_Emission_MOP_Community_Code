;+
; IPT_emiss_MOP_community_code.pro
;
; Io Plasma Torus (IPT) UV/Optical Emission Raytracer Library
; ============================================================
;
; This library provides comprehensive tools for calculating emission spectra from
; the Io Plasma Torus using line-of-sight integration through a 3D plasma model
; with pre-calculated CHIANTI atomic database emission tables. Supports both single
; and double Maxwellian electron distributions with proper species-by-species
; interpolation and vectorized operations for optimal performance.
;
; PHYSICAL MODEL:
; - 3D plasma model on non-uniform rectilinear Cartesian grid
; - Trilinear interpolation of densities and temperatures along LOS
; - Bilinear interpolation of single Maxwellian photon emission rates (2D tables)
; - Quadlinear interpolation of double Maxwellian photon emission rates (4D tables)
; - Per-ion photon emission rate calculation (photons/s/ion)
; - Simpson's rule integration over line of sight for brightness (Rayleighs)
; - Gaussian instrumental response convolution via analytic ERF formulation
;
; LINE-OF-SIGHT INTEGRATION METHOD:
; - Ray tracing through 3D plasma model with configurable step size
; - Trilinear interpolation of plasma parameters at each LOS point
; - Bilinear/quadlinear interpolation of photon emission rates from tables
; - Per-ion emission rate: emissivity [photons/s/cm^3] / n_e [cm^-3] = [photons/s/ion]
; - Brightness integrand: (emission rate per ion) * n_ion [cm^-3] = [photons/s/cm^3]
; - Simpson's rule integration with R_J to cm conversion and Rayleigh factor
; - Final brightness in Rayleighs: integral * (R_J in cm) * 1e-6
;
; APPLICABLE WAVELENGTH RANGES:
; - UV instruments: 550-2100 Angstroms (Europa-UVS, JUICE-UVS, HST/STIS)
; - Optical instruments: 3000-10000 Angstroms (ground-based telescopes)
;
; COORDINATE SYSTEMS AND UNITS:
; - Positions: Jupiter radii [R_J], System III right-handed Cartesian
; - Temperatures: electron volts [eV]
; - Densities: particles per cubic centimeter [cm^-3]
; - Wavelengths: Angstroms [Angstroms]
; - Photon emission rates: photons per second per ion [photons s^-1 ion^-1]
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
; VERSION: 1.0
; DATE: November 2025
;
; CHIANTI ACKNOWLEDGMENT:
; CHIANTI is a collaborative project involving George Mason University, the University
; of Michigan (USA), University of Cambridge (UK), and NASA Goddard Space Flight Center (USA).
;-

;==============================================================================
FUNCTION IPT_SIMPSON_INTEGRATE, x, y
  ;+
  ; NAME:
  ;   IPT_SIMPSON_INTEGRATE
  ;
  ; PURPOSE:
  ;   Perform numerical integration using composite Simpson's rule for
  ;   irregularly spaced data points. Provides O(h^4) accuracy for smooth
  ;   functions, significantly more accurate than trapezoidal rule.
  ;
  ; CATEGORY:
  ;   Numerical Methods, Integration
  ;
  ; CALLING SEQUENCE:
  ;   integral = IPT_SIMPSON_INTEGRATE(x, y)
  ;
  ; INPUTS:
  ;   x - Array of x values (independent variable), must be monotonically increasing
  ;   y - Array of y values (dependent variable), function values at each x position
  ;
  ; OUTPUTS:
  ;   Returns the approximate integral (scalar) of y(x) dx over [x[0], x[n-1]]
  ;
  ; METHOD:
  ;   Uses composite Simpson's rule with adaptive formulation for irregular spacing.
  ;   For each pair of consecutive intervals, fits a parabola through three points
  ;   and integrates exactly. Handles both odd and even numbers of points.
  ;
  ; PERFORMANCE:
  ;   Fully vectorized implementation processes all interval pairs simultaneously.
  ;-
  COMPILE_OPT IDL2, HIDDEN

  n = N_ELEMENTS(x)

  ; Handle edge cases
  IF n LT 2 THEN RETURN, 0.0d
  IF n EQ 2 THEN RETURN, 0.5d * (y[0] + y[1]) * (x[1] - x[0])

  ; Number of complete Simpson's rule applications (pairs of intervals)
  n_pairs = (n - 1) / 2

  ; Extract point triplets for vectorized computation
  indices = LINDGEN(n_pairs) * 2L

  ; Extract x and y values for all triplets simultaneously
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

  ; Handle remaining interval if number of intervals is odd (n even)
  IF (n MOD 2) EQ 0 THEN BEGIN
    integral += 0.5d * (y[n-2] + y[n-1]) * (x[n-1] - x[n-2])
  ENDIF

  RETURN, integral
END

;==============================================================================
FUNCTION IPT_INTERPOLATE_EMISSIVITY_2D_VECTORIZED, temp_eV, dens_cm3, $
  temp_arr, dens_arr, emiss_table_species
  ;+
  ; NAME:
  ;   IPT_INTERPOLATE_EMISSIVITY_2D_VECTORIZED
  ;
  ; PURPOSE:
  ;   Interpolate 2D emissivity rates for a single species from single
  ;   Maxwellian emission tables at a specific temperature and density.
  ;   Performs bilinear interpolation in log10(Te)-log10(ne) space.
  ;   Fully vectorized to interpolate all emission lines simultaneously.
  ;
  ; CATEGORY:
  ;   Io Plasma Torus, Emission Modeling, Atomic Physics
  ;
  ; CALLING SEQUENCE:
  ;   emiss_interp = IPT_INTERPOLATE_EMISSIVITY_2D_VECTORIZED(temp_eV, dens_cm3, $
  ;                    temp_arr, dens_arr, emiss_table_species)
  ;
  ; INPUTS:
  ;   temp_eV             - Electron temperature [eV] to interpolate at (scalar)
  ;   dens_cm3            - Electron density [cm^-3] to interpolate at (scalar)
  ;   temp_arr            - Temperature grid from tables [eV], shape (n_temp,)
  ;   dens_arr            - Density grid from tables [cm^-3], shape (n_dens,)
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
  ;   Bilinear interpolation in log10(Te)-log10(ne) space, fully vectorized
  ;   to process all emission lines simultaneously.
  ;
  ; PERFORMANCE:
  ;   Approximately 10-50x faster than looping over emission lines.
  ;-
  COMPILE_OPT IDL2, HIDDEN

  ; Get dimensions of the emission table
  dims = SIZE(emiss_table_species, /DIMENSIONS)
  n_temp = dims[0]
  n_dens = dims[1]
  n_lines = dims[2]

  ; Check for valid input
  IF temp_eV LE 0 OR dens_cm3 LE 0 THEN RETURN, DBLARR(n_lines)

  ; Convert to log10 space for interpolation
  log10_temp = ALOG10(temp_eV)
  log10_dens = ALOG10(dens_cm3)
  log10_temp_arr = ALOG10(temp_arr)
  log10_dens_arr = ALOG10(dens_arr)

  ; Check bounds - return zeros if outside table range
  IF log10_temp LT log10_temp_arr[0] OR log10_temp GT log10_temp_arr[n_temp-1] THEN $
    RETURN, DBLARR(n_lines)
  IF log10_dens LT log10_dens_arr[0] OR log10_dens GT log10_dens_arr[n_dens-1] THEN $
    RETURN, DBLARR(n_lines)

  ; Find bracketing indices using VALUE_LOCATE
  i_temp = VALUE_LOCATE(log10_temp_arr, log10_temp)
  i_dens = VALUE_LOCATE(log10_dens_arr, log10_dens)

  ; Clamp indices to valid interpolation range [0, n-2]
  i_temp = 0 > i_temp < (n_temp - 2)
  i_dens = 0 > i_dens < (n_dens - 2)

  ; Calculate fractional position between grid points
  frac_temp = (log10_temp - log10_temp_arr[i_temp]) / $
    (log10_temp_arr[i_temp + 1] - log10_temp_arr[i_temp])
  frac_dens = (log10_dens - log10_dens_arr[i_dens]) / $
    (log10_dens_arr[i_dens + 1] - log10_dens_arr[i_dens])

  ; Clamp fractional positions to [0, 1]
  frac_temp = 0.0d > frac_temp < 1.0d
  frac_dens = 0.0d > frac_dens < 1.0d

  ; Extract 4 corner values for all emission lines
  v00 = REFORM(emiss_table_species[i_temp,   i_dens,   *], n_lines)
  v01 = REFORM(emiss_table_species[i_temp,   i_dens+1, *], n_lines)
  v10 = REFORM(emiss_table_species[i_temp+1, i_dens,   *], n_lines)
  v11 = REFORM(emiss_table_species[i_temp+1, i_dens+1, *], n_lines)

  ; Vectorized bilinear interpolation
  emiss_interp = v00 * (1.0d - frac_temp) * (1.0d - frac_dens) + $
    v01 * (1.0d - frac_temp) * frac_dens + $
    v10 * frac_temp * (1.0d - frac_dens) + $
    v11 * frac_temp * frac_dens

  ; Set any negative values to zero
  neg_idx = WHERE(emiss_interp LT 0.0d, n_neg)
  IF n_neg GT 0 THEN emiss_interp[neg_idx] = 0.0d

  RETURN, emiss_interp
END

;==============================================================================
FUNCTION IPT_INTERPOLATE_EMISSIVITY_4D_VECTORIZED, Tec_eV, Teh_eV, ne_cm3, feh, $
  tec_arr, teh_arr, ne_arr, feh_arr, emiss_table_species
  ;+
  ; NAME:
  ;   IPT_INTERPOLATE_EMISSIVITY_4D_VECTORIZED
  ;
  ; PURPOSE:
  ;   Interpolate 4D emissivity rates for a single species from double
  ;   Maxwellian emission tables. Performs quadralinear (4D linear)
  ;   interpolation in log10 parameter space. Fully vectorized using
  ;   matrix operations for all emission lines simultaneously.
  ;
  ; CATEGORY:
  ;   Io Plasma Torus, Emission Modeling, Atomic Physics
  ;
  ; CALLING SEQUENCE:
  ;   emiss_interp = IPT_INTERPOLATE_EMISSIVITY_4D_VECTORIZED(Tec_eV, Teh_eV, $
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
  ;-
  COMPILE_OPT IDL2, HIDDEN

  ; Get dimensions
  dims = SIZE(emiss_table_species, /DIMENSIONS)
  n_tec = dims[0]
  n_teh = dims[1]
  n_ne = dims[2]
  n_feh = dims[3]
  n_lines = dims[4]

  ; Check for valid input
  IF Tec_eV LE 0 OR Teh_eV LE 0 OR ne_cm3 LE 0 OR feh LE 0 THEN RETURN, DBLARR(n_lines)

  ; Convert to log10 space (feh stays linear or log depending on table)
  log10_tec = ALOG10(Tec_eV)
  log10_teh = ALOG10(Teh_eV)
  log10_ne = ALOG10(ne_cm3)
  log10_feh = ALOG10(feh)

  log10_tec_arr = ALOG10(tec_arr)
  log10_teh_arr = ALOG10(teh_arr)
  log10_ne_arr = ALOG10(ne_arr)
  log10_feh_arr = ALOG10(feh_arr)

  ; Check bounds - return zeros if outside table range
  IF log10_tec LT log10_tec_arr[0] OR log10_tec GT log10_tec_arr[n_tec-1] THEN $
    RETURN, DBLARR(n_lines)
  IF log10_teh LT log10_teh_arr[0] OR log10_teh GT log10_teh_arr[n_teh-1] THEN $
    RETURN, DBLARR(n_lines)
  IF log10_ne LT log10_ne_arr[0] OR log10_ne GT log10_ne_arr[n_ne-1] THEN $
    RETURN, DBLARR(n_lines)
  IF log10_feh LT log10_feh_arr[0] OR log10_feh GT log10_feh_arr[n_feh-1] THEN $
    RETURN, DBLARR(n_lines)

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
  emiss_interp = REFORM(corner_weights # corners, n_lines)

  ; Clamp negative values to zero
  neg_idx = WHERE(emiss_interp LT 0.0d, n_neg)
  IF n_neg GT 0 THEN emiss_interp[neg_idx] = 0.0d

  RETURN, emiss_interp
END
;==============================================================================
FUNCTION IPT_TRILINEAR_INTERPOLATE_VECTORIZED, positions_x, positions_y, positions_z, $
  x_axis, y_axis, z_axis, field_3d
  ;+
  ; NAME:
  ;   IPT_TRILINEAR_INTERPOLATE_VECTORIZED
  ;
  ; PURPOSE:
  ;   Perform trilinear interpolation on a 3D field at multiple positions
  ;   along a line of sight. Fully vectorized for optimal performance.
  ;
  ; CATEGORY:
  ;   Numerical Methods, Interpolation
  ;
  ; CALLING SEQUENCE:
  ;   values = IPT_TRILINEAR_INTERPOLATE_VECTORIZED(positions_x, positions_y, $
  ;              positions_z, x_axis, y_axis, z_axis, field_3d)
  ;
  ; INPUTS:
  ;   positions_x, positions_y, positions_z - Arrays of positions [R_J]
  ;   x_axis, y_axis, z_axis - 1D coordinate arrays defining the grid
  ;   field_3d - 3D field array with shape [n_z, n_y, n_x]
  ;
  ; OUTPUTS:
  ;   Returns array of interpolated values at each position
  ;   Invalid positions (outside grid or NaN) return 0.0
  ;
  ; METHOD:
  ;   Uses vectorized index finding and weight calculation for all positions
  ;   simultaneously, then loops only for the final interpolation step.
  ;-
  COMPILE_OPT IDL2, HIDDEN

  n_points = N_ELEMENTS(positions_x)
  result = DBLARR(n_points)

  ; Get grid bounds
  x_min = x_axis[0]
  x_max = x_axis[N_ELEMENTS(x_axis)-1]
  y_min = y_axis[0]
  y_max = y_axis[N_ELEMENTS(y_axis)-1]
  z_min = z_axis[0]
  z_max = z_axis[N_ELEMENTS(z_axis)-1]

  n_x = N_ELEMENTS(x_axis)
  n_y = N_ELEMENTS(y_axis)
  n_z = N_ELEMENTS(z_axis)

  ; Find valid points (inside grid bounds)
  valid = (positions_x GE x_min) AND (positions_x LE x_max) AND $
    (positions_y GE y_min) AND (positions_y LE y_max) AND $
    (positions_z GE z_min) AND (positions_z LE z_max)

  valid_idx = WHERE(valid, n_valid)
  IF n_valid EQ 0 THEN RETURN, result

  ; Process valid points
  FOR ii = 0L, n_valid-1L DO BEGIN
    i = valid_idx[ii]
    x = positions_x[i]
    y = positions_y[i]
    z = positions_z[i]

    ; Find indices
    ix = VALUE_LOCATE(x_axis, x)
    iy = VALUE_LOCATE(y_axis, y)
    iz = VALUE_LOCATE(z_axis, z)

    ; Clamp indices
    ix = 0 > ix < (n_x - 2)
    iy = 0 > iy < (n_y - 2)
    iz = 0 > iz < (n_z - 2)

    ; Calculate weights
    wx = (x - x_axis[ix]) / (x_axis[ix+1] - x_axis[ix])
    wy = (y - y_axis[iy]) / (y_axis[iy+1] - y_axis[iy])
    wz = (z - z_axis[iz]) / (z_axis[iz+1] - z_axis[iz])

    ; Get 8 corner values (field_3d is [n_z, n_y, n_x])
    v000 = field_3d[iz,   iy,   ix]
    v001 = field_3d[iz,   iy,   ix+1]
    v010 = field_3d[iz,   iy+1, ix]
    v011 = field_3d[iz,   iy+1, ix+1]
    v100 = field_3d[iz+1, iy,   ix]
    v101 = field_3d[iz+1, iy,   ix+1]
    v110 = field_3d[iz+1, iy+1, ix]
    v111 = field_3d[iz+1, iy+1, ix+1]

    ; Handle NaN or negative corner values
    IF ~FINITE(v000) OR v000 LT 0 THEN v000 = 0.0d
    IF ~FINITE(v001) OR v001 LT 0 THEN v001 = 0.0d
    IF ~FINITE(v010) OR v010 LT 0 THEN v010 = 0.0d
    IF ~FINITE(v011) OR v011 LT 0 THEN v011 = 0.0d
    IF ~FINITE(v100) OR v100 LT 0 THEN v100 = 0.0d
    IF ~FINITE(v101) OR v101 LT 0 THEN v101 = 0.0d
    IF ~FINITE(v110) OR v110 LT 0 THEN v110 = 0.0d
    IF ~FINITE(v111) OR v111 LT 0 THEN v111 = 0.0d

    ; Trilinear interpolation
    interp_val = $
      v000 * (1-wx) * (1-wy) * (1-wz) + $
      v001 * wx * (1-wy) * (1-wz) + $
      v010 * (1-wx) * wy * (1-wz) + $
      v011 * wx * wy * (1-wz) + $
      v100 * (1-wx) * (1-wy) * wz + $
      v101 * wx * (1-wy) * wz + $
      v110 * (1-wx) * wy * wz + $
      v111 * wx * wy * wz

    IF FINITE(interp_val) AND interp_val GE 0 THEN result[i] = interp_val
  ENDFOR

  RETURN, result
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
  ;     Flux_bin = (Brightness_line / 2) * [erf(z_upper) - erf(z_lower)]
  ;
  ; PERFORMANCE:
  ;   Fully vectorized implementation using 2D array operations.
  ;-
  COMPILE_OPT IDL2, HIDDEN

  ; Set default FWHM
  IF ~KEYWORD_SET(fwhm) THEN fwhm = 6.0d

  ; Convert FWHM to ERF scaling parameter
  rootc = 2.0d * SQRT(ALOG(2.0d)) / fwhm

  ; Get dimensions
  n_bins = N_ELEMENTS(wavelength_grid)
  n_lines = N_ELEMENTS(line_brightnesses)

  ; Find lines with positive brightness
  good_idx = WHERE(line_brightnesses GT 0.0d, n_good)

  ; Return zero spectrum if no positive lines
  IF n_good EQ 0 THEN RETURN, DBLARR(n_bins)

  ; Extract active lines
  active_wavelengths = line_wavelengths[good_idx]
  active_brightnesses = line_brightnesses[good_idx]

  ; Create 2D arrays for vectorized computation: (n_good, n_bins)
  line_wl_2d = REBIN(active_wavelengths, n_good, n_bins)
  grid_2d = REBIN(REFORM(wavelength_grid, 1, n_bins), n_good, n_bins)
  bright_2d = REBIN(active_brightnesses, n_good, n_bins)

  ; Calculate wavelength difference
  delta_lambda = grid_2d - line_wl_2d

  ; Calculate ERF arguments for bin edges
  half_bin = bin_width / 2.0d
  z_upper = (delta_lambda + half_bin) * rootc
  z_lower = (delta_lambda - half_bin) * rootc

  ; Compute integrated Gaussian contribution
  contributions = bright_2d * 0.5d * (ERF(z_upper) - ERF(z_lower))

  ; Sum contributions and normalize
  spectrum = TOTAL(contributions, 1) / bin_width

  RETURN, spectrum
END

;==============================================================================
FUNCTION IPT_LOAD_PLASMA_MODEL, filename
  ;+
  ; NAME:
  ;   IPT_LOAD_PLASMA_MODEL
  ;
  ; PURPOSE:
  ;   Load 3D plasma model from HDF5 file containing densities, temperatures,
  ;   and hot fractions on a rectilinear Cartesian grid.
  ;
  ; CALLING SEQUENCE:
  ;   plasma_model = IPT_LOAD_PLASMA_MODEL(filename)
  ;
  ; INPUTS:
  ;   filename - Path to HDF5 file containing plasma model
  ;
  ; OUTPUTS:
  ;   plasma_model - Structure containing:
  ;     .x_axis, .y_axis, .z_axis - 1D coordinate arrays [R_J]
  ;     .nec  - Core electron density [cm^-3], 3D array [n_z, n_y, n_x]
  ;     .neh  - Hot electron density [cm^-3], 3D array [n_z, n_y, n_x]
  ;     .ne_total - Total electron density [cm^-3], 3D array [n_z, n_y, n_x]
  ;     .Tec  - Core electron temperature [eV], 3D array [n_z, n_y, n_x]
  ;     .Teh  - Hot electron temperature [eV], 3D array [n_z, n_y, n_x]
  ;     .feh  - Hot electron fraction, 3D array [n_z, n_y, n_x]
  ;     .nsp, .ns2p, .ns3p - Sulfur ion densities [cm^-3], 3D arrays
  ;     .nop, .no2p - Oxygen ion densities [cm^-3], 3D arrays
  ;-
  PRINT, 'Loading 3D plasma model from HDF5 file...'
  PRINT, '  File: ' + filename

  file_id = H5F_OPEN(filename)

  ; Read coordinate axes from coordinates group
  x_axis = H5D_READ(H5D_OPEN(file_id, 'coordinates/x'))
  y_axis = H5D_READ(H5D_OPEN(file_id, 'coordinates/y'))
  z_axis = H5D_READ(H5D_OPEN(file_id, 'coordinates/z'))

  PRINT, '  Grid shape: ', STRTRIM(STRING(N_ELEMENTS(x_axis)),2), ' x ', $
    STRTRIM(STRING(N_ELEMENTS(y_axis)),2), ' x ', STRTRIM(STRING(N_ELEMENTS(z_axis)),2)
  PRINT, '  X range: ', STRTRIM(STRING(MIN(x_axis), FORMAT='(F6.1)'),2), ' to ', $
    STRTRIM(STRING(MAX(x_axis), FORMAT='(F6.1)'),2), ' R_J'
  PRINT, '  Y range: ', STRTRIM(STRING(MIN(y_axis), FORMAT='(F6.1)'),2), ' to ', $
    STRTRIM(STRING(MAX(y_axis), FORMAT='(F6.1)'),2), ' R_J'
  PRINT, '  Z range: ', STRTRIM(STRING(MIN(z_axis), FORMAT='(F6.1)'),2), ' to ', $
    STRTRIM(STRING(MAX(z_axis), FORMAT='(F6.1)'),2), ' R_J'

  ; Read 3D fields from data group
  ; HDF5 Python stores as (n_x, n_y, n_z) C-order
  ; IDL reads as (n_z, n_y, n_x) due to row-major vs column-major difference
  ; We keep this ordering for trilinear interpolation

  ; Electron densities [cm^-3]
  nec = H5D_READ(H5D_OPEN(file_id, 'data/ne_c'))   ; Cold/core electron density
  neh = H5D_READ(H5D_OPEN(file_id, 'data/ne_h'))   ; Hot electron density

  ; Electron temperatures [eV]
  Tec = H5D_READ(H5D_OPEN(file_id, 'data/Te_c'))   ; Cold electron temperature
  Teh = H5D_READ(H5D_OPEN(file_id, 'data/Te_h'))   ; Hot electron temperature

  ; Ion densities [cm^-3]
  nsp = H5D_READ(H5D_OPEN(file_id, 'data/nsp'))    ; S+ density
  ns2p = H5D_READ(H5D_OPEN(file_id, 'data/ns2p'))  ; S++ density
  ns3p = H5D_READ(H5D_OPEN(file_id, 'data/ns3p'))  ; S+++ density
  nop = H5D_READ(H5D_OPEN(file_id, 'data/nop'))    ; O+ density
  no2p = H5D_READ(H5D_OPEN(file_id, 'data/no2p'))  ; O++ density

  H5F_CLOSE, file_id

  ; Calculate derived quantities
  ne_total = nec + neh
  feh = DBLARR(SIZE(ne_total, /DIMENSIONS))
  idx_valid = WHERE(ne_total GT 0, count_valid)
  IF count_valid GT 0 THEN feh[idx_valid] = neh[idx_valid] / ne_total[idx_valid]

  ; Handle NaN or invalid values in all fields
  all_fields = LIST(nec, neh, ne_total, Tec, Teh, feh, nsp, ns2p, ns3p, nop, no2p)
  FOR i = 0, N_ELEMENTS(all_fields)-1 DO BEGIN
    field = all_fields[i]
    idx_bad = WHERE(~FINITE(field) OR field LT 0, count_bad)
    IF count_bad GT 0 THEN field[idx_bad] = 0.0d
    all_fields[i] = field
  ENDFOR

  nec = all_fields[0]
  neh = all_fields[1]
  ne_total = all_fields[2]
  Tec = all_fields[3]
  Teh = all_fields[4]
  feh = all_fields[5]
  nsp = all_fields[6]
  ns2p = all_fields[7]
  ns3p = all_fields[8]
  nop = all_fields[9]
  no2p = all_fields[10]

  plasma_model = { $
    x_axis: x_axis, $
    y_axis: y_axis, $
    z_axis: z_axis, $
    nec: nec, $
    neh: neh, $
    ne_total: ne_total, $
    Tec: Tec, $
    Teh: Teh, $
    feh: feh, $
    nsp: nsp, $
    ns2p: ns2p, $
    ns3p: ns3p, $
    nop: nop, $
    no2p: no2p $
  }

  PRINT, 'Plasma model loaded successfully'
  RETURN, plasma_model
END

;==============================================================================
FUNCTION IPT_LOAD_EMISSION_TABLES_SINGLE, filename
  ;+
  ; NAME:
  ;   IPT_LOAD_EMISSION_TABLES_SINGLE
  ;
  ; PURPOSE:
  ;   Load single Maxwellian emission tables from HDF5 file and organize
  ;   by species for efficient interpolation.
  ;
  ; CALLING SEQUENCE:
  ;   tables_single = IPT_LOAD_EMISSION_TABLES_SINGLE(filename)
  ;
  ; INPUTS:
  ;   filename - Path to HDF5 file containing CHIANTI emission tables
  ;
  ; OUTPUTS:
  ;   tables_single - Structure containing emission tables by species
  ;-
  PRINT, 'Loading single Maxwellian emission tables...'
  PRINT, '  File: ' + filename

  file_id = H5F_OPEN(filename)

  ; Read parameter grids
  temp_arr = H5D_READ(H5D_OPEN(file_id, 'T'))
  dens_arr = H5D_READ(H5D_OPEN(file_id, 'n'))

  ; Read wavelengths and species
  wavelength_all = H5D_READ(H5D_OPEN(file_id, 'wavelength'))
  species_all = H5D_READ(H5D_OPEN(file_id, 'species'))

  ; Read volume emission rates and transpose
  emissivity_all = H5D_READ(H5D_OPEN(file_id, 'emiss'))
  emissivity_all = TRANSPOSE(emissivity_all, [2, 1, 0])

  H5F_CLOSE, file_id

  ; Organize by species
  species_names = ['S', 'SP', 'S2P', 'S3P', 'S4P', 'O', 'OP', 'O2P', 'NAP']

  wavelengths = CREATE_STRUCT('DUMMY', 0.0d)
  emissivities = CREATE_STRUCT('DUMMY', DBLARR(1,1,1))
  first_valid = 1

  PRINT, '  Organizing by species:'
  FOR i = 0, N_ELEMENTS(species_names)-1 DO BEGIN
    species_name = species_names[i]
    idx = WHERE(species_all EQ species_name, count)

    IF count GT 0 THEN BEGIN
      waves = wavelength_all[idx]
      emiss = emissivity_all[*, *, idx]

      field_name = STRUPCASE(species_name)
      IF first_valid EQ 1 THEN BEGIN
        wavelengths = CREATE_STRUCT(field_name, waves)
        emissivities = CREATE_STRUCT(field_name, emiss)
        first_valid = 0
      ENDIF ELSE BEGIN
        wavelengths = CREATE_STRUCT(wavelengths, field_name, waves)
        emissivities = CREATE_STRUCT(emissivities, field_name, emiss)
      ENDELSE

      PRINT, '    ', species_name, ': ', STRTRIM(STRING(count),2), ' lines'
    ENDIF
  ENDFOR

  tables_single = { $
    temp_arr: temp_arr, $
    dens_arr: dens_arr, $
    wavelengths: wavelengths, $
    emissivities: emissivities $
  }

  PRINT, 'Single Maxwellian tables loaded successfully'
  PRINT, '  Temperature range: ', STRTRIM(STRING(MIN(temp_arr), FORMAT='(F6.3)'),2), $
    ' - ', STRTRIM(STRING(MAX(temp_arr), FORMAT='(F6.1)'),2), ' eV'
  PRINT, '  Density range: ', STRTRIM(STRING(MIN(dens_arr), FORMAT='(F6.1)'),2), $
    ' - ', STRTRIM(STRING(MAX(dens_arr), FORMAT='(F8.1)'),2), ' cm^-3'

  RETURN, tables_single
END

;==============================================================================
FUNCTION IPT_LOAD_EMISSION_TABLES_DOUBLE, filename
  ;+
  ; NAME:
  ;   IPT_LOAD_EMISSION_TABLES_DOUBLE
  ;
  ; PURPOSE:
  ;   Load double Maxwellian emission tables from HDF5 file and organize
  ;   by species for efficient 4D interpolation.
  ;
  ; CALLING SEQUENCE:
  ;   tables_double = IPT_LOAD_EMISSION_TABLES_DOUBLE(filename)
  ;
  ; INPUTS:
  ;   filename - Path to HDF5 file containing CHIANTI double Maxwellian tables
  ;
  ; OUTPUTS:
  ;   tables_double - Structure containing emission tables by species
  ;-
  PRINT, 'Loading double Maxwellian emission tables...'
  PRINT, '  File: ' + filename

  file_id = H5F_OPEN(filename)

  ; Read parameter grids
  tec_arr = H5D_READ(H5D_OPEN(file_id, 'T_cold'))
  teh_arr = H5D_READ(H5D_OPEN(file_id, 'T_hot'))
  ne_arr = H5D_READ(H5D_OPEN(file_id, 'n'))
  feh_arr = H5D_READ(H5D_OPEN(file_id, 'feh'))

  ; Read wavelengths and species
  wavelength_all = H5D_READ(H5D_OPEN(file_id, 'wavelength'))
  species_all = H5D_READ(H5D_OPEN(file_id, 'species'))

  ; Read volume emission rates and transpose
  emissivity_all = H5D_READ(H5D_OPEN(file_id, 'emiss'))
  emissivity_all = TRANSPOSE(emissivity_all, [4, 3, 2, 1, 0])

  H5F_CLOSE, file_id

  ; Organize by species
  species_names = ['S', 'SP', 'S2P', 'S3P', 'S4P', 'O', 'OP', 'O2P', 'NAP']

  wavelengths = CREATE_STRUCT('DUMMY', 0.0d)
  emissivities = CREATE_STRUCT('DUMMY', DBLARR(1,1,1,1,1))
  first_valid = 1

  PRINT, '  Organizing by species:'
  FOR i = 0, N_ELEMENTS(species_names)-1 DO BEGIN
    species_name = species_names[i]
    idx = WHERE(species_all EQ species_name, count)

    IF count GT 0 THEN BEGIN
      waves = wavelength_all[idx]
      emiss = emissivity_all[*, *, *, *, idx]

      field_name = STRUPCASE(species_name)
      IF first_valid EQ 1 THEN BEGIN
        wavelengths = CREATE_STRUCT(field_name, waves)
        emissivities = CREATE_STRUCT(field_name, emiss)
        first_valid = 0
      ENDIF ELSE BEGIN
        wavelengths = CREATE_STRUCT(wavelengths, field_name, waves)
        emissivities = CREATE_STRUCT(emissivities, field_name, emiss)
      ENDELSE

      PRINT, '    ', species_name, ': ', STRTRIM(STRING(count),2), ' lines'
    ENDIF
  ENDFOR

  tables_double = { $
    tec_arr: tec_arr, $
    teh_arr: teh_arr, $
    ne_arr: ne_arr, $
    feh_arr: feh_arr, $
    wavelengths: wavelengths, $
    emissivities: emissivities $
  }

  PRINT, 'Double Maxwellian tables loaded successfully'
  PRINT, '  Core temp range: ', STRTRIM(STRING(MIN(tec_arr), FORMAT='(F6.3)'),2), $
    ' - ', STRTRIM(STRING(MAX(tec_arr), FORMAT='(F6.1)'),2), ' eV'
  PRINT, '  Hot temp range: ', STRTRIM(STRING(MIN(teh_arr), FORMAT='(F6.1)'),2), $
    ' - ', STRTRIM(STRING(MAX(teh_arr), FORMAT='(F7.0)'),2), ' eV'
  PRINT, '  Density range: ', STRTRIM(STRING(MIN(ne_arr), FORMAT='(F8.1)'),2), $
    ' - ', STRTRIM(STRING(MAX(ne_arr), FORMAT='(F10.1)'),2), ' cm^-3'
  PRINT, '  Hot fraction range: ', STRTRIM(STRING(MIN(feh_arr), FORMAT='(F8.6)'),2), $
    ' - ', STRTRIM(STRING(MAX(feh_arr), FORMAT='(F6.4)'),2)

  RETURN, tables_double
END

;==============================================================================
FUNCTION IPT_GET_ION_DENSITY_FIELD, plasma_model, species_key
  ;+
  ; NAME:
  ;   IPT_GET_ION_DENSITY_FIELD
  ;
  ; PURPOSE:
  ;   Return the appropriate ion density 3D field from the plasma model
  ;   based on species key.
  ;
  ; INPUTS:
  ;   plasma_model - Plasma model structure
  ;   species_key  - Species identifier ('SP', 'S2P', 'S3P', 'OP', 'O2P')
  ;
  ; OUTPUTS:
  ;   Returns the 3D ion density field [cm^-3]
  ;-
  COMPILE_OPT IDL2, HIDDEN

  CASE STRUPCASE(species_key) OF
    'SP': RETURN, plasma_model.nsp
    'S2P': RETURN, plasma_model.ns2p
    'S3P': RETURN, plasma_model.ns3p
    'S4P': RETURN, plasma_model.ns3p * 0.0d   ; S4P not typically in model
    'OP': RETURN, plasma_model.nop
    'O2P': RETURN, plasma_model.no2p
    ELSE: RETURN, plasma_model.nec * 0.0d
  ENDCASE
END

;==============================================================================
FUNCTION IPT_INTEGRATE_SPECIES_EMISSION_SINGLE, plasma_model, tables_single, $
  positions_x, positions_y, positions_z, s_values, species_key
  ;+
  ; NAME:
  ;   IPT_INTEGRATE_SPECIES_EMISSION_SINGLE
  ;
  ; PURPOSE:
  ;   Integrate all emission lines for a single species along line of sight
  ;   using single Maxwellian distribution. Returns brightness in Rayleighs
  ;   for each emission line of this species.
  ;
  ; CATEGORY:
  ;   Io Plasma Torus, Emission Modeling, Line-of-Sight Integration
  ;
  ; CALLING SEQUENCE:
  ;   brightnesses = IPT_INTEGRATE_SPECIES_EMISSION_SINGLE(plasma_model, $
  ;                    tables_single, positions_x, positions_y, positions_z, $
  ;                    s_values, species_key)
  ;
  ; INPUTS:
  ;   plasma_model  - Plasma model structure from IPT_LOAD_PLASMA_MODEL
  ;   tables_single - Single Maxwellian tables from IPT_LOAD_EMISSION_TABLES_SINGLE
  ;   positions_x/y/z - Arrays of positions along LOS [R_J]
  ;   s_values      - Path length values [R_J]
  ;   species_key   - Species name ('SP', 'S2P', 'S3P', 'OP', 'O2P')
  ;
  ; OUTPUTS:
  ;   Returns 1D array of line brightnesses in Rayleighs for all emission
  ;   lines of this species.
  ;
  ; METHOD:
  ;   1. Interpolate plasma parameters (Te, ne, n_ion) along LOS
  ;   2. For each valid LOS point, interpolate per-ion photon emission rate
  ;      Per-ion rate = volume_emissivity [photons/s/cm^3] / n_e [cm^-3]
  ;   3. Calculate brightness integrand: per_ion_rate * n_ion [photons/s/cm^3]
  ;   4. Integrate over LOS with Simpson's rule
  ;   5. Convert to Rayleighs: integral * R_J_cm * 1e-6
  ;-
  COMPILE_OPT IDL2, HIDDEN

  ; Physical constants
  R_J_CM = 7.1492d9  ; Jupiter radius in cm
  RAYLEIGH_FACTOR = 1.0d-6 * R_J_CM  ; Conversion to Rayleighs

  n_points = N_ELEMENTS(s_values)

  ; Get species tag index
  species_tags = TAG_NAMES(tables_single.wavelengths)
  species_idx = WHERE(species_tags EQ STRUPCASE(species_key), n_match)
  IF n_match EQ 0 THEN BEGIN
    PRINT, '  Warning: Species ', species_key, ' not found in tables'
    RETURN, [0.0d]
  ENDIF

  ; Get number of lines for this species
  n_lines = N_ELEMENTS(tables_single.wavelengths.(species_idx))

  ; Interpolate plasma parameters along LOS
  ne_los = IPT_TRILINEAR_INTERPOLATE_VECTORIZED(positions_x, positions_y, positions_z, $
    plasma_model.x_axis, plasma_model.y_axis, plasma_model.z_axis, plasma_model.nec)

  Te_los = IPT_TRILINEAR_INTERPOLATE_VECTORIZED(positions_x, positions_y, positions_z, $
    plasma_model.x_axis, plasma_model.y_axis, plasma_model.z_axis, plasma_model.Tec)

  ion_field = IPT_GET_ION_DENSITY_FIELD(plasma_model, species_key)
  n_ion_los = IPT_TRILINEAR_INTERPOLATE_VECTORIZED(positions_x, positions_y, positions_z, $
    plasma_model.x_axis, plasma_model.y_axis, plasma_model.z_axis, ion_field)

  ; Find valid points (positive densities and temperatures)
  valid = (ne_los GT 0) AND (Te_los GT 0) AND (n_ion_los GT 0)
  valid_idx = WHERE(valid, n_valid)

  IF n_valid EQ 0 THEN RETURN, DBLARR(n_lines)

  ; Get emission table for this species
  emiss_table_species = tables_single.emissivities.(species_idx)

  ; Initialize integrand array: (n_points, n_lines)
  integrand = DBLARR(n_points, n_lines)

  ; Calculate per-ion photon emission rates and brightness integrand at each valid point
  FOR ii = 0L, n_valid-1L DO BEGIN
    i = valid_idx[ii]

    ; Interpolate per ion photon emission rate for all lines at this (Te, ne)
    ; Units: photons s^-1 ion^-1
    per_ion_rate = IPT_INTERPOLATE_EMISSIVITY_2D_VECTORIZED(Te_los[i], ne_los[i], $
      tables_single.temp_arr, tables_single.dens_arr, emiss_table_species)

    ; Tables contain per-ion photon emission rates [photons s^-1 ion^-1]
    ; Brightness integrand = per_ion_rate * n_ion [photons s^-1 cm^-3]
    integrand[i, *] = per_ion_rate * n_ion_los[i]
  ENDFOR

  ; Integrate each line over LOS using Simpson's rule
  brightnesses = DBLARR(n_lines)
  FOR line = 0L, n_lines-1L DO BEGIN
    brightnesses[line] = RAYLEIGH_FACTOR * $
      IPT_SIMPSON_INTEGRATE(s_values, integrand[*, line])
  ENDFOR

  RETURN, brightnesses
END

;==============================================================================
FUNCTION IPT_INTEGRATE_SPECIES_EMISSION_DOUBLE, plasma_model, tables_single, $
  tables_double, positions_x, positions_y, positions_z, s_values, species_key
  ;+
  ; NAME:
  ;   IPT_INTEGRATE_SPECIES_EMISSION_DOUBLE
  ;
  ; PURPOSE:
  ;   Integrate all emission lines for a single species along line of sight
  ;   using double Maxwellian distribution with automatic fallback to single
  ;   Maxwellian for low densities or hot fractions. Returns brightness in
  ;   Rayleighs for each emission line of this species.
  ;
  ; CATEGORY:
  ;   Io Plasma Torus, Emission Modeling, Line-of-Sight Integration
  ;
  ; CALLING SEQUENCE:
  ;   brightnesses = IPT_INTEGRATE_SPECIES_EMISSION_DOUBLE(plasma_model, $
  ;                    tables_single, tables_double, positions_x, positions_y, $
  ;                    positions_z, s_values, species_key)
  ;
  ; INPUTS:
  ;   plasma_model  - Plasma model structure from IPT_LOAD_PLASMA_MODEL
  ;   tables_single - Single Maxwellian tables from IPT_LOAD_EMISSION_TABLES_SINGLE
  ;   tables_double - Double Maxwellian tables from IPT_LOAD_EMISSION_TABLES_DOUBLE
  ;   positions_x/y/z - Arrays of positions along LOS [R_J]
  ;   s_values      - Path length values [R_J]
  ;   species_key   - Species name ('SP', 'S2P', 'S3P', 'OP', 'O2P')
  ;
  ; OUTPUTS:
  ;   Returns 1D array of line brightnesses in Rayleighs for all emission
  ;   lines of this species.
  ;
  ; NOTES:
  ;   Falls back to single Maxwellian for points where feh < feh_min or
  ;   ne_total < ne_total_min or Teh< Teh_min from interp grids
  ;-
  COMPILE_OPT IDL2, HIDDEN

  ; Physical constants
  R_J_CM = 7.1492d9
  RAYLEIGH_FACTOR = 1.0d-6 * R_J_CM

  n_points = N_ELEMENTS(s_values)

  ; Get species tag index for both tables
  species_tags_double = TAG_NAMES(tables_double.wavelengths)
  species_idx_double = WHERE(species_tags_double EQ STRUPCASE(species_key), n_match_double)

  species_tags_single = TAG_NAMES(tables_single.wavelengths)
  species_idx_single = WHERE(species_tags_single EQ STRUPCASE(species_key), n_match_single)

  IF n_match_double EQ 0 OR n_match_single EQ 0 THEN BEGIN
    PRINT, '  Warning: Species ', species_key, ' not found in tables'
    RETURN, [0.0d]
  ENDIF

  n_lines = N_ELEMENTS(tables_double.wavelengths.(species_idx_double))

  ; Interpolate all plasma parameters along LOS
  ne_total_los = IPT_TRILINEAR_INTERPOLATE_VECTORIZED(positions_x, positions_y, positions_z, $
    plasma_model.x_axis, plasma_model.y_axis, plasma_model.z_axis, plasma_model.ne_total)

  nec_los = IPT_TRILINEAR_INTERPOLATE_VECTORIZED(positions_x, positions_y, positions_z, $
    plasma_model.x_axis, plasma_model.y_axis, plasma_model.z_axis, plasma_model.nec)

  feh_los = IPT_TRILINEAR_INTERPOLATE_VECTORIZED(positions_x, positions_y, positions_z, $
    plasma_model.x_axis, plasma_model.y_axis, plasma_model.z_axis, plasma_model.feh)

  Tec_los = IPT_TRILINEAR_INTERPOLATE_VECTORIZED(positions_x, positions_y, positions_z, $
    plasma_model.x_axis, plasma_model.y_axis, plasma_model.z_axis, plasma_model.Tec)

  Teh_los = IPT_TRILINEAR_INTERPOLATE_VECTORIZED(positions_x, positions_y, positions_z, $
    plasma_model.x_axis, plasma_model.y_axis, plasma_model.z_axis, plasma_model.Teh)

  ion_field = IPT_GET_ION_DENSITY_FIELD(plasma_model, species_key)
  n_ion_los = IPT_TRILINEAR_INTERPOLATE_VECTORIZED(positions_x, positions_y, positions_z, $
    plasma_model.x_axis, plasma_model.y_axis, plasma_model.z_axis, ion_field)

  ; Determine which points should use double vs single Maxwellian
  use_single = (feh_los LT tables_double.feh_arr[0]) OR (ne_total_los LT tables_double.ne_arr[0]) OR (Teh_los LT tables_double.Teh_arr[0])

  ; Valid points for each distribution type
  valid_double = (~use_single) AND (ne_total_los GT 0) AND (Tec_los GT 0) AND $
    (Teh_los GT 0) AND (feh_los GE 0) AND (feh_los LE 1) AND (n_ion_los GT 0)

  valid_single = use_single AND (nec_los GT 0) AND (Tec_los GT 0) AND (n_ion_los GT 0)

  ; Get emission tables
  emiss_table_4d = tables_double.emissivities.(species_idx_double)
  emiss_table_2d = tables_single.emissivities.(species_idx_single)

  ; Initialize integrand array
  integrand = DBLARR(n_points, n_lines)

  ; Process double Maxwellian points
  double_idx = WHERE(valid_double, n_double)
  IF n_double GT 0 THEN BEGIN
    FOR ii = 0L, n_double-1L DO BEGIN
      i = double_idx[ii]

      ; Interpolate Per-ion photon emission rate from 4D tables
      per_ion_rate = IPT_INTERPOLATE_EMISSIVITY_4D_VECTORIZED(Tec_los[i], Teh_los[i], $
        ne_total_los[i], feh_los[i], tables_double.tec_arr, tables_double.teh_arr, $
        tables_double.ne_arr, tables_double.feh_arr, emiss_table_4d)

      
      ; Tables contain per-ion emission rates [photons s^-1 ion^-1]
      integrand[i, *] = per_ion_rate  * n_ion_los[i]
    ENDFOR
  ENDIF

  ; Process single Maxwellian points
  single_idx = WHERE(valid_single, n_single)
  IF n_single GT 0 THEN BEGIN
    FOR ii = 0L, n_single-1L DO BEGIN
      i = single_idx[ii]

      ; Interpolate Per-ion photon emission rate from 2D tables
      per_ion_rate  = IPT_INTERPOLATE_EMISSIVITY_2D_VECTORIZED(Tec_los[i], nec_los[i], $
        tables_single.temp_arr, tables_single.dens_arr, emiss_table_2d)

 

      ; Brightness integrand
      ; Tables contain per-ion emission rates [photons s^-1 ion^-1]
      integrand[i, *] = per_ion_rate * n_ion_los[i]
    ENDFOR
  ENDIF

  ; Integrate each line over LOS
  brightnesses = DBLARR(n_lines)
  FOR line = 0L, n_lines-1L DO BEGIN
    brightnesses[line] = RAYLEIGH_FACTOR * $
      IPT_SIMPSON_INTEGRATE(s_values, integrand[*, line])
  ENDFOR

  RETURN, brightnesses
END

;==============================================================================
PRO IPT_CALCULATE_SPECTRUM_SINGLE, plasma_model, tables_single, $
  slit_pos_vec, norm_vec, wave_bins, spectrum, line_list, $
  WAVELENGTH_RANGE=wavelength_range, BIN_WIDTH=bin_width, FWHM=fwhm, DS=ds
  ;+
  ; NAME:
  ;   IPT_CALCULATE_SPECTRUM_SINGLE
  ;
  ; PURPOSE:
  ;   Calculate complete UV/optical spectrum for single Maxwellian distribution
  ;   with proper line-of-sight integration through 3D plasma model.
  ;
  ; CATEGORY:
  ;   Io Plasma Torus, Spectral Modeling, Ray Tracing
  ;
  ; CALLING SEQUENCE:
  ;   IPT_CALCULATE_SPECTRUM_SINGLE, plasma_model, tables_single, $
  ;     slit_pos_vec, norm_vec, wave_bins, spectrum, line_list, $
  ;     WAVELENGTH_RANGE=[min, max], BIN_WIDTH=bin_width, FWHM=fwhm, DS=ds
  ;
  ; INPUTS:
  ;   plasma_model  - Plasma model structure from IPT_LOAD_PLASMA_MODEL
  ;   tables_single - Single Maxwellian tables from IPT_LOAD_EMISSION_TABLES_SINGLE
  ;   slit_pos_vec  - Starting position [x, y, z] in R_J
  ;   norm_vec      - Direction vector [nx, ny, nz]
  ;
  ; KEYWORDS:
  ;   WAVELENGTH_RANGE - [min, max] wavelength in Angstroms (default: [550, 2100])
  ;   BIN_WIDTH        - Spectral bin width in Angstroms (default: 1.0)
  ;   FWHM             - Instrumental FWHM in Angstroms (default: 6.0)
  ;   DS               - Integration step size in R_J (default: 0.01)
  ;
  ; OUTPUTS:
  ;   wave_bins - Wavelength bin centers [Angstroms]
  ;   spectrum  - Spectrum in Rayleighs/Angstrom
  ;   line_list - List structure array with .wavelength, .brightness, .species
  ;-
  ; Set defaults
  IF N_ELEMENTS(wavelength_range) EQ 0 THEN wavelength_range = [550.0d, 2100.0d]
  IF N_ELEMENTS(bin_width) EQ 0 THEN bin_width = 1.0d
  IF N_ELEMENTS(fwhm) EQ 0 THEN fwhm = 6.0d
  IF N_ELEMENTS(ds) EQ 0 THEN ds = 0.01d

  ; Create wavelength bins
  n_bins = LONG((wavelength_range[1] - wavelength_range[0]) / bin_width)
  wave_bins = wavelength_range[0] + bin_width * DINDGEN(n_bins)
  spectrum = DBLARR(n_bins)

  ; Normalize direction vector
  norm = SQRT(TOTAL(norm_vec^2))
  direction = norm_vec / norm

  ; Create ray points along LOS
  max_distance = 40.0d
  n_points = LONG(max_distance / ds) + 1L
  s_values = ds * DINDGEN(n_points)

  ; Calculate positions along ray
  positions_x = slit_pos_vec[0] + s_values * direction[0]
  positions_y = slit_pos_vec[1] + s_values * direction[1]
  positions_z = slit_pos_vec[2] + s_values * direction[2]

  ; Collect all emission lines
  all_wavelengths = []
  all_brightnesses = []
  all_species = []

  ; Process each species
  species_list = ['SP', 'S2P', 'S3P', 'OP', 'O2P']

  FOR s = 0, N_ELEMENTS(species_list)-1 DO BEGIN
    species_key = species_list[s]

    ; Get wavelengths for this species
    species_tags = TAG_NAMES(tables_single.wavelengths)
    species_idx = WHERE(species_tags EQ STRUPCASE(species_key), n_match)
    IF n_match EQ 0 THEN CONTINUE

    wavelengths = tables_single.wavelengths.(species_idx)

    ; Filter to wavelength range
    in_range = WHERE(wavelengths GE wavelength_range[0] AND $
      wavelengths LE wavelength_range[1], n_in_range)

    IF n_in_range EQ 0 THEN CONTINUE

    PRINT, '  Processing ', species_key, ': ', STRTRIM(STRING(n_in_range),2), ' lines in range'

    ; Calculate brightness for all lines of this species
    brightnesses = IPT_INTEGRATE_SPECIES_EMISSION_SINGLE(plasma_model, tables_single, $
      positions_x, positions_y, positions_z, s_values, species_key)

    ; Add to collection
    FOR i = 0L, n_in_range-1L DO BEGIN
      line_idx = in_range[i]
      IF brightnesses[line_idx] GT 1d-10 THEN BEGIN
        all_wavelengths = [all_wavelengths, wavelengths[line_idx]]
        all_brightnesses = [all_brightnesses, brightnesses[line_idx]]
        all_species = [all_species, species_key]
      ENDIF
    ENDFOR
  ENDFOR

  n_lines_found = N_ELEMENTS(all_wavelengths)
  PRINT, '  Total lines with non-zero brightness: ', STRTRIM(STRING(n_lines_found),2)

  ; Build line list
  IF n_lines_found GT 0 THEN BEGIN
    line_list = REPLICATE({wavelength: 0.0d, brightness: 0.0d, species: ''}, n_lines_found)
    FOR i = 0L, n_lines_found-1L DO BEGIN
      line_list[i].wavelength = all_wavelengths[i]
      line_list[i].brightness = all_brightnesses[i]
      line_list[i].species = all_species[i]
    ENDFOR

    ; Convolve with instrumental response using ERF
    spectrum = IPT_SIMULATE_SPECTRUM_ERF(wave_bins, bin_width, $
      all_wavelengths, all_brightnesses, FWHM=fwhm)
  ENDIF ELSE BEGIN
    line_list = [{wavelength: 0.0d, brightness: 0.0d, species: ''}]
  ENDELSE

END

;==============================================================================
PRO IPT_CALCULATE_SPECTRUM_DOUBLE, plasma_model, tables_single, tables_double, $
  slit_pos_vec, norm_vec, wave_bins, spectrum, line_list, $
  WAVELENGTH_RANGE=wavelength_range, BIN_WIDTH=bin_width, FWHM=fwhm, DS=ds
  ;+
  ; NAME:
  ;   IPT_CALCULATE_SPECTRUM_DOUBLE
  ;
  ; PURPOSE:
  ;   Calculate complete UV/optical spectrum for double Maxwellian distribution
  ;   with proper line-of-sight integration through 3D plasma model. Automatically
  ;   falls back to single Maxwellian for low densities or hot fractions.
  ;
  ; CATEGORY:
  ;   Io Plasma Torus, Spectral Modeling, Ray Tracing
  ;
  ; CALLING SEQUENCE:
  ;   IPT_CALCULATE_SPECTRUM_DOUBLE, plasma_model, tables_single, tables_double, $
  ;     slit_pos_vec, norm_vec, wave_bins, spectrum, line_list, $
  ;     WAVELENGTH_RANGE=[min, max], BIN_WIDTH=bin_width, FWHM=fwhm, DS=ds
  ;
  ; INPUTS:
  ;   plasma_model  - Plasma model structure from IPT_LOAD_PLASMA_MODEL
  ;   tables_single - Single Maxwellian tables from IPT_LOAD_EMISSION_TABLES_SINGLE
  ;   tables_double - Double Maxwellian tables from IPT_LOAD_EMISSION_TABLES_DOUBLE
  ;   slit_pos_vec  - Starting position [x, y, z] in R_J
  ;   norm_vec      - Direction vector [nx, ny, nz]
  ;
  ; KEYWORDS:
  ;   WAVELENGTH_RANGE - [min, max] wavelength in Angstroms (default: [550, 2100])
  ;   BIN_WIDTH        - Spectral bin width in Angstroms (default: 1.0)
  ;   FWHM             - Instrumental FWHM in Angstroms (default: 6.0)
  ;   DS               - Integration step size in R_J (default: 0.01)
  ;
  ; OUTPUTS:
  ;   wave_bins - Wavelength bin centers [Angstroms]
  ;   spectrum  - Spectrum in Rayleighs/Angstrom
  ;   line_list - List structure array with .wavelength, .brightness, .species
  ;-
  ; Set defaults
  IF N_ELEMENTS(wavelength_range) EQ 0 THEN wavelength_range = [550.0d, 2100.0d]
  IF N_ELEMENTS(bin_width) EQ 0 THEN bin_width = 1.0d
  IF N_ELEMENTS(fwhm) EQ 0 THEN fwhm = 6.0d
  IF N_ELEMENTS(ds) EQ 0 THEN ds = 0.01d

  ; Create wavelength bins
  n_bins = LONG((wavelength_range[1] - wavelength_range[0]) / bin_width)
  wave_bins = wavelength_range[0] + bin_width * DINDGEN(n_bins)
  spectrum = DBLARR(n_bins)

  ; Normalize direction vector
  norm = SQRT(TOTAL(norm_vec^2))
  direction = norm_vec / norm

  ; Create ray points along LOS
  max_distance = 40.0d
  n_points = LONG(max_distance / ds) + 1L
  s_values = ds * DINDGEN(n_points)

  ; Calculate positions along ray
  positions_x = slit_pos_vec[0] + s_values * direction[0]
  positions_y = slit_pos_vec[1] + s_values * direction[1]
  positions_z = slit_pos_vec[2] + s_values * direction[2]

  ; Collect all emission lines
  all_wavelengths = []
  all_brightnesses = []
  all_species = []

  ; Process each species
  species_list = ['SP', 'S2P', 'S3P', 'OP', 'O2P']

  FOR s = 0, N_ELEMENTS(species_list)-1 DO BEGIN
    species_key = species_list[s]

    ; Get wavelengths for this species
    species_tags = TAG_NAMES(tables_double.wavelengths)
    species_idx = WHERE(species_tags EQ STRUPCASE(species_key), n_match)
    IF n_match EQ 0 THEN CONTINUE

    wavelengths = tables_double.wavelengths.(species_idx)

    ; Filter to wavelength range
    in_range = WHERE(wavelengths GE wavelength_range[0] AND $
      wavelengths LE wavelength_range[1], n_in_range)

    IF n_in_range EQ 0 THEN CONTINUE

    PRINT, '  Processing ', species_key, ': ', STRTRIM(STRING(n_in_range),2), ' lines in range'

    ; Calculate brightness for all lines of this species
    brightnesses = IPT_INTEGRATE_SPECIES_EMISSION_DOUBLE(plasma_model, tables_single, $
      tables_double, positions_x, positions_y, positions_z, s_values, species_key)

    ; Add to collection
    FOR i = 0L, n_in_range-1L DO BEGIN
      line_idx = in_range[i]
      IF brightnesses[line_idx] GT 1d-10 THEN BEGIN
        all_wavelengths = [all_wavelengths, wavelengths[line_idx]]
        all_brightnesses = [all_brightnesses, brightnesses[line_idx]]
        all_species = [all_species, species_key]
      ENDIF
    ENDFOR
  ENDFOR

  n_lines_found = N_ELEMENTS(all_wavelengths)
  PRINT, '  Total lines with non-zero brightness: ', STRTRIM(STRING(n_lines_found),2)

  ; Build line list
  IF n_lines_found GT 0 THEN BEGIN
    line_list = REPLICATE({wavelength: 0.0d, brightness: 0.0d, species: ''}, n_lines_found)
    FOR i = 0L, n_lines_found-1L DO BEGIN
      line_list[i].wavelength = all_wavelengths[i]
      line_list[i].brightness = all_brightnesses[i]
      line_list[i].species = all_species[i]
    ENDFOR

    ; Convolve with instrumental response using ERF
    spectrum = IPT_SIMULATE_SPECTRUM_ERF(wave_bins, bin_width, $
      all_wavelengths, all_brightnesses, FWHM=fwhm)
  ENDIF ELSE BEGIN
    line_list = [{wavelength: 0.0d, brightness: 0.0d, species: ''}]
  ENDELSE

END

;==============================================================================
PRO IPT_emiss_MOP_community_code
  ;+
  ; NAME:
  ;   IPT_emiss_MOP_community_code
  ;
  ; PURPOSE:
  ;   Dummy procedure to load the module. This allows the use of
  ;   @IPT_emiss_MOP_community_code syntax to compile all functions in this file.
  ;
  ; DESCRIPTION:
  ;   Prints a banner indicating that the IPT emission raytracer library
  ;   has been loaded and lists available functions.
  ;-

  PRINT, ''
  PRINT, '=========================================================================='
  PRINT, 'IPT UV/Optical Emission Raytracer Library Loaded'
  PRINT, 'MOP Community Code v1.0 - Proper Line-of-Sight Integration'
  PRINT, '=========================================================================='
  PRINT, ''
  PRINT, 'Available functions:'
  PRINT, '  IPT_LOAD_PLASMA_MODEL               - Load 3D plasma model from HDF5'
  PRINT, '  IPT_LOAD_EMISSION_TABLES_SINGLE     - Load single Maxwellian tables'
  PRINT, '  IPT_LOAD_EMISSION_TABLES_DOUBLE     - Load double Maxwellian tables'
  PRINT, '  IPT_CALCULATE_SPECTRUM_SINGLE       - Calculate spectrum (single Maxwellian)'
  PRINT, '  IPT_CALCULATE_SPECTRUM_DOUBLE       - Calculate spectrum (double Maxwellian)'
  PRINT, '  IPT_SIMULATE_SPECTRUM_ERF           - Convolve with instrument response'
  PRINT, '  IPT_SIMPSON_INTEGRATE               - Numerical integration'
  PRINT, ''
  PRINT, 'Line-of-Sight Integration Method:'
  PRINT, '  1. Ray tracing through 3D plasma model'
  PRINT, '  2. Trilinear interpolation of Te, ne, n_ion along LOS'
  PRINT, '  3. Per-ion photon emission rate from CHIANTI tables'
  PRINT, '  4. Simpson''s rule integration for brightness (Rayleighs)'
  PRINT, '  5. Gaussian instrumental response via analytic ERF'
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