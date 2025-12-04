function calculate_IPT_emiss_te_same_size_ne_all_species, Tel, nel, xwavi=xwavi
  ;+
  ; AUTHOR: Edward G. Nerney
  ; NAME:
  ;   calculate_IPT_emiss_te_same_size_ne_all_species
  ;
  ; PURPOSE:
  ;   Calculate UV emission line intensities for multiple temperature/density pairs
  ;   simultaneously. This is the vectorized version for efficiency when computing
  ;   multiple plasma conditions.
  ;
  ; INPUTS:
  ;   Tel    - Electron temperature array [eV]
  ;   nel    - Electron number density array [#/cm^3]
  ;   
  ;
  ; KEYWORDS:
  ;   xwavi  - OUTPUT: Array of emission line wavelengths [Angstroms]
  ;
 
  ; PHYSICS:
  ;   Uses CHIANTI emissivities (photons/s) 
  ;  
  ;-

  ; Get number of temperature/density conditions to calculate
  num_tel=n_elements(tel)

  ; UNIT CONVERSION: eV to Kelvin
  ; From NIST: 1 eV = 11604.51812 K (https://physics.nist.gov/cuu/Constants)
  conversion=11604.51812d

  ; Convert to log10 scale as required by CHIANTI routines
  log10tel=alog10(Tel*conversion)  ; Log10 of temperature in Kelvin
  log10nel=alog10(nel)              ; Log10 of electron density in cm^-3

  ; ============================================================================
  ; CALCULATE EMISSIVITIES USING CHIANTI
  ; ============================================================================
  ; emiss_calc returns emissivities in units of photons/s
  ; These are volume emissivities that must be integrated along line of sight
  ;
  ; Key parameters used:
  ;   /no_de     - Drops the hc/lambda factor in the computation of the
  ;   emissivities. Useful for emission measure analyses involving
  ;   photon fluxes
  ;   radt=1.d   - Specify background radiation temperature (default: 6000 K) set to 1 to neglect this
  ;   /quiet     - Suppress informational messages
  ;   /NOPROT    - Exclude proton collision rates
  ;   /NOIONREC  - Exclude ionization/recombination (use fixed ion fractions)
  ;   /NO_RREC   - Exclude radiative recombination


   num_species = 9 
  ; Calculate emissivities for Sulfur ions (S, S+, S++, S+++, S++++)
  ; Atomic number for Sulfur = 16

  s1em = emiss_calc(16, 1, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)
  s2em = emiss_calc(16, 2, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)
  s3em = emiss_calc(16, 3, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)
  s4em = emiss_calc(16, 4, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)
  s5em = emiss_calc(16, 5, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)

  ; Calculate emissivities for Oxygen ions (O, O+, O++)
  ; Atomic number for Oxygen = 8
  o1em = emiss_calc(8, 1, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)
  o2em = emiss_calc(8, 2, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)
  o3em = emiss_calc(8, 3, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)

  ; Na+
  na2em = emiss_calc(11, 2, temp = log10tel, dens = log10nel, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)

  ; ============================================================================
  ; ORGANIZE WAVELENGTHS AND ION NAMES
  ; ============================================================================
  ; Combine all wavelengths from different ions
  xwavi_S = s1em.lambda
  xwavi_Sp = s2em.lambda
  xwavi_S2p = s3em.lambda
  xwavi_S3p = s4em.lambda
  xwavi_S4p = s5em.lambda
  
  xwavi_O = o1em.lambda
  xwavi_Op = o2em.lambda
  xwavi_O2p = o3em.lambda
  
  xwavi_Nap = na2em.lambda


  
  xwavi = {s:dblarr(n_elements(xwavi_s)),sp:dblarr(n_elements(xwavi_sp)),s2p:dblarr(n_elements(xwavi_s2p)),s3p:dblarr(n_elements(xwavi_s3p)),$
  s4p:dblarr(n_elements(xwavi_s4p)),o:dblarr(n_elements(xwavi_o)),op:dblarr(n_elements(xwavi_op)),o2p:dblarr(n_elements(xwavi_o2p)),nap:dblarr(n_elements(xwavi_nap))}

  xwavi.s = xwavi_s
  xwavi.sp = xwavi_sp
  xwavi.s2p = xwavi_s2p
  xwavi.s3p = xwavi_s3p
  xwavi.s4p = xwavi_s4p

  xwavi.o = xwavi_o
  xwavi.op = xwavi_op
  xwavi.o2p = xwavi_o2p

  xwavi.nap = xwavi_nap


  ; ============================================================================

  ; ============================================================================
  ; CALCULATE LINE INTENSITIES
  ; ============================================================================
  ; Initialize output array [n_conditions, n_lines]
  yptsi_S=dblarr(num_tel,n_elements(xwavi_S))
  yptsi_Sp=dblarr(num_tel,n_elements(xwavi_Sp))
  yptsi_S2p=dblarr(num_tel,n_elements(xwavi_S2p))
  yptsi_S3p=dblarr(num_tel,n_elements(xwavi_S3p))
  yptsi_S4p=dblarr(num_tel,n_elements(xwavi_S4p))
  yptsi_O=dblarr(num_tel,n_elements(xwavi_O))
  yptsi_Op=dblarr(num_tel,n_elements(xwavi_Op))
  yptsi_O2p=dblarr(num_tel,n_elements(xwavi_O2p))
  yptsi_nap=dblarr(num_tel,n_elements(xwavi_Nap))
 
 
  yptsi = {s:dblarr(num_tel,n_elements(xwavi_S)),sp:dblarr(num_tel,n_elements(xwavi_Sp)),s2p:dblarr(num_tel,n_elements(xwavi_S2p)),s3p:dblarr(num_tel,n_elements(xwavi_S3p)),$
    s4p:dblarr(num_tel,n_elements(xwavi_S4p)),o:dblarr(num_tel,n_elements(xwavi_O)),op:dblarr(num_tel,n_elements(xwavi_op)),o2p:dblarr(num_tel,n_elements(xwavi_o2p)),nap:dblarr(num_tel,n_elements(xwavi_nap))}


  ; Loop through each temperature/density condition
  for i=0, num_tel - 1 do begin
    ; Extract emissivities for this condition and multiply 
    ; reform() extracts the [i,i] diagonal element (same Te and ne index)


    yptsi_s[i,*] =reform(s1em.em[i,i,*])
    yptsi_sp[i,*] =reform(s2em.em[i,i,*])
    yptsi_s2p[i,*] =reform(s3em.em[i,i,*])
    yptsi_s3p[i,*] =reform(s4em.em[i,i,*])
    yptsi_s4p[i,*] =reform(s5em.em[i,i,*])
    
    yptsi_o[i,*] =reform(o1em.em[i,i,*])
    yptsi_op[i,*] =reform(o2em.em[i,i,*])
    yptsi_o2p[i,*] =reform(o3em.em[i,i,*])
    
    yptsi_nap[i,*] =reform(na2em.em[i,i,*])
    
  endfor
  
  yptsi.s = yptsi_s
  yptsi.sp = yptsi_sp
  yptsi.s2p = yptsi_s2p
  yptsi.s3p = yptsi_s3p
  yptsi.s4p = yptsi_s4p
  
  yptsi.o = yptsi_o
  yptsi.op = yptsi_op
  yptsi.o2p = yptsi_o2p
  
  yptsi.nap = yptsi_nap
  
  
  

  return, yptsi
  ; Note: xwavi is returned as a keyword output
end

;==============================================================================
function calculate_IPT_emiss_citep2_double_max_for_tables, Tec,Teh, ne_total,feh, $
   xwavi=xwavi
  ;+
  ; NAME:
  ;   calculate_IPT_emiss_citep2_double_max_for_tables
  ;
  ; PURPOSE:
  ;   Calculate UV emission for a DOUBLE Maxwellian electron distribution.
  ;   This represents plasma with both cold/core and hot electron populations,
  ;   common in magnetospheric plasmas.
  ;
  ; INPUTS:
  ;   Tec      - Core (cold) electron temperature [eV]
  ;   Teh      - Hot electron temperature [eV]
  ;   ne_total - Total electron density [#/cm^3]
  ;   feh      - Fraction of hot electrons (0 to 1)
  ;
  ; KEYWORDS:
  ;   xwavi    - OUTPUT: Emission line wavelengths
  ;
  ; PHYSICS:
  ;   The double Maxwellian represents:
  ;   f(v) = fec * f_Maxwell(v,Tec) + feh * f_Maxwell(v,Teh)
  ;   where fec + feh = 1
  ;
  ;   This is important for:
  ;   - Magnetospheric plasmas with wave-particle heating
  ;   - Pickup ion populations
  ;   - Reconnection regions
  ;-


  ; Convert temperatures to Kelvin (log scale for CHIANTI)
  conversion=11604.51812d
  log10tec=alog10(Tec*conversion)
  log10teh=alog10(Teh*conversion)
  log10ne_total=alog10(ne_total)

  ; Calculate fraction of cold electrons
  fec = 1.d - feh

  ; ============================================================================
  ; SETUP FOR CHIANTI'S MULTI-TEMPERATURE CAPABILITY
  ; ============================================================================
  ; CHIANTI can handle multiple Maxwellian components using sum_mwl_coeff
  ; This properly weights the contribution from each temperature component
  log10_temps = [log10tec, log10teh]  ; Array of temperatures
  fracs = [fec, feh]                   ; Fractional weights (must sum to 1)

  ; ============================================================================
  ; CALCULATE WEIGHTED EMISSIVITIES
  ; ============================================================================
  ; The sum_mwl_coeff keyword tells CHIANTI to calculate:
  ; emiss_total = fec * emiss(Tec) + feh * emiss(Teh)
  ; This accounts for different excitation rates at different temperatures

  ; Sulfur ions with double Maxwellian
  s1em = emiss_calc(16, 1, temp = log10_temps, sum_mwl_coeff=fracs, dens = log10ne_total, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)
  s2em = emiss_calc(16, 2, temp = log10_temps, sum_mwl_coeff=fracs, dens = log10ne_total, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)
  s3em = emiss_calc(16, 3, temp = log10_temps, sum_mwl_coeff=fracs, dens = log10ne_total, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)
  s4em = emiss_calc(16, 4, temp = log10_temps, sum_mwl_coeff=fracs, dens = log10ne_total, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)
  s5em = emiss_calc(16, 5, temp = log10_temps, sum_mwl_coeff=fracs, dens = log10ne_total, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)

  ; Oxygen ions with double Maxwellian
  o1em = emiss_calc(8, 1, temp = log10_temps, sum_mwl_coeff=fracs, dens = log10ne_total, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)
  o2em = emiss_calc(8, 2, temp = log10_temps, sum_mwl_coeff=fracs, dens = log10ne_total, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)
  o3em = emiss_calc(8, 3, temp = log10_temps, sum_mwl_coeff=fracs, dens = log10ne_total, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)

  ; Na+
  na2em = emiss_calc(11, 2, temp = log10_temps, sum_mwl_coeff=fracs, dens = log10ne_total, /no_de, radt = 1.d, /quiet,/NOPROT,/NOIONREC,/NO_RREC)

  
  

  ; ============================================================================
  ; ORGANIZE WAVELENGTHS AND ION NAMES
  ; ============================================================================
  ; Combine all wavelengths from different ions
  xwavi_S = s1em.lambda
  xwavi_Sp = s2em.lambda
  xwavi_S2p = s3em.lambda
  xwavi_S3p = s4em.lambda
  xwavi_S4p = s5em.lambda

  xwavi_O = o1em.lambda
  xwavi_Op = o2em.lambda
  xwavi_O2p = o3em.lambda

  xwavi_Nap = na2em.lambda



  xwavi = {s:dblarr(n_elements(xwavi_s)),sp:dblarr(n_elements(xwavi_sp)),s2p:dblarr(n_elements(xwavi_s2p)),s3p:dblarr(n_elements(xwavi_s3p)),$
    s4p:dblarr(n_elements(xwavi_s4p)),o:dblarr(n_elements(xwavi_o)),op:dblarr(n_elements(xwavi_op)),o2p:dblarr(n_elements(xwavi_o2p)),nap:dblarr(n_elements(xwavi_nap))}

  xwavi.s = xwavi_s
  xwavi.sp = xwavi_sp
  xwavi.s2p = xwavi_s2p
  xwavi.s3p = xwavi_s3p
  xwavi.s4p = xwavi_s4p

  xwavi.o = xwavi_o
  xwavi.op = xwavi_op
  xwavi.o2p = xwavi_o2p

  xwavi.nap = xwavi_nap


  ; ============================================================================

  

  yptsi = {s:dblarr(n_elements(xwavi_S)),sp:dblarr(n_elements(xwavi_Sp)),s2p:dblarr(n_elements(xwavi_S2p)),s3p:dblarr(n_elements(xwavi_S3p)),$
    s4p:dblarr(n_elements(xwavi_S4p)),o:dblarr(n_elements(xwavi_O)),op:dblarr(n_elements(xwavi_op)),o2p:dblarr(n_elements(xwavi_o2p)),nap:dblarr(n_elements(xwavi_nap))}


  

  yptsi.s = reform(s1em.em)
  yptsi.sp = reform(s2em.em)
  yptsi.s2p = reform(s3em.em)
  yptsi.s3p = reform(s4em.em)
  yptsi.s4p = reform(s5em.em)

  yptsi.o = reform(o1em.em)
  yptsi.op = reform(o2em.em)
  yptsi.o2p = reform(o3em.em)

  yptsi.nap = reform(na2em.em)


  return, yptsi
end


pro make_emission_tables_chianti11_double_max

  
  ; Tec: Cold electron temperature [eV]
  ; 16 points covering 0.5-30 eV with denser sampling in 2-15 eV range
  n_tec = 16
  tec_arr = [0.5d, 1.0d, 1.5d, 2.0d, 2.8d, 3.8d, 5.0d, 6.5d, $
    8.5d, 11.0d, 14.0d, 18.0d, 22.0d, 26.0d, 30.0d, 35.0d]

  ; Teh: Hot electron temperature [eV]
  ; 8 points covering 25-500 eV
  n_teh = 8
  teh_arr = [25.0d, 45.0d, 80.0d, 140.0d, 220.0d, 320.0d, 420.0d, 500.0d]

  ; n_e: Total electron density [cm^-3]
  ; 12 points log-spaced from 1 to 6000 cm^-3
  ; Denser sampling in 100-3000 range where most torus emission occurs
  n_ne_total = 12
  ne_total_arr = [1.0d, 5.0d, 20.0d, 60.0d, 150.0d, 400.0d, $
    800.0d, 1400.0d, 2200.0d, 3200.0d, 4500.0d, 6000.0d]
  ne_total_log10_arr = alog10(ne_total_arr)

  ; feh: Hot electron fraction
  ; 10 points log-spaced from 0.0003 to 0.2 (extended for off-equator)
  n_feh = 10
  feh_min = 0.0003d
  feh_max = 0.2d
  feh_arr = 10.0d^(dindgen(n_feh)/(n_feh-1) * alog10(feh_max/feh_min) + alog10(feh_min))
  feh_log10_arr = alog10(feh_arr)





yptsi_temp = calculate_IPT_emiss_citep2_double_max_for_tables( 5d,270d, 2200d,0.0025, xwavi=xwavi)

n_xwavi_s = n_elements(xwavi.s)
n_xwavi_sp = n_elements(xwavi.sp)
n_xwavi_s2p = n_elements(xwavi.s2p)
n_xwavi_s3p = n_elements(xwavi.s3p)
n_xwavi_s4p = n_elements(xwavi.s4p)


n_xwavi_o = n_elements(xwavi.o)
n_xwavi_op = n_elements(xwavi.op)
n_xwavi_o2p = n_elements(xwavi.o2p)

n_xwavi_nap = n_elements(xwavi.nap)


 yptsi = {s:dblarr(n_tec,n_teh,n_ne_total,n_feh,n_xwavi_s),sp:dblarr(n_tec,n_teh,n_ne_total,n_feh,n_xwavi_sp),s2p:dblarr(n_tec,n_teh,n_ne_total,n_feh,n_xwavi_s2p),s3p:dblarr(n_tec,n_teh,n_ne_total,n_feh,n_xwavi_s3p),$
  s4p:dblarr(n_tec,n_teh,n_ne_total,n_feh,n_xwavi_s4p),o:dblarr(n_tec,n_teh,n_ne_total,n_feh,n_xwavi_o),op:dblarr(n_tec,n_teh,n_ne_total,n_feh,n_xwavi_op),o2p:dblarr(n_tec,n_teh,n_ne_total,n_feh,n_xwavi_o2p),nap:dblarr(n_tec,n_teh,n_ne_total,n_feh,n_xwavi_nap)}


for ii = 0, n_tec -1 do begin
  for jj = 0, n_teh -1 do begin
    print,ii+1,,' of ', n_tec,' and ',jj+1, ' of ', n_teh
    for kk = 0, n_ne_total -1 do begin
      for ll = 0, n_feh -1 do begin

  
  
   Tec =  tec_arr[ii]
   Teh =  teh_arr[jj]
   ne_total =  ne_total_arr[kk]
   feh =  feh_arr[ll]
  
  
  yptsi_temp = calculate_IPT_emiss_citep2_double_max_for_tables( Tec,Teh, ne_total,feh, xwavi=xwavi)
  yptsi.s[ii,jj,kk,ll,*] = yptsi_temp.s
  yptsi.sp[ii,jj,kk,ll,*] = yptsi_temp.sp
  yptsi.s2p[ii,jj,kk,ll,*] = yptsi_temp.s2p
  yptsi.s3p[ii,jj,kk,ll,*] = yptsi_temp.s3p
  yptsi.s4p[ii,jj,kk,ll,*] = yptsi_temp.s4p
  
  yptsi.o[ii,jj,kk,ll,*] = yptsi_temp.o
  yptsi.op[ii,jj,kk,ll,*] = yptsi_temp.op
  yptsi.o2p[ii,jj,kk,ll,*] = yptsi_temp.o2p
  
  yptsi.nap[ii,jj,kk,ll,*] = yptsi_temp.nap
      endfor
    endfor
  endfor
endfor

;emiss_table - 5D table [n_tec, n_teh, n_ne, n_feh, n_lines]
SAVE, ne_total_ARR,TEC_ARR,TEh_ARR,feh_ARR,XWAVI,YPTSI,/verbose, FILENAME = 'CHIANTI_11.0.2_emiss_arrays_all_species_all_wavelengths_' + strtrim(n_tec,1) + 'x' + strtrim(n_teh,1) + 'x' + strtrim(n_ne_total,1) + 'x' + strtrim(n_feh,1) + '_hote_logspaced.sav'

stop
end