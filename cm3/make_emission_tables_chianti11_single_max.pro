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

pro make_emission_tables_chianti11_single_max

  ; Density grid (cm^-3) - logarithmically spaced
  n_dens = 50
  dens_min = 0.1d
  dens_max = 10000.0d
  dens_arr = 10.0d^(dindgen(n_dens)/(n_dens-1) * alog10(dens_max/dens_min) + alog10(dens_min))

  ; Temperature grid (eV) - logarithmically spaced
  n_temp = 50
  temp_min = 0.1d
  temp_max = 500.0d
  temp_arr = 10.0d^(dindgen(n_temp)/(n_temp-1) * alog10(temp_max/temp_min) + alog10(temp_min))


nels = dens_arr

 yptsi = {s:dblarr(n_temp,n_dens,9),sp:dblarr(n_temp,n_dens,619),s2p:dblarr(n_temp,n_dens,355),s3p:dblarr(n_temp,n_dens,205),$
  s4p:dblarr(n_temp,n_dens,3157),o:dblarr(n_temp,n_dens,16),op:dblarr(n_temp,n_dens,4166),o2p:dblarr(n_temp,n_dens,3164),nap:dblarr(n_temp,n_dens,5054)}



for ii = 0, n_temp -1 do begin
  
  print," ii = ", ii + 1, " out of ", n_temp
   Tels = replicate( temp_arr[ii], n_dens)
  
  
  yptsi_temp = calculate_IPT_emiss_te_same_size_ne_all_species( Tels, nels, xwavi=xwavi)
  yptsi.s[ii,*,*] = yptsi_temp.s
  yptsi.sp[ii,*,*] = yptsi_temp.sp
  yptsi.s2p[ii,*,*] = yptsi_temp.s2p
  yptsi.s3p[ii,*,*] = yptsi_temp.s3p
  yptsi.s4p[ii,*,*] = yptsi_temp.s4p
  
  yptsi.o[ii,*,*] = yptsi_temp.o
  yptsi.op[ii,*,*] = yptsi_temp.op
  yptsi.o2p[ii,*,*] = yptsi_temp.o2p
  
  yptsi.nap[ii,*,*] = yptsi_temp.nap
endfor


SAVE, DENS_ARR,TEMP_ARR,XWAVI,YPTSI,/verbose, FILENAME = 'CHIANTI_11.0.2_emiss_arrays_all_species_all_wavelengths_50x50_logspaced.sav'

stop
end