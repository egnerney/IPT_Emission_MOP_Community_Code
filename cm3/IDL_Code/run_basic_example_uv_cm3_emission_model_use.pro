pro run_basic_example_uv_cm3_emission_model_use
  ;+
  ; AUTHOR: Edward (Eddie) G. Nerney
  ; NAME:
  ;   run_basic_example_uv_cm3_emission_model_use
  ;
  ; PURPOSE:
  ;   Setup script to configure CHIANTI atomic database paths and run UV emission
  ;   simulations for the Io Plasma Torus (IPT) in Jupiter's magnetosphere.
  ;   This is the main driver routine that sets up the environment and executes
  ;   the spectral simulation.
  ;
  ; DESCRIPTION:
  ;   This procedure:
  ;   1. Configures CHIANTI database paths for your system (Windows/Mac/Linux)
  ;   2. Sets up IDL system variables needed by CHIANTI
  ;   3. Compiles and runs the main emission model simulation
  ;
  ; REQUIREMENTS:
  ;   - CHIANTI atomic database (version 11.0.2 or compatible)
  ;   - CHIANTI IDL procedures must be in the specified path
  ;   - IDL 8.0 or later recommended
  ;
  ; 

  ; ============================================================================
  ; STEP 1: DETERMINE CURRENT DIRECTORY AND OPERATING SYSTEM
  ; ============================================================================
  ; Get the current working directory where this script is located
  cd, current = current

  ; Determine the operating system to use appropriate path delimiters
  ; Windows uses backslash (\) and semicolon (;) for path separation
  ; Unix-like systems (Mac/Linux) use forward slash (/) and colon (:)
  if !version.os_family eq 'Windows' then begin
    ; Windows path construction
    chiantipath = current + '\chianti\'
    path_delimiter = ';'
  endif else begin
    ; macOS/Linux path construction 
    chiantipath = current + '/chianti/'
    path_delimiter = ':'
  endelse

  ; ============================================================================
  ; ALTERNATIVE CHIANTI LOCATION (if not in current directory)
  ; ============================================================================
  ; If your CHIANTI installation is located elsewhere, uncomment and modify:
  ;
  ; For Windows example:
  ; chiantipath = 'C:\Users\username\pro\chianti\'
  ;
  ; For macOS/Linux example:
  ; chiantipath = '/Users/username/pro/chianti/'
  ; chiantipath = '/home/username/chianti/'

  ; ============================================================================
  ; STEP 2: CONFIGURE CHIANTI SYSTEM VARIABLES
  ; ============================================================================
  ; Check if CHIANTI system variables have already been defined
  ; This prevents redundant initialization if running multiple times
  defsysv, '!chianti', exists = compile

  ; Define the !chianti system variable if not already compiled
  ; This variable indicates that CHIANTI is available in the IDL session
  if ~compile then defsysv, '!chianti', 1
  
  ; ============================================================================
  ; STEP 3: ADD CHIANTI PROCEDURES TO IDL PATH
  ; ============================================================================
  ; Create path string for CHIANTI IDL procedures
  ; The '+' prefix tells IDL to recursively search subdirectories
  if ~compile then chiantipros = path_delimiter + expand_path('+' + chiantipath + 'idl')

  ; Add CHIANTI procedures to the IDL search path
  ; This allows IDL to find and execute CHIANTI routines
  if ~compile then !path = !path + chiantipros

  ; ============================================================================
  ; STEP 4: SET UP CHIANTI DATABASE LOCATION
  ; ============================================================================
  ; Tell CHIANTI where to find the atomic database files
  ; The 'dbase' directory contains all the atomic data files needed for calculations
  if ~compile then use_chianti, chiantipath + 'dbase'

  ; ============================================================================
  ; STEP 5: COMPILE AND RUN THE MAIN SIMULATION
  ; ============================================================================
  ; Compile the main emission model procedure
  ; This ensures all dependencies are loaded before execution
  resolve_routine, 'basic_example_uv_cm3_emission_model_use'

  ; Execute the main UV emission simulation
  ; This will calculate and plot UV spectra for the Io Plasma Torus
  basic_example_uv_cm3_emission_model_use
  
end