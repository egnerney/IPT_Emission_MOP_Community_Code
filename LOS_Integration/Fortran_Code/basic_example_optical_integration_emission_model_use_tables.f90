!==============================================================================
! basic_example_optical_integration_emission_model_use_tables.f90
!
! Example Optical Emission Calculation Using Line-of-Sight Integration
! =====================================================================
!
! This program demonstrates optical emission calculations through the Io Plasma
! Torus using ray tracing with species-by-species interpolation from
! CHIANTI emission tables.
!
! TEST CASES:
!   1. Single Maxwellian: equatorial line of sight (z=0)
!   2. Single Maxwellian: off-equator line of sight (z=0.5 R_J)
!   3. Double Maxwellian: equatorial line of sight
!   4. Double Maxwellian: off-equator line of sight
!
! WAVELENGTH RANGE:
!   3000-10000 Angstroms (ground-based optical telescope spectroscopy)
!
! KEY DIAGNOSTIC LINES:
!   - [S II] 6716/6731 A: Density diagnostic doublet
!   - [S III] 6312 A: Temperature diagnostic
!   - [O III] 4959/5007 A: Nebular doublet
!   - [O II] 3726/3729 A: Density diagnostic doublet
!
! DIRECTORY STRUCTURE EXPECTED:
!   IPT_Emission_MOP_Community_Code/
!     LOS_Integration/Fortran_Code/           (executable runs here)
!     Emiss_tables/                           (CHIANTI emission tables)
!     3D_Torus_Model/                         (3D plasma model)
!
! 
!
! OUTPUT FILES:
!   - optical_spectrum_single_equatorial.dat, .gp, .png
!   - optical_spectrum_single_off_equator.dat, .gp, .png
!   - optical_spectrum_single_comparison.dat, .gp, .png
!   - optical_spectrum_double_equatorial.dat, .gp, .png (if available)
!   - optical_spectrum_double_off_equator.dat, .gp, .png (if available)
!   - optical_spectrum_single_vs_double.dat, .gp, .png (if available)
!   - optical_sii_doublet_region.dat, .gp, .png
!   - optical_red_line_region.dat, .gp, .png
!   - optical_oiii_doublet_region.dat, .gp, .png
!
! AUTHOR: Edward (Eddie) G. Nerney
! INSTITUTION: Laboratory for Atmospheric and Space Physics, CU Boulder
! LICENSE: Open source for academic and research use (MOP Community Code)
! DATE: November 2025
!==============================================================================

program optical_emission_example
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use ipt_emission_constants
    use ipt_emission_utils
    use ipt_emission_types
    use ipt_emission_raytracer
    implicit none
    
    type(jovian_uv_raytracer) :: raytracer
    
    ! Observation geometry
    real(real64) :: slit_pos1(3), slit_pos2(3), norm_vec(3)
    real(real64) :: wave_range(2), bin_width, fwhm, ds
    
    ! Results - single Maxwellian
    real(real64), allocatable :: wave1_s(:), spectrum1_s(:)
    real(real64), allocatable :: wave2_s(:), spectrum2_s(:)
    type(emission_line_result), allocatable :: lines1_s(:), lines2_s(:)
    integer(int32) :: n_bins1_s, n_lines1_s, n_bins2_s, n_lines2_s
    
    ! Results - double Maxwellian
    real(real64), allocatable :: wave1_d(:), spectrum1_d(:)
    real(real64), allocatable :: wave2_d(:), spectrum2_d(:)
    type(emission_line_result), allocatable :: lines1_d(:), lines2_d(:)
    integer(int32) :: n_bins1_d, n_lines1_d, n_bins2_d, n_lines2_d
    
    ! Timing and statistics
    real(real64) :: total_bright1_s, total_bright2_s
    real(real64) :: total_bright1_d, total_bright2_d
    real(real64) :: enhancement1, enhancement2
    real(real64) :: start_time, end_time
    integer(int32) :: i, ierr
    
    ! Diagnostic line brightnesses
    real(real64) :: sii_6716, sii_6731, sii_ratio
    real(real64) :: oiii_5007, oiii_4959, oiii_ratio
    
    ! Get start time
    call cpu_time(start_time)
    
    ! Print header
    print '(A)', ''
    print '(A)', '======================================================================'
    print '(A)', 'OPTICAL EMISSION MODEL FOR IO PLASMA TORUS'
    print '(A)', 'Ground-Based Telescope Spectroscopy'
    print '(A)', 'Line-of-Sight Integration with CHIANTI 11.0.2'
    print '(A)', '======================================================================'
    print '(A)', ''
    print '(A)', '======================================================================'
    print '(A)', 'Jovian Optical Emission Raytracer - Line-of-Sight Integration'
    print '(A)', '======================================================================'
    print '(A)', ''
    print '(A)', 'This example demonstrates optical emission calculations through'
    print '(A)', 'the Io Plasma Torus using proper line-of-sight integration'
    print '(A)', 'with species-by-species interpolation from CHIANTI tables.'
    print '(A)', ''
    print '(A)', 'Wavelength range: 3000-10000 A (ground-based optical)'
    print '(A)', ''
    
    ! Initialize raytracer with default paths
    print '(A)', 'Loading plasma model and emission tables...'
    call raytracer%initialize(ierr=ierr)
    
    if (ierr /= 0) then
        print '(A)', ''
        print '(A)', 'ERROR: Failed to initialize raytracer.'
        print '(A)', 'Please ensure:'
        print '(A)', '  1. HDF5 files exist in expected directories'
        print '(A)', '  2. Run: python convert_hdf5_strings.py'
        stop 1
    end if
    
    ! Set observation parameters for optical wavelength range
    wave_range = [3000.0_real64, 10000.0_real64]  ! Optical range [Angstroms]
    bin_width = 0.6_real64                         ! Spectral bin width [Angstroms] - R~3000
    fwhm = 1.2_real64                              ! Instrumental FWHM [Angstroms] - ground-based
    ds = 0.1_real64                               ! Integration step size [R_J]
    
    ! Test case 1: Equatorial LOS (z=0)
    slit_pos1 = [6.0_real64, -20.0_real64, 0.0_real64]
    norm_vec = [0.0_real64, 1.0_real64, 0.0_real64]
    
    ! Test case 2: Off-equator LOS (z=0.5 R_J)
    slit_pos2 = [6.0_real64, -20.0_real64, 0.5_real64]
    
    !=========================================================================
    ! SINGLE MAXWELLIAN CALCULATIONS
    !=========================================================================
    
    print '(A)', ''
    print '(A)', '======================================================================'
    print '(A)', 'Test Case 1: Equatorial Line of Sight (Single Maxwellian)'
    print '(A)', '----------------------------------------------------------------------'
    write(*, '(A,3(F8.1,A))') 'Starting position: (', slit_pos1(1), ',', slit_pos1(2), ',', slit_pos1(3), ') R_J'
    write(*, '(A,3(F6.1,A))') 'Direction vector: (', norm_vec(1), ',', norm_vec(2), ',', norm_vec(3), ')'
    write(*, '(A,F7.0,A,F8.0,A)') 'Wavelength range: ', wave_range(1), ' - ', wave_range(2), ' A'
    print '(A)', ''
    
    call raytracer%calculate_spectrum_single(slit_pos1, norm_vec, wave_range, &
        bin_width, fwhm, ds, wave1_s, spectrum1_s, lines1_s, n_bins1_s, n_lines1_s)
    
    total_bright1_s = simpson_integrate_simple(spectrum1_s, wave1_s, n_bins1_s)
    
    print '(A)', ''
    print '(A)', 'Single Maxwellian Results:'
    write(*, '(A,F12.1,A)') '  Total brightness: ', total_bright1_s, ' Rayleighs'
    write(*, '(A,I0)') '  Number of emission lines: ', n_lines1_s
    
    ! Test case 2: Off-equator (single Maxwellian)
    print '(A)', ''
    print '(A)', '======================================================================'
    print '(A)', 'Test Case 2: Off-Equator Line of Sight (Single Maxwellian)'
    print '(A)', '----------------------------------------------------------------------'
    write(*, '(A,3(F8.1,A))') 'Starting position: (', slit_pos2(1), ',', slit_pos2(2), ',', slit_pos2(3), ') R_J'
    write(*, '(A,3(F6.1,A))') 'Direction vector: (', norm_vec(1), ',', norm_vec(2), ',', norm_vec(3), ')'
    print '(A)', ''
    
    call raytracer%calculate_spectrum_single(slit_pos2, norm_vec, wave_range, &
        bin_width, fwhm, ds, wave2_s, spectrum2_s, lines2_s, n_bins2_s, n_lines2_s)
    
    total_bright2_s = simpson_integrate_simple(spectrum2_s, wave2_s, n_bins2_s)
    
    print '(A)', ''
    print '(A)', 'Single Maxwellian Results:'
    write(*, '(A,F12.1,A)') '  Total brightness: ', total_bright2_s, ' Rayleighs'
    
    !=========================================================================
    ! DOUBLE MAXWELLIAN CALCULATIONS
    !=========================================================================
    
    enhancement1 = 1.0_real64
    enhancement2 = 1.0_real64
    total_bright1_d = 0.0_real64
    total_bright2_d = 0.0_real64
    
    if (raytracer%double_maxwellian_loaded) then
        ! Test case 3: Equatorial (double Maxwellian)
        print '(A)', ''
        print '(A)', '======================================================================'
        print '(A)', 'Test Case 3: Equatorial Line of Sight (Double Maxwellian)'
        print '(A)', '----------------------------------------------------------------------'
        
        call raytracer%calculate_spectrum_double(slit_pos1, norm_vec, wave_range, &
            bin_width, fwhm, ds, wave1_d, spectrum1_d, lines1_d, n_bins1_d, n_lines1_d)
        
        total_bright1_d = simpson_integrate_simple(spectrum1_d, wave1_d, n_bins1_d)
        
        if (total_bright1_s > 0.0_real64) then
            enhancement1 = total_bright1_d / total_bright1_s
        end if
        
        print '(A)', ''
        print '(A)', 'Double Maxwellian Results:'
        write(*, '(A,F12.1,A)') '  Total brightness: ', total_bright1_d, ' Rayleighs'
        write(*, '(A,F8.3)') '  Enhancement factor: ', enhancement1
        
        ! Test case 4: Off-equator (double Maxwellian)
        print '(A)', ''
        print '(A)', '======================================================================'
        print '(A)', 'Test Case 4: Off-Equator Line of Sight (Double Maxwellian)'
        print '(A)', '----------------------------------------------------------------------'
        
        call raytracer%calculate_spectrum_double(slit_pos2, norm_vec, wave_range, &
            bin_width, fwhm, ds, wave2_d, spectrum2_d, lines2_d, n_bins2_d, n_lines2_d)
        
        total_bright2_d = simpson_integrate_simple(spectrum2_d, wave2_d, n_bins2_d)
        
        if (total_bright2_s > 0.0_real64) then
            enhancement2 = total_bright2_d / total_bright2_s
        end if
        
        print '(A)', ''
        print '(A)', 'Double Maxwellian Results:'
        write(*, '(A,F12.1,A)') '  Total brightness: ', total_bright2_d, ' Rayleighs'
        write(*, '(A,F8.3)') '  Enhancement factor: ', enhancement2
    end if
    
    !=========================================================================
    ! PRINT STRONGEST EMISSION LINES
    !=========================================================================
    
    print '(A)', ''
    print '(A)', '======================================================================'
    print '(A)', 'Strongest Optical Emission Lines (Equatorial, Single Maxwellian)'
    print '(A)', '----------------------------------------------------------------------'
    
    if (n_lines1_s > 0) then
        call print_top_lines(lines1_s, n_lines1_s, 20)
    else
        print '(A)', '  No lines with non-zero brightness'
    end if
    
    !=========================================================================
    ! KEY DIAGNOSTIC LINE RATIOS
    !=========================================================================
    
    print '(A)', ''
    print '(A)', '======================================================================'
    print '(A)', 'Key Diagnostic Line Ratios'
    print '(A)', '----------------------------------------------------------------------'
    
    ! Find [S II] 6716/6731 ratio (density diagnostic)
    sii_6716 = find_line_brightness(lines1_s, n_lines1_s, 6716.4_real64, 3.0_real64)
    sii_6731 = find_line_brightness(lines1_s, n_lines1_s, 6730.8_real64, 3.0_real64)
    
    if (sii_6731 > 0.0_real64) then
        sii_ratio = sii_6716 / sii_6731
        print '(A)', '[S II] 6716/6731 ratio (density diagnostic):'
        write(*, '(A,F8.3)') '  Single Maxwellian: ', sii_ratio
        print '(A)', '  (Low density limit ~1.5, high density limit ~0.44)'
    end if
    
    ! Find [O III] 5007/4959 ratio (should be ~3)
    oiii_5007 = find_line_brightness(lines1_s, n_lines1_s, 5006.8_real64, 3.0_real64)
    oiii_4959 = find_line_brightness(lines1_s, n_lines1_s, 4958.9_real64, 3.0_real64)
    
    if (oiii_4959 > 0.0_real64) then
        oiii_ratio = oiii_5007 / oiii_4959
        print '(A)', ''
        print '(A)', '[O III] 5007/4959 ratio (should be ~3):'
        write(*, '(A,F8.3)') '  Single Maxwellian: ', oiii_ratio
    end if
    
    !=========================================================================
    ! GENERATE OUTPUT FILES
    !=========================================================================
    
    print '(A)', ''
    print '(A)', '======================================================================'
    print '(A)', 'Generating Optical Spectrum Plots'
    print '(A)', '----------------------------------------------------------------------'
    
    ! Write spectrum data files and gnuplot scripts
    call write_spectrum_file('optical_spectrum_single_equatorial', wave1_s, spectrum1_s, n_bins1_s, &
        'Single Maxwellian: Equatorial LOS - Optical (x=6 R_J, z=0)', '#1f77b4', total_bright1_s)
    
    call write_spectrum_file('optical_spectrum_single_off_equator', wave2_s, spectrum2_s, n_bins2_s, &
        'Single Maxwellian: Off-Equator LOS - Optical (x=6 R_J, z=0.5 R_J)', '#1f77b4', total_bright2_s)
    
    call write_comparison_file('optical_spectrum_single_comparison', &
        wave1_s, spectrum1_s, n_bins1_s, total_bright1_s, &
        wave2_s, spectrum2_s, n_bins2_s, total_bright2_s, &
        'Equatorial (z=0)', 'Off-equator (z=0.5 R_J)', '#1f77b4', '#ff7f0e', &
        'Single Maxwellian: Equatorial vs Off-Equator Optical Spectra')
    
    if (raytracer%double_maxwellian_loaded .and. allocated(spectrum1_d)) then
        call write_spectrum_file('optical_spectrum_double_equatorial', wave1_d, spectrum1_d, n_bins1_d, &
            'Double Maxwellian: Equatorial LOS - Optical', '#ff7f0e', total_bright1_d)
        
        call write_spectrum_file('optical_spectrum_double_off_equator', wave2_d, spectrum2_d, n_bins2_d, &
            'Double Maxwellian: Off-Equator LOS - Optical', '#ff7f0e', total_bright2_d)
        
        call write_single_vs_double_file('optical_spectrum_single_vs_double', &
            wave1_s, spectrum1_s, n_bins1_s, total_bright1_s, &
            wave1_d, spectrum1_d, n_bins1_d, total_bright1_d, &
            'Equatorial LOS: Single vs Double Maxwellian - Optical')
    end if
    
    !=========================================================================
    ! ZOOM PLOTS FOR KEY DIAGNOSTIC LINE REGIONS
    !=========================================================================
    
    print '(A)', ''
    print '(A)', 'Generating zoom plots for key diagnostic regions...'
    
    ! [S II] doublet region (6700-6750 A)
    call write_zoom_region_file('optical_sii_doublet_region', &
        wave1_s, spectrum1_s, n_bins1_s, &
        wave1_d, spectrum1_d, n_bins1_d, &
        6700.0_real64, 6750.0_real64, &
        '[S II] 6716/6731 A Doublet Region - Density Diagnostic', &
        raytracer%double_maxwellian_loaded .and. allocated(spectrum1_d))
    
    ! Red line region (6280-6330 A) - [S III] 6312
    call write_zoom_region_file('optical_red_line_region', &
        wave1_s, spectrum1_s, n_bins1_s, &
        wave1_d, spectrum1_d, n_bins1_d, &
        6280.0_real64, 6330.0_real64, &
        '[S III] 6312 A Region - Temperature Diagnostic', &
        raytracer%double_maxwellian_loaded .and. allocated(spectrum1_d))
    
    ! [O III] nebular doublet region (4940-5020 A)
    call write_zoom_region_file('optical_oiii_doublet_region', &
        wave1_s, spectrum1_s, n_bins1_s, &
        wave1_d, spectrum1_d, n_bins1_d, &
        4940.0_real64, 5020.0_real64, &
        '[O III] 4959/5007 A Nebular Doublet Region', &
        raytracer%double_maxwellian_loaded .and. allocated(spectrum1_d))
    
    !=========================================================================
    ! PHYSICAL INTERPRETATION
    !=========================================================================
    
    print '(A)', ''
    print '(A)', '======================================================================'
    print '(A)', 'Physical Interpretation - Optical Observations'
    print '(A)', '----------------------------------------------------------------------'
    print '(A)', 'The plasma torus optical emission shows:'
    print '(A)', '- Peak emission near 6 R_J (Ios orbital radius)'
    print '(A)', '- Scale height of ~0.5-1 R_J'
    print '(A)', '- Dominant optical emission from [S II] and [O II] forbidden lines'
    print '(A)', '- [S II] 6716/6731 doublet ratio sensitive to electron density'
    print '(A)', '- [O III] 5007/4959 ratio fixed by atomic physics (~3:1)'
    print '(A)', '- [S III] 6312 A line useful for temperature diagnostics'
    print '(A)', '- Temperature ~5-10 eV in the cold torus'
    print '(A)', '- Electron density ~100-2000 cm^-3 at peak'
    if (raytracer%double_maxwellian_loaded) then
        print '(A)', '- Hot electrons enhance high-excitation optical lines'
        write(*, '(A,F6.1,A)') '- Overall brightness enhancement: ', (enhancement1-1.0_real64)*100.0_real64, '%'
    end if
    
    !=========================================================================
    ! SUMMARY
    !=========================================================================
    
    call cpu_time(end_time)
    
    print '(A)', ''
    print '(A)', '======================================================================'
    print '(A)', 'Optical Emission Calculation Complete'
    print '(A)', '======================================================================'
    print '(A)', ''
    print '(A)', 'Total Integrated Brightness [Rayleighs]:'
    print '(A)', '                         Equatorial    Off-Equator'
    write(*, '(A,F14.1,F14.1)') '  Single Maxwellian: ', total_bright1_s, total_bright2_s
    
    if (raytracer%double_maxwellian_loaded .and. allocated(spectrum1_d)) then
        write(*, '(A,F14.1,F14.1)') '  Double Maxwellian: ', total_bright1_d, total_bright2_d
    end if
    
    print '(A)', ''
    write(*, '(A,F8.2,A)') 'Execution time: ', end_time - start_time, ' seconds'
    print '(A)', ''
    print '(A)', 'Output files generated:'
    print '(A)', '  - optical_spectrum_single_equatorial.png'
    print '(A)', '  - optical_spectrum_single_off_equator.png'
    print '(A)', '  - optical_spectrum_single_comparison.png'
    if (raytracer%double_maxwellian_loaded) then
        print '(A)', '  - optical_spectrum_double_equatorial.png'
        print '(A)', '  - optical_spectrum_double_off_equator.png'
        print '(A)', '  - optical_spectrum_single_vs_double.png'
    end if
    print '(A)', '  - optical_sii_doublet_region.png'
    print '(A)', '  - optical_red_line_region.png'
    print '(A)', '  - optical_oiii_doublet_region.png'
    print '(A)', ''
    print '(A)', 'Calculation methodology:'
    print '(A)', '  - Ray tracing through 3D plasma model with trilinear interpolation'
    print '(A)', '  - Vectorized bilinear interpolation of single Maxwellian emission rates (2D)'
    if (raytracer%double_maxwellian_loaded) then
        print '(A)', '  - Vectorized quadrilinear interpolation of double Maxwellian rates (4D)'
    end if
    print '(A)', '  - Per-ion photon emission rates x ion density x 1e-6 -> Rayleighs'
    print '(A)', '  - Simpsons rule for line-of-sight integration'
    print '(A)', '  - ERF-based Gaussian convolution for instrument response'
    write(*, '(A,F6.3,A)') '  - Integration step size: ds = ', ds, ' R_J'
    write(*, '(A,F5.1,A)') '  - Instrumental FWHM: ', fwhm, ' A (R~3000 ground-based)'
    write(*, '(A,F7.0,A,F8.0,A)') '  - Wavelength range: ', wave_range(1), '-', wave_range(2), ' A'
    print '(A)', ''
    
    ! Clean up
    call raytracer%cleanup()
    if (allocated(wave1_s)) deallocate(wave1_s)
    if (allocated(spectrum1_s)) deallocate(spectrum1_s)
    if (allocated(lines1_s)) deallocate(lines1_s)
    if (allocated(wave2_s)) deallocate(wave2_s)
    if (allocated(spectrum2_s)) deallocate(spectrum2_s)
    if (allocated(lines2_s)) deallocate(lines2_s)
    if (allocated(wave1_d)) deallocate(wave1_d)
    if (allocated(spectrum1_d)) deallocate(spectrum1_d)
    if (allocated(lines1_d)) deallocate(lines1_d)
    if (allocated(wave2_d)) deallocate(wave2_d)
    if (allocated(spectrum2_d)) deallocate(spectrum2_d)
    if (allocated(lines2_d)) deallocate(lines2_d)

contains

    function simpson_integrate_simple(y, x, n) result(integral)
        !----------------------------------------------------------------------
        ! Simple Simpson's rule integration for total brightness calculation.
        ! Matches scipy.integrate.simpson behavior with composite Simpson's rule
        ! and trapezoidal correction for even number of panels.
        !----------------------------------------------------------------------
        integer(int32), intent(in) :: n
        real(real64), intent(in) :: y(n), x(n)
        real(real64) :: integral
        integer(int32) :: i
        real(real64) :: h
        
        integral = 0.0_real64
        if (n < 2) return
        
        if (n == 2) then
            integral = 0.5_real64 * (y(1) + y(2)) * (x(2) - x(1))
            return
        end if
        
        do i = 1, n - 2, 2
            h = (x(i+2) - x(i)) / 2.0_real64
            integral = integral + (h / 3.0_real64) * (y(i) + 4.0_real64*y(i+1) + y(i+2))
        end do
        
        if (mod(n, 2) == 0) then
            h = x(n) - x(n-1)
            integral = integral + 0.5_real64 * (y(n-1) + y(n)) * h
        end if
        
    end function simpson_integrate_simple
    
    
    function find_line_brightness(lines, n_lines, target_wav, tolerance) result(brightness)
        !----------------------------------------------------------------------
        ! Find brightness of emission line closest to target wavelength.
        ! Returns 0 if no line found within tolerance.
        !----------------------------------------------------------------------
        type(emission_line_result), intent(in) :: lines(:)
        integer(int32), intent(in) :: n_lines
        real(real64), intent(in) :: target_wav, tolerance
        real(real64) :: brightness
        
        integer(int32) :: i
        real(real64) :: min_diff, diff
        
        brightness = 0.0_real64
        min_diff = tolerance + 1.0_real64
        
        do i = 1, n_lines
            diff = abs(lines(i)%wavelength - target_wav)
            if (diff < min_diff) then
                min_diff = diff
                if (diff < tolerance) then
                    brightness = lines(i)%brightness
                end if
            end if
        end do
        
    end function find_line_brightness
    
    
    subroutine print_top_lines(lines, n_lines, n_top)
        !----------------------------------------------------------------------
        ! Print the brightest emission lines sorted by brightness.
        !----------------------------------------------------------------------
        type(emission_line_result), intent(in) :: lines(:)
        integer(int32), intent(in) :: n_lines, n_top
        
        integer(int32) :: i, j, n_print
        integer(int32), allocatable :: sorted_idx(:)
        real(real64) :: max_bright
        integer(int32) :: max_idx
        logical, allocatable :: used(:)
        character(len=8) :: ion_name
        
        n_print = min(n_top, n_lines)
        allocate(sorted_idx(n_print), used(n_lines))
        used = .false.
        
        ! Simple selection sort for top n_print lines
        do i = 1, n_print
            max_bright = -1.0_real64
            max_idx = 0
            do j = 1, n_lines
                if (.not. used(j) .and. lines(j)%brightness > max_bright) then
                    max_bright = lines(j)%brightness
                    max_idx = j
                end if
            end do
            if (max_idx > 0) then
                sorted_idx(i) = max_idx
                used(max_idx) = .true.
            end if
        end do
        
        print '(A8,A12,A12)', 'Ion', 'Wavelength', 'Brightness'
        print '(A)', '-----------------------------------'
        
        do i = 1, n_print
            if (sorted_idx(i) > 0) then
                ion_name = get_display_name(lines(sorted_idx(i))%species_idx)
                write(*, '(A8,F10.2,A,F10.2,A)') trim(ion_name), &
                    lines(sorted_idx(i))%wavelength, ' A', &
                    lines(sorted_idx(i))%brightness, ' R'
            end if
        end do
        
        deallocate(sorted_idx, used)
        
    end subroutine print_top_lines
    
    
    subroutine write_spectrum_file(basename, wavelength, spectrum, n_bins, title, color, total)
        !----------------------------------------------------------------------
        ! Write spectrum data file and gnuplot script for PNG generation.
        !----------------------------------------------------------------------
        character(len=*), intent(in) :: basename, title, color
        real(real64), intent(in) :: wavelength(:), spectrum(:), total
        integer(int32), intent(in) :: n_bins
        
        integer :: unit_dat, unit_gp, i
        character(len=256) :: dat_file, gp_file, png_file
        
        dat_file = trim(basename) // '.dat'
        gp_file = trim(basename) // '.gp'
        png_file = trim(basename) // '.png'
        
        ! Write data file
        open(newunit=unit_dat, file=dat_file, status='replace')
        write(unit_dat, '(A)') '# Wavelength [A]    Brightness [R/A]'
        do i = 1, n_bins
            write(unit_dat, '(F12.3,ES16.6)') wavelength(i), spectrum(i)
        end do
        close(unit_dat)
        
        ! Write gnuplot script
        open(newunit=unit_gp, file=gp_file, status='replace')
        write(unit_gp, '(A)') 'set terminal pngcairo size 1200,600 enhanced font "Arial,12"'
        write(unit_gp, '(A,A,A)') 'set output "', trim(png_file), '"'
        write(unit_gp, '(A,A,A)') 'set title "', trim(title), '"'
        write(unit_gp, '(A)') 'set xlabel "Wavelength [Å]"'
        write(unit_gp, '(A)') 'set ylabel "Brightness [Rayleighs/Å]"'
        write(unit_gp, '(A)') 'set grid'
        write(unit_gp, '(A,F10.1,A,F10.1,A)') 'set xrange [', wavelength(1), ':', wavelength(n_bins), ']'
        write(unit_gp, '(A)') 'set yrange [0:*]'
        write(unit_gp, '(A,F10.1,A)') 'set label "Total: ', total, ' R" at graph 0.02, 0.95'
        write(unit_gp, '(A,A,A,A,A)') 'plot "', trim(dat_file), '" using 1:2 with lines lw 1 lc rgb "', &
            trim(color), '" notitle'
        close(unit_gp)
        
        ! Execute gnuplot
        call execute_command_line('gnuplot ' // trim(gp_file) // ' 2>/dev/null', wait=.true.)
        
        print '(A,A)', '  Generated: ', trim(png_file)
        
    end subroutine write_spectrum_file
    
    
    subroutine write_comparison_file(basename, wave1, spec1, n1, total1, &
            wave2, spec2, n2, total2, label1, label2, color1, color2, plot_title)
        !----------------------------------------------------------------------
        ! Write comparison plot with two spectra.
        !----------------------------------------------------------------------
        character(len=*), intent(in) :: basename, label1, label2, color1, color2, plot_title
        real(real64), intent(in) :: wave1(:), spec1(:), wave2(:), spec2(:)
        real(real64), intent(in) :: total1, total2
        integer(int32), intent(in) :: n1, n2
        
        integer :: unit_dat, unit_gp, i
        character(len=256) :: dat_file, gp_file, png_file
        
        dat_file = trim(basename) // '.dat'
        gp_file = trim(basename) // '.gp'
        png_file = trim(basename) // '.png'
        
        ! Write data file (two columns for each spectrum)
        open(newunit=unit_dat, file=dat_file, status='replace')
        write(unit_dat, '(A)') '# Wavelength [A]    Spec1 [R/A]    Spec2 [R/A]'
        do i = 1, min(n1, n2)
            write(unit_dat, '(F12.3,2ES16.6)') wave1(i), spec1(i), spec2(i)
        end do
        close(unit_dat)
        
        ! Write gnuplot script
        open(newunit=unit_gp, file=gp_file, status='replace')
        write(unit_gp, '(A)') 'set terminal pngcairo size 1200,600 enhanced font "Arial,12"'
        write(unit_gp, '(A,A,A)') 'set output "', trim(png_file), '"'
        write(unit_gp, '(A,A,A)') 'set title "', trim(plot_title), '"'
        write(unit_gp, '(A)') 'set xlabel "Wavelength [Å]"'
        write(unit_gp, '(A)') 'set ylabel "Brightness [Rayleighs/Å]"'
        write(unit_gp, '(A)') 'set grid'
        write(unit_gp, '(A,F10.1,A,F10.1,A)') 'set xrange [', wave1(1), ':', wave1(n1), ']'
        write(unit_gp, '(A)') 'set yrange [0:*]'
        write(unit_gp, '(A)') 'set key top right'
        write(unit_gp, '(A,A,A,F10.1,A)') 'set label "', trim(label1), ': ', total1, ' R" at graph 0.02, 0.95'
        write(unit_gp, '(A,A,A,F10.1,A)') 'set label "', trim(label2), ': ', total2, ' R" at graph 0.02, 0.90'
        write(unit_gp, '(A,A,A,A,A,A,A)') 'plot "', trim(dat_file), '" using 1:2 with lines lw 1 lc rgb "', &
            trim(color1), '" title "', trim(label1), '", \'
        write(unit_gp, '(A,A,A,A,A,A,A)') '     "', trim(dat_file), '" using 1:3 with lines lw 1 lc rgb "', &
            trim(color2), '" title "', trim(label2), '"'
        close(unit_gp)
        
        ! Execute gnuplot
        call execute_command_line('gnuplot ' // trim(gp_file) // ' 2>/dev/null', wait=.true.)
        
        print '(A,A)', '  Generated: ', trim(png_file)
        
    end subroutine write_comparison_file


    subroutine write_single_vs_double_file(basename, wave1, spec1, n1, total1, &
            wave2, spec2, n2, total2, plot_title)
        !----------------------------------------------------------------------
        ! Write single vs double Maxwellian comparison plot.
        !----------------------------------------------------------------------
        character(len=*), intent(in) :: basename, plot_title
        real(real64), intent(in) :: wave1(:), spec1(:), wave2(:), spec2(:)
        real(real64), intent(in) :: total1, total2
        integer(int32), intent(in) :: n1, n2
        
        integer :: unit_dat, unit_gp, i
        character(len=256) :: dat_file, gp_file, png_file
        real(real64) :: enhancement
        
        dat_file = trim(basename) // '.dat'
        gp_file = trim(basename) // '.gp'
        png_file = trim(basename) // '.png'
        
        ! Calculate enhancement factor
        if (total1 > 0.0_real64) then
            enhancement = total2 / total1
        else
            enhancement = 0.0_real64
        end if
        
        ! Write data file
        open(newunit=unit_dat, file=dat_file, status='replace')
        write(unit_dat, '(A)') '# Wavelength [A]    Single [R/A]    Double [R/A]'
        do i = 1, min(n1, n2)
            write(unit_dat, '(F12.3,2ES16.6)') wave1(i), spec1(i), spec2(i)
        end do
        close(unit_dat)
        
        ! Write gnuplot script
        open(newunit=unit_gp, file=gp_file, status='replace')
        write(unit_gp, '(A)') 'set terminal pngcairo size 1200,600 enhanced font "Arial,12"'
        write(unit_gp, '(A,A,A)') 'set output "', trim(png_file), '"'
        write(unit_gp, '(A,A,A)') 'set title "', trim(plot_title), '"'
        write(unit_gp, '(A)') 'set xlabel "Wavelength [Å]"'
        write(unit_gp, '(A)') 'set ylabel "Brightness [Rayleighs/Å]"'
        write(unit_gp, '(A)') 'set grid'
        write(unit_gp, '(A,F10.1,A,F10.1,A)') 'set xrange [', wave1(1), ':', wave1(n1), ']'
        write(unit_gp, '(A)') 'set yrange [0:*]'
        write(unit_gp, '(A)') 'set key top right'
        write(unit_gp, '(A,F10.1,A)') 'set label "Single: ', total1, ' R" at graph 0.02, 0.95'
        write(unit_gp, '(A,F10.1,A)') 'set label "Double: ', total2, ' R" at graph 0.02, 0.90'
        write(unit_gp, '(A,F6.3,A)') 'set label "Enhancement: ', enhancement, 'x" at graph 0.02, 0.85'
        write(unit_gp, '(A,A,A)') 'plot "', trim(dat_file), '" using 1:2 with lines lw 1 lc rgb "#1f77b4" title "Single Maxwellian", \'
        write(unit_gp, '(A,A,A)') '     "', trim(dat_file), '" using 1:3 with lines lw 1 lc rgb "#ff7f0e" title "Double Maxwellian"'
        close(unit_gp)
        
        ! Execute gnuplot
        call execute_command_line('gnuplot ' // trim(gp_file) // ' 2>/dev/null', wait=.true.)
        
        print '(A,A)', '  Generated: ', trim(png_file)
        
    end subroutine write_single_vs_double_file


    subroutine write_zoom_region_file(basename, wave_s, spec_s, n_s, &
            wave_d, spec_d, n_d, wav_min, wav_max, plot_title, have_double)
        !----------------------------------------------------------------------
        ! Write zoom plot for a specific wavelength region.
        ! Overlays single and double Maxwellian spectra if available.
        !----------------------------------------------------------------------
        character(len=*), intent(in) :: basename, plot_title
        real(real64), intent(in) :: wave_s(:), spec_s(:)
        real(real64), intent(in) :: wave_d(:), spec_d(:)
        real(real64), intent(in) :: wav_min, wav_max
        integer(int32), intent(in) :: n_s, n_d
        logical, intent(in) :: have_double
        
        integer :: unit_dat, unit_gp, i
        character(len=256) :: dat_file, gp_file, png_file
        integer(int32) :: i_start, i_end, n_zoom
        
        dat_file = trim(basename) // '.dat'
        gp_file = trim(basename) // '.gp'
        png_file = trim(basename) // '.png'
        
        ! Find indices for zoom region in single Maxwellian spectrum
        i_start = 1
        i_end = n_s
        do i = 1, n_s
            if (wave_s(i) >= wav_min .and. i_start == 1) i_start = i
            if (wave_s(i) <= wav_max) i_end = i
        end do
        n_zoom = i_end - i_start + 1
        
        if (n_zoom < 2) return
        
        ! Write data file
        open(newunit=unit_dat, file=dat_file, status='replace')
        if (have_double .and. n_d > 0) then
            write(unit_dat, '(A)') '# Wavelength [A]    Single [R/A]    Double [R/A]'
            do i = i_start, i_end
                if (i <= n_d) then
                    write(unit_dat, '(F12.3,2ES16.6)') wave_s(i), spec_s(i), spec_d(i)
                else
                    write(unit_dat, '(F12.3,2ES16.6)') wave_s(i), spec_s(i), 0.0_real64
                end if
            end do
        else
            write(unit_dat, '(A)') '# Wavelength [A]    Brightness [R/A]'
            do i = i_start, i_end
                write(unit_dat, '(F12.3,ES16.6)') wave_s(i), spec_s(i)
            end do
        end if
        close(unit_dat)
        
        ! Write gnuplot script
        open(newunit=unit_gp, file=gp_file, status='replace')
        write(unit_gp, '(A)') 'set terminal pngcairo size 1000,600 enhanced font "Arial,12"'
        write(unit_gp, '(A,A,A)') 'set output "', trim(png_file), '"'
        write(unit_gp, '(A,A,A)') 'set title "', trim(plot_title), '" font "Arial,14"'
        write(unit_gp, '(A)') 'set xlabel "Wavelength [Å]"'
        write(unit_gp, '(A)') 'set ylabel "Brightness [Rayleighs/Å]"'
        write(unit_gp, '(A)') 'set grid'
        write(unit_gp, '(A,F10.1,A,F10.1,A)') 'set xrange [', wav_min, ':', wav_max, ']'
        write(unit_gp, '(A)') 'set yrange [0:*]'
        write(unit_gp, '(A)') 'set key top right font "Arial,11"'
        
        ! Add vertical lines for key diagnostic wavelengths
        if (wav_min <= 6716.4_real64 .and. wav_max >= 6716.4_real64) then
            write(unit_gp, '(A)') 'set arrow from 6716.4, graph 0 to 6716.4, graph 1 nohead lt 0 lw 0.8 lc rgb "gray"'
            write(unit_gp, '(A)') 'set label "[S II]\n6716 Å" at 6716.4, graph 0.9 center font "Arial,9"'
        end if
        if (wav_min <= 6730.8_real64 .and. wav_max >= 6730.8_real64) then
            write(unit_gp, '(A)') 'set arrow from 6730.8, graph 0 to 6730.8, graph 1 nohead lt 0 lw 0.8 lc rgb "gray"'
            write(unit_gp, '(A)') 'set label "[S II]\n6731 Å" at 6730.8, graph 0.9 center font "Arial,9"'
        end if
        if (wav_min <= 6312.1_real64 .and. wav_max >= 6312.1_real64) then
            write(unit_gp, '(A)') 'set arrow from 6312.1, graph 0 to 6312.1, graph 1 nohead lt 0 lw 0.8 lc rgb "gray"'
            write(unit_gp, '(A)') 'set label "[S III]\n6312 Å" at 6312.1, graph 0.9 center font "Arial,9"'
        end if
        if (wav_min <= 4958.9_real64 .and. wav_max >= 4958.9_real64) then
            write(unit_gp, '(A)') 'set arrow from 4958.9, graph 0 to 4958.9, graph 1 nohead lt 0 lw 0.8 lc rgb "gray"'
            write(unit_gp, '(A)') 'set label "[O III]\n4959 Å" at 4958.9, graph 0.9 center font "Arial,9"'
        end if
        if (wav_min <= 5006.8_real64 .and. wav_max >= 5006.8_real64) then
            write(unit_gp, '(A)') 'set arrow from 5006.8, graph 0 to 5006.8, graph 1 nohead lt 0 lw 0.8 lc rgb "gray"'
            write(unit_gp, '(A)') 'set label "[O III]\n5007 Å" at 5006.8, graph 0.9 center font "Arial,9"'
        end if
        
        if (have_double .and. n_d > 0) then
            write(unit_gp, '(A,A,A)') 'plot "', trim(dat_file), &
                '" using 1:2 with lines lw 1.5 lc rgb "#1f77b4" title "Single Maxwellian", \'
            write(unit_gp, '(A,A,A)') '     "', trim(dat_file), &
                '" using 1:3 with lines lw 1.5 lc rgb "#ff7f0e" title "Double Maxwellian"'
        else
            write(unit_gp, '(A,A,A)') 'plot "', trim(dat_file), &
                '" using 1:2 with lines lw 1.5 lc rgb "#1f77b4" title "Single Maxwellian"'
        end if
        close(unit_gp)
        
        ! Execute gnuplot
        call execute_command_line('gnuplot ' // trim(gp_file) // ' 2>/dev/null', wait=.true.)
        
        print '(A,A)', '  Generated: ', trim(png_file)
        
    end subroutine write_zoom_region_file

end program optical_emission_example
