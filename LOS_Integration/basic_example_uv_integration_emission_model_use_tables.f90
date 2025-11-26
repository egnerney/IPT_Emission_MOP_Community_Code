!==============================================================================
! basic_example_uv_integration_emission_model_use_tables.f90
!
! Example program demonstrating UV/optical emission calculations from the
! Io Plasma Torus using line-of-sight integration through a 3D plasma model.
!
! This example shows:
! - Loading plasma model and CHIANTI emission tables
! - Calculating UV spectra using single and double Maxwellian distributions
! - Comparing equatorial and off-equatorial lines of sight
! - Generating diagnostic plots with gnuplot
!
! AUTHOR: Edward (Eddie) G. Nerney
! INSTITUTION: Laboratory for Atmospheric and Space Physics, University of Colorado Boulder
! VERSION: 2.0
! DATE: November 2025
!==============================================================================

program uv_emission_example
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use ipt_emission_constants
    use ipt_emission_types
    use ipt_emission_utils
    use ipt_emission_raytracer
    use ipt_emission_plotting
    implicit none
    
    ! Raytracer object
    type(jovian_uv_raytracer) :: raytracer
    
    ! File paths
    character(len=256) :: plasma_file, emission_file_single, emission_file_double
    
    ! LOS configuration
    real(real64) :: slit_pos_equatorial(3), slit_pos_off_equator(3)
    real(real64) :: los_direction(3)
    real(real64) :: wavelength_range(2)
    real(real64) :: bin_width, fwhm, ds
    
    ! Output arrays
    real(real64), allocatable :: wave_bins(:), spectrum_single_eq(:), spectrum_single_off(:)
    real(real64), allocatable :: spectrum_double_eq(:), spectrum_double_off(:)
    type(emission_line), allocatable :: line_list(:)
    integer(int32) :: n_bins, n_lines
    
    ! Results
    real(real64) :: total_single_eq, total_single_off, total_double_eq, total_double_off
    real(real64) :: enhancement_eq, enhancement_off
    
    ! Timing
    real(real64) :: cpu_start, cpu_end
    
    ! Loop variables
    integer(int32) :: i, j, k
    type(emission_line), allocatable :: sorted_lines(:)
    real(real64) :: temp_br
    integer(int32) :: temp_sp, max_idx
    real(real64) :: temp_wav
    
    ! Plasma diagnostics
    real(real64), allocatable :: s_values(:), positions(:,:)
    real(real64), allocatable :: ne_los(:), Te_los(:)
    real(real64), allocatable :: nsp_los(:), ns2p_los(:), ns3p_los(:), nop_los(:), no2p_los(:)
    real(real64), allocatable :: feh_los(:), Teh_los(:)
    integer(int32) :: n_pos
    real(real64) :: rho
    
    call cpu_time(cpu_start)
    
    ! Print header
    print '(A)', '======================================================================'
    print '(A)', 'Jovian UV Emission Raytracer - Line-of-Sight Integration'
    print '(A)', '======================================================================'
    print '(A)', ''
    print '(A)', 'This example demonstrates UV emission calculations through'
    print '(A)', 'the Io Plasma Torus using proper line-of-sight integration'
    print '(A)', 'with species-by-species interpolation from CHIANTI tables.'
    print '(A)', ''
    
    ! Set file paths
    plasma_file = '../3D_Torus_Model/jovian_plasma_interpolated_381x381x231.h5'
    emission_file_single = '../Emiss_tables/CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5'
    emission_file_double = '../Emiss_tables/CHIANTI_11.0.2_emiss_tables_double_maxwellian_24x10x24x12.h5'
    
    ! Initialize raytracer
    print '(A)', 'Loading plasma model and emission tables...'
    call raytracer%initialize(plasma_file, emission_file_single, emission_file_double)
    
    if (.not. raytracer%plasma_loaded .or. .not. raytracer%single_maxwellian_loaded) then
        print '(A)', 'ERROR: Failed to initialize raytracer'
        stop 1
    end if
    
    ! Configure observation parameters
    slit_pos_equatorial = [6.0_real64, -20.0_real64, 0.0_real64]
    slit_pos_off_equator = [6.0_real64, -20.0_real64, 0.5_real64]
    los_direction = [0.0_real64, 1.0_real64, 0.0_real64]
    
    wavelength_range = [550.0_real64, 1450.0_real64]
    bin_width = 1.0_real64
    fwhm = 6.0_real64
    ds = 0.01_real64
    
    ! Allocate output arrays
    n_bins = int((wavelength_range(2) - wavelength_range(1)) / bin_width) + 1
    allocate(wave_bins(n_bins))
    allocate(spectrum_single_eq(n_bins), spectrum_single_off(n_bins))
    allocate(spectrum_double_eq(n_bins), spectrum_double_off(n_bins))
    allocate(line_list(MAX_LINES), sorted_lines(100))
    
    ! -------------------------------------------------------------------------
    ! Test Case 1: Equatorial LOS - Single Maxwellian
    ! -------------------------------------------------------------------------
    print '(A)', ''
    print '(A)', '======================================================================'
    print '(A)', 'Test Case 1: Equatorial Line of Sight (Single Maxwellian)'
    print '(A)', '----------------------------------------------------------------------'
    print '(A,F6.1,A,F6.1,A,F6.1,A)', 'Starting position: (', &
        slit_pos_equatorial(1), ', ', slit_pos_equatorial(2), ', ', slit_pos_equatorial(3), ') R_J'
    print '(A,F5.2,A,F5.2,A,F5.2,A)', 'Direction vector: (', &
        los_direction(1), ', ', los_direction(2), ', ', los_direction(3), ')'
    print '(A)', ''
    
    call raytracer%calculate_spectrum_single(slit_pos_equatorial, los_direction, &
                                              wavelength_range, bin_width, fwhm, ds, &
                                              wave_bins, spectrum_single_eq, line_list, n_bins, n_lines)
    
    total_single_eq = simpson_integrate(spectrum_single_eq(1:n_bins), wave_bins(1:n_bins))
    
    print '(A)', ''
    print '(A)', 'Single Maxwellian Results:'
    print '(A,F12.1,A)', '  Total brightness: ', total_single_eq, ' Rayleighs'
    print '(A,F10.2,A,F8.1,A)', '  Peak brightness: ', maxval(spectrum_single_eq), ' R/A at ', &
        wave_bins(maxloc(spectrum_single_eq, 1)), ' A'
    
    ! -------------------------------------------------------------------------
    ! Test Case 2: Off-Equator LOS - Single Maxwellian
    ! -------------------------------------------------------------------------
    print '(A)', ''
    print '(A)', '======================================================================'
    print '(A)', 'Test Case 2: Off-Equator Line of Sight (Single Maxwellian)'
    print '(A)', '----------------------------------------------------------------------'
    print '(A,F6.1,A,F6.1,A,F6.1,A)', 'Starting position: (', &
        slit_pos_off_equator(1), ', ', slit_pos_off_equator(2), ', ', slit_pos_off_equator(3), ') R_J'
    print '(A,F5.2,A,F5.2,A,F5.2,A)', 'Direction vector: (', &
        los_direction(1), ', ', los_direction(2), ', ', los_direction(3), ')'
    print '(A)', ''
    
    call raytracer%calculate_spectrum_single(slit_pos_off_equator, los_direction, &
                                              wavelength_range, bin_width, fwhm, ds, &
                                              wave_bins, spectrum_single_off, line_list, n_bins, n_lines)
    
    total_single_off = simpson_integrate(spectrum_single_off(1:n_bins), wave_bins(1:n_bins))
    
    print '(A)', ''
    print '(A)', 'Single Maxwellian Results:'
    print '(A,F12.1,A)', '  Total brightness: ', total_single_off, ' Rayleighs'
    print '(A,F10.2,A,F8.1,A)', '  Peak brightness: ', maxval(spectrum_single_off), ' R/A at ', &
        wave_bins(maxloc(spectrum_single_off, 1)), ' A'
    
    ! -------------------------------------------------------------------------
    ! Test Case 3: Equatorial LOS - Double Maxwellian
    ! -------------------------------------------------------------------------
    print '(A)', ''
    print '(A)', '======================================================================'
    print '(A)', 'Test Case 3: Equatorial LOS (Double Maxwellian)'
    print '(A)', '----------------------------------------------------------------------'
    
    if (raytracer%double_maxwellian_loaded) then
        call raytracer%calculate_spectrum_double(slit_pos_equatorial, los_direction, &
                                                  wavelength_range, bin_width, fwhm, ds, &
                                                  wave_bins, spectrum_double_eq, line_list, n_bins, n_lines)
        
        total_double_eq = simpson_integrate(spectrum_double_eq(1:n_bins), wave_bins(1:n_bins))
        if (total_single_eq > 0.0_real64) then
            enhancement_eq = total_double_eq / total_single_eq
        else
            enhancement_eq = 0.0_real64
        end if
        
        print '(A)', ''
        print '(A)', 'Double Maxwellian Results:'
        print '(A,F12.1,A)', '  Total brightness: ', total_double_eq, ' Rayleighs'
        print '(A,F10.2,A,F8.1,A)', '  Peak brightness: ', maxval(spectrum_double_eq), ' R/A at ', &
            wave_bins(maxloc(spectrum_double_eq, 1)), ' A'
        print '(A,F8.3)', '  Enhancement factor: ', enhancement_eq
    else
        spectrum_double_eq = 0.0_real64
        total_double_eq = 0.0_real64
        enhancement_eq = 0.0_real64
        print '(A)', 'Double Maxwellian tables not loaded - skipping'
    end if
    
    ! -------------------------------------------------------------------------
    ! Test Case 4: Off-Equator LOS - Double Maxwellian
    ! -------------------------------------------------------------------------
    print '(A)', ''
    print '(A)', '======================================================================'
    print '(A)', 'Test Case 4: Off-Equator LOS (Double Maxwellian)'
    print '(A)', '----------------------------------------------------------------------'
    
    if (raytracer%double_maxwellian_loaded) then
        call raytracer%calculate_spectrum_double(slit_pos_off_equator, los_direction, &
                                                  wavelength_range, bin_width, fwhm, ds, &
                                                  wave_bins, spectrum_double_off, line_list, n_bins, n_lines)
        
        total_double_off = simpson_integrate(spectrum_double_off(1:n_bins), wave_bins(1:n_bins))
        if (total_single_off > 0.0_real64) then
            enhancement_off = total_double_off / total_single_off
        else
            enhancement_off = 0.0_real64
        end if
        
        print '(A)', ''
        print '(A)', 'Double Maxwellian Results:'
        print '(A,F12.1,A)', '  Total brightness: ', total_double_off, ' Rayleighs'
        print '(A,F10.2,A,F8.1,A)', '  Peak brightness: ', maxval(spectrum_double_off), ' R/A at ', &
            wave_bins(maxloc(spectrum_double_off, 1)), ' A'
        print '(A,F8.3)', '  Enhancement factor: ', enhancement_off
    else
        spectrum_double_off = 0.0_real64
        total_double_off = 0.0_real64
        enhancement_off = 0.0_real64
        print '(A)', 'Double Maxwellian tables not loaded - skipping'
    end if
    
    ! -------------------------------------------------------------------------
    ! Generate Plots
    ! -------------------------------------------------------------------------
    print '(A)', ''
    print '(A)', '======================================================================'
    print '(A)', 'Generating Plots'
    print '(A)', '----------------------------------------------------------------------'
    
    call plot_spectrum_gnuplot(wave_bins, spectrum_single_eq, n_bins, line_list, n_lines, &
                               'Single Maxwellian - Equatorial (z=0)', 0, 'spectrum_single_equatorial')
    
    call plot_spectrum_gnuplot(wave_bins, spectrum_single_off, n_bins, line_list, n_lines, &
                               'Single Maxwellian - Off-Equator (z=0.5 R_J)', 1, 'spectrum_single_off_equator')
    
    call plot_spectrum_comparison_gnuplot(wave_bins, spectrum_single_eq, spectrum_single_off, n_bins, &
                                          'Equatorial (z=0)', 'Off-Equator (z=0.5)', &
                                          total_single_eq, total_single_off, &
                                          'Single Maxwellian: Equatorial vs Off-Equator', &
                                          'spectrum_single_comparison')
    
    if (raytracer%double_maxwellian_loaded) then
        call plot_spectrum_gnuplot(wave_bins, spectrum_double_eq, n_bins, line_list, n_lines, &
                                   'Double Maxwellian - Equatorial (z=0)', 0, 'spectrum_double_equatorial')
        
        call plot_spectrum_gnuplot(wave_bins, spectrum_double_off, n_bins, line_list, n_lines, &
                                   'Double Maxwellian - Off-Equator (z=0.5 R_J)', 1, 'spectrum_double_off_equator')
        
        call plot_single_vs_double_gnuplot(wave_bins, spectrum_single_eq, spectrum_double_eq, n_bins, &
                                           total_single_eq, total_double_eq, enhancement_eq, &
                                           'Single vs Double Maxwellian - Equatorial', &
                                           'spectrum_single_vs_double')
    end if
    
    ! -------------------------------------------------------------------------
    ! Sort and print strongest lines (from last calculation)
    ! -------------------------------------------------------------------------
    print '(A)', ''
    print '(A)', '======================================================================'
    print '(A)', 'Strongest UV Emission Lines (Equatorial, Single Maxwellian)'
    print '(A)', '----------------------------------------------------------------------'
    
    ! Recalculate to get line list
    call raytracer%calculate_spectrum_single(slit_pos_equatorial, los_direction, &
                                              wavelength_range, bin_width, fwhm, ds, &
                                              wave_bins, spectrum_single_eq, line_list, n_bins, n_lines)
    
    if (n_lines > 0) then
        ! Copy to sorted array
        do i = 1, min(n_lines, 100)
            sorted_lines(i) = line_list(i)
        end do
        
        ! Simple selection sort for top 20 by brightness
        do i = 1, min(20, min(n_lines, 100))
            max_idx = i
            do j = i+1, min(n_lines, 100)
                if (sorted_lines(j)%brightness > sorted_lines(max_idx)%brightness) max_idx = j
            end do
            if (max_idx /= i) then
                temp_wav = sorted_lines(i)%wavelength
                temp_br = sorted_lines(i)%brightness
                temp_sp = sorted_lines(i)%species_idx
                sorted_lines(i)%wavelength = sorted_lines(max_idx)%wavelength
                sorted_lines(i)%brightness = sorted_lines(max_idx)%brightness
                sorted_lines(i)%species_idx = sorted_lines(max_idx)%species_idx
                sorted_lines(max_idx)%wavelength = temp_wav
                sorted_lines(max_idx)%brightness = temp_br
                sorted_lines(max_idx)%species_idx = temp_sp
            end if
        end do
        
        print '(A)', ''
        print '(A8, A12, A12)', 'Species', 'Lambda [A]', 'Brightness [R]'
        print '(A)', '--------------------------------'
        do i = 1, min(20, n_lines)
            if (sorted_lines(i)%brightness > 0.01_real64) then
                print '(A8, F12.2, F12.2)', trim(get_display_name(sorted_lines(i)%species_idx)), &
                    sorted_lines(i)%wavelength, sorted_lines(i)%brightness
            end if
        end do
    end if
    
    ! -------------------------------------------------------------------------
    ! Diagnostic: Plasma parameters along ray
    ! -------------------------------------------------------------------------
    print '(A)', ''
    print '(A)', '======================================================================'
    print '(A)', 'Diagnostic Information: Plasma Parameters Along Ray (Equatorial)'
    print '(A)', '----------------------------------------------------------------------'
    
    allocate(s_values(MAX_LOS_POINTS), positions(MAX_LOS_POINTS, 3))
    allocate(ne_los(MAX_LOS_POINTS), Te_los(MAX_LOS_POINTS))
    allocate(nsp_los(MAX_LOS_POINTS), ns2p_los(MAX_LOS_POINTS), ns3p_los(MAX_LOS_POINTS))
    allocate(nop_los(MAX_LOS_POINTS), no2p_los(MAX_LOS_POINTS))
    allocate(feh_los(MAX_LOS_POINTS), Teh_los(MAX_LOS_POINTS))
    
    call raytracer%trace_ray(slit_pos_equatorial, los_direction, ds, 40.0_real64, s_values, positions, n_pos)
    call raytracer%interpolate_plasma_trilinear(positions, n_pos, ne_los, Te_los, &
                                                 nsp_los, ns2p_los, ns3p_los, nop_los, no2p_los, &
                                                 feh_los, Teh_los)
    
    print '(A)', ''
    print '(A)', 'Sampling plasma parameters at key positions:'
    print '(A)', '                Position (R_J) rho (R_J)  ne (cm^-3)   Te (eV)'
    print '(A)', '-----------------------------------------------------------------'
    
    do i = 1, n_pos, n_pos/5
        rho = sqrt(positions(i,1)**2 + positions(i,2)**2)
        print '(A,F6.1,A,F6.1,A,F6.1,A,F6.2,F12.2,F10.2)', '[', positions(i,1), ', ', &
            positions(i,2), ', ', positions(i,3), ']', rho, ne_los(i), Te_los(i)
    end do
    
    ! -------------------------------------------------------------------------
    ! Summary
    ! -------------------------------------------------------------------------
    print '(A)', ''
    print '(A)', '======================================================================'
    print '(A)', 'Physical Interpretation'
    print '(A)', '----------------------------------------------------------------------'
    print '(A)', 'The plasma torus shows:'
    print '(A)', '- Peak emission near 6 R_J (Io''s orbital radius)'
    print '(A)', '- Scale height of ~0.5-1 R_J'
    print '(A)', '- Dominant emission from S III and O II ions'
    print '(A)', '- Temperature ~5-10 eV in the cold torus'
    print '(A)', '- Electron density ~100-2000 cm^-3 at peak'
    print '(A)', '- Hot electrons enhance high-excitation lines'
    if (raytracer%double_maxwellian_loaded .and. enhancement_eq > 0.0_real64) then
        print '(A,F8.1,A)', '- Overall brightness enhancement: ', (enhancement_eq - 1.0_real64) * 100.0_real64, '%'
    end if
    
    ! Cleanup
    call raytracer%cleanup()
    deallocate(wave_bins, spectrum_single_eq, spectrum_single_off)
    deallocate(spectrum_double_eq, spectrum_double_off)
    deallocate(line_list, sorted_lines)
    deallocate(s_values, positions, ne_los, Te_los)
    deallocate(nsp_los, ns2p_los, ns3p_los, nop_los, no2p_los, feh_los, Teh_los)
    
    call cpu_time(cpu_end)
    
    print '(A)', ''
    print '(A)', '======================================================================'
    print '(A)', 'UV Emission Calculation Complete'
    print '(A)', '======================================================================'
    print '(A,F8.2,A)', 'Execution time: ', cpu_end - cpu_start, ' seconds'
    print '(A)', ''
    print '(A)', 'Output files generated:'
    print '(A)', '  - spectrum_single_equatorial.png'
    print '(A)', '  - spectrum_single_off_equator.png'
    print '(A)', '  - spectrum_single_comparison.png'
    print '(A)', '  - spectrum_double_equatorial.png'
    print '(A)', '  - spectrum_double_off_equator.png'
    print '(A)', '  - spectrum_single_vs_double.png'
    print '(A)', ''
    print '(A)', 'Calculation methodology:'
    print '(A)', '  - Ray tracing through 3D plasma model with trilinear interpolation'
    print '(A)', '  - Vectorized bilinear interpolation of single Maxwellian emission rates (2D)'
    print '(A)', '  - Vectorized quadrilinear interpolation of double Maxwellian emission rates (4D)'
    print '(A)', '  - Per-ion photon emission rates x ion density x 1e-6 -> Rayleighs'
    print '(A)', '  - Simpson''s rule for line-of-sight integration'
    print '(A)', '  - ERF-based Gaussian convolution for instrument response'
    print '(A,F6.2,A)', '  - Integration step size: ds = ', ds, ' R_J'
    print '(A,F6.1,A)', '  - Instrumental FWHM: ', fwhm, ' A'
    
end program uv_emission_example
