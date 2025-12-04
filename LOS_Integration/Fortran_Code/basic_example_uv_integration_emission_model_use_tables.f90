!==============================================================================
! basic_example_uv_integration_emission_model_use_tables.f90
!
! Example UV Emission Calculation Using Line-of-Sight Integration
! ================================================================
!
! This program demonstrates UV emission calculations through the Io Plasma
! Torus using ray tracing with species-by-species interpolation from
! CHIANTI emission tables.
!
! TEST CASES:
!   1. Single Maxwellian: equatorial line of sight (z=0)
!   2. Single Maxwellian: off-equator line of sight (z=0.5 R_J)
!   3. Double Maxwellian: equatorial line of sight
!   4. Double Maxwellian: off-equator line of sight
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
!   - spectrum_single_equatorial.dat, .gp, .png
!   - spectrum_single_off_equator.dat, .gp, .png
!   - spectrum_single_comparison.dat, .gp, .png
!   - spectrum_double_equatorial.dat, .gp, .png (if available)
!   - spectrum_double_off_equator.dat, .gp, .png (if available)
!
! AUTHOR: Edward (Eddie) G. Nerney
! INSTITUTION: Laboratory for Atmospheric and Space Physics, CU Boulder
! LICENSE: Open source for academic and research use (MOP Community Code)
! DATE: November 2025
!==============================================================================

program uv_emission_example
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
    real(real64) :: start_time, end_time
    integer(int32) :: i, ierr
    
    ! Get start time
    call cpu_time(start_time)
    
    ! Print header
    print '(A)', ''
    print '(A)', '======================================================================'
    print '(A)', 'Jovian UV Emission Raytracer - Line-of-Sight Integration'
    print '(A)', '======================================================================'
    print '(A)', ''
    print '(A)', 'This example demonstrates UV emission calculations through'
    print '(A)', 'the Io Plasma Torus using proper line-of-sight integration'
    print '(A)', 'with species-by-species interpolation from CHIANTI tables.'
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
    
    ! Set observation parameters
    wave_range = [550.0_real64, 2100.0_real64]  ! UV range [Angstroms]
    bin_width = 1.0_real64                       ! Spectral bin width [Angstroms]
    fwhm = 6.0_real64                            ! Instrumental FWHM [Angstroms]
    ds = 0.1_real64                             ! Integration step size [R_J]
    
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
    print '(A)', ''
    
    call raytracer%calculate_spectrum_single(slit_pos1, norm_vec, wave_range, &
        bin_width, fwhm, ds, wave1_s, spectrum1_s, lines1_s, n_bins1_s, n_lines1_s)
    
    total_bright1_s = simpson_integrate_simple(spectrum1_s, wave1_s, n_bins1_s)
    
    print '(A)', ''
    print '(A)', 'Single Maxwellian Results:'
    write(*, '(A,F12.1,A)') '  Total brightness: ', total_bright1_s, ' Rayleighs'
    
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
    
    if (raytracer%double_maxwellian_loaded) then
        ! Test case 3: Equatorial (double Maxwellian)
        print '(A)', ''
        print '(A)', '======================================================================'
        print '(A)', 'Test Case 3: Equatorial Line of Sight (Double Maxwellian)'
        print '(A)', '----------------------------------------------------------------------'
        
        call raytracer%calculate_spectrum_double(slit_pos1, norm_vec, wave_range, &
            bin_width, fwhm, ds, wave1_d, spectrum1_d, lines1_d, n_bins1_d, n_lines1_d)
        
        total_bright1_d = simpson_integrate_simple(spectrum1_d, wave1_d, n_bins1_d)
        
        print '(A)', ''
        print '(A)', 'Double Maxwellian Results:'
        write(*, '(A,F12.1,A)') '  Total brightness: ', total_bright1_d, ' Rayleighs'
        
        ! Test case 4: Off-equator (double Maxwellian)
        print '(A)', ''
        print '(A)', '======================================================================'
        print '(A)', 'Test Case 4: Off-Equator Line of Sight (Double Maxwellian)'
        print '(A)', '----------------------------------------------------------------------'
        
        call raytracer%calculate_spectrum_double(slit_pos2, norm_vec, wave_range, &
            bin_width, fwhm, ds, wave2_d, spectrum2_d, lines2_d, n_bins2_d, n_lines2_d)
        
        total_bright2_d = simpson_integrate_simple(spectrum2_d, wave2_d, n_bins2_d)
        
        print '(A)', ''
        print '(A)', 'Double Maxwellian Results:'
        write(*, '(A,F12.1,A)') '  Total brightness: ', total_bright2_d, ' Rayleighs'
    end if
    
    !=========================================================================
    ! PRINT STRONGEST EMISSION LINES
    !=========================================================================
    
    print '(A)', ''
    print '(A)', '======================================================================'
    print '(A)', 'Strongest UV Emission Lines (Equatorial, Single Maxwellian)'
    print '(A)', '----------------------------------------------------------------------'
    
    if (n_lines1_s > 0) then
        call print_top_lines(lines1_s, n_lines1_s, 15)
    else
        print '(A)', '  No lines with non-zero brightness'
    end if
    
    !=========================================================================
    ! GENERATE OUTPUT FILES
    !=========================================================================
    
    print '(A)', ''
    print '(A)', '======================================================================'
    print '(A)', 'Generating Output Files'
    print '(A)', '----------------------------------------------------------------------'
    
    ! Write spectrum data files and gnuplot scripts
    call write_spectrum_file('spectrum_single_equatorial', wave1_s, spectrum1_s, n_bins1_s, &
        'Single Maxwellian: Equatorial LOS', '#1f77b4', total_bright1_s)
    
    call write_spectrum_file('spectrum_single_off_equator', wave2_s, spectrum2_s, n_bins2_s, &
        'Single Maxwellian: Off-Equator LOS', '#1f77b4', total_bright2_s)
    
    call write_comparison_file('spectrum_single_comparison', &
        wave1_s, spectrum1_s, n_bins1_s, total_bright1_s, &
        wave2_s, spectrum2_s, n_bins2_s, total_bright2_s, &
        'Equatorial', 'Off-Equator', '#1f77b4', '#ff7f0e')
    

    if (raytracer%double_maxwellian_loaded .and. allocated(spectrum1_d)) then
        call write_spectrum_file('spectrum_double_equatorial', wave1_d, spectrum1_d, n_bins1_d, &
            'Double Maxwellian: Equatorial LOS', '#ff7f0e', total_bright1_d)
        
        call write_spectrum_file('spectrum_double_off_equator', wave2_d, spectrum2_d, n_bins2_d, &
            'Double Maxwellian: Off-Equator LOS', '#ff7f0e', total_bright2_d)
        
        ! ADD THIS LINE:
        call write_single_vs_double_file('spectrum_single_vs_double', &
            wave1_s, spectrum1_s, n_bins1_s, total_bright1_s, &
            wave1_d, spectrum1_d, n_bins1_d, total_bright1_d)
    end if
    
    !=========================================================================
    ! SUMMARY
    !=========================================================================
    
    call cpu_time(end_time)
    
    print '(A)', ''
    print '(A)', '======================================================================'
    print '(A)', 'Summary'
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
    print '(A)', 'Calculation methodology:'
    print '(A)', '  - Ray tracing through 3D plasma model with trilinear interpolation'
    print '(A)', '  - Vectorized bilinear interpolation of single Maxwellian emission rates'
    print '(A)', '  - Vectorized quadrilinear interpolation of double Maxwellian rates'
    print '(A)', '  - Per-ion photon emission rates x ion density x 1e-6 -> Rayleighs'
    print '(A)', '  - Simpsons rule for line-of-sight integration'
    print '(A)', '  - ERF-based Gaussian convolution for instrument response'
    write(*, '(A,F6.3,A)') '  - Integration step size: ds = ', ds, ' R_J'
    write(*, '(A,F5.1,A)') '  - Instrumental FWHM: ', fwhm, ' A'
    print '(A)', ''
    print '(A)', 'UV Emission Calculation Complete'
    print '(A)', '======================================================================'
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
            wave2, spec2, n2, total2, label1, label2, color1, color2)
        !----------------------------------------------------------------------
        ! Write comparison plot with two spectra.
        !----------------------------------------------------------------------
        character(len=*), intent(in) :: basename, label1, label2, color1, color2
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
        write(unit_gp, '(A)') 'set title "Single Maxwellian: Equatorial vs Off-Equator"'
        write(unit_gp, '(A)') 'set xlabel "Wavelength [Å]"'
        write(unit_gp, '(A)') 'set ylabel "Brightness [Rayleighs/Å]"'
        write(unit_gp, '(A)') 'set grid'
        write(unit_gp, '(A,F10.1,A,F10.1,A)') 'set xrange [', wave1(1), ':', wave1(n1), ']'
        write(unit_gp, '(A)') 'set yrange [0:*]'
        write(unit_gp, '(A)') 'set key top right'
        write(unit_gp, '(A,F10.1,A)') 'set label "Equatorial: ', total1, ' R" at graph 0.02, 0.95'
        write(unit_gp, '(A,F10.1,A)') 'set label "Off-Equator: ', total2, ' R" at graph 0.02, 0.90'
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
            wave2, spec2, n2, total2)
        !----------------------------------------------------------------------
        ! Write single vs double Maxwellian comparison plot.
        !----------------------------------------------------------------------
        character(len=*), intent(in) :: basename
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
        write(unit_gp, '(A)') 'set title "Equatorial LOS: Single vs Double Maxwellian"'
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

end program uv_emission_example
