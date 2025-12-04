!==============================================================================
! IPT_emiss_MOP_community_code.f90
!
! Io Plasma Torus (IPT) UV/Optical Emission Line-of-Sight Integration Module
! ==========================================================================
!
! This module provides comprehensive functionality for calculating UV and 
! optical emission spectra from the Io Plasma Torus (IPT) using ray tracing 
! through a 3D plasma model with CHIANTI atomic emission tables.
!
! PHYSICAL MODEL:
!   - Plasma torus extends from ~5-10 R_J in cylindrical radius
!   - Peak emission near 6 R_J at centrifugal equator
!   - Per-ion photon emission rates from CHIANTI 11.0.2 atomic database
!   - Line-of-sight integration using Simpson's rule
!   - Optically thin plasma approximation (valid for IPT)
!
! COORDINATE SYSTEMS AND UNITS:
!   - Positions: Jupiter radii [R_J], System III right-handed Cartesian
!   - Temperatures: electron volts [eV]
!   - Densities: particles per cubic centimeter [cm^-3]
!   - Wavelengths: Angstroms [A]
!   - Brightnesses: Rayleighs [R] = 10^6 photons s^-1 cm^-2 (4pi sr)^-1
!
! SPECIES SUPPORTED:
!   - S+  (S II)  - Singly ionized sulfur       (index 1)
!   - S++ (S III) - Doubly ionized sulfur       (index 2)
!   - S3+ (S IV)  - Triply ionized sulfur       (index 3)
!   - O+  (O II)  - Singly ionized oxygen       (index 4)
!   - O++ (O III) - Doubly ionized oxygen       (index 5)
!
!
! AUTHOR: Edward (Eddie) G. Nerney
! INSTITUTION: Laboratory for Atmospheric and Space Physics, CU Boulder
! LICENSE: Open source for academic and research use (MOP Community Code)
! VERSION: 1.0
! DATE: November 2025
!==============================================================================

module ipt_emission_constants
    !--------------------------------------------------------------------------
    ! Physical constants and array size limits for IPT emission calculations.
    !--------------------------------------------------------------------------
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none
    
    ! Physical constants
    real(real64), parameter :: R_J_KM = 71492.0_real64          ! Jupiter radius [km]
    real(real64), parameter :: R_J_CM = 7.1492e9_real64         ! Jupiter radius [cm]
    real(real64), parameter :: RAYLEIGH_FACTOR = 1.0e-6_real64 * R_J_CM  ! Conversion to Rayleighs
    real(real64), parameter :: PI = 3.14159265358979323846_real64
    real(real64), parameter :: SQRT2 = 1.41421356237309504880_real64
    real(real64), parameter :: LN2 = 0.69314718055994530942_real64
    real(real64), parameter :: FWHM_TO_SIGMA = 1.0_real64 / (2.0_real64 * sqrt(2.0_real64 * LN2))
    
    ! Array size limits
    integer(int32), parameter :: MAX_LINES = 100000
    integer(int32), parameter :: MAX_LOS_POINTS = 10000
    integer(int32), parameter :: MAX_WAVE_BINS = 10000
    integer(int32), parameter :: N_SPECIES = 5
    
    ! Species indices
    integer(int32), parameter :: IDX_SP = 1    ! S+
    integer(int32), parameter :: IDX_S2P = 2   ! S++
    integer(int32), parameter :: IDX_S3P = 3   ! S+++
    integer(int32), parameter :: IDX_OP = 4    ! O+
    integer(int32), parameter :: IDX_O2P = 5   ! O++
    
    ! Species key strings
    character(len=4), parameter :: SPECIES_KEYS(N_SPECIES) = ['SP  ', 'S2P ', 'S3P ', 'OP  ', 'O2P ']
    
end module ipt_emission_constants


module ipt_emission_utils
    !--------------------------------------------------------------------------
    ! Utility functions for numerical integration and species handling.
    !--------------------------------------------------------------------------
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use ipt_emission_constants
    implicit none
    
contains

    pure function simpson_integrate(y, x, n) result(integral)
        !----------------------------------------------------------------------
        ! Simpson's rule integration matching scipy.integrate.simpson behavior.
        ! Uses composite Simpson's rule with trapezoidal correction for even n.
        !----------------------------------------------------------------------
        integer(int32), intent(in) :: n
        real(real64), intent(in) :: y(n), x(n)
        real(real64) :: integral
        integer(int32) :: i
        real(real64) :: h
        
        integral = 0.0_real64
        if (n < 2) return
        
        ! For n=2, use trapezoidal rule
        if (n == 2) then
            integral = 0.5_real64 * (y(1) + y(2)) * (x(2) - x(1))
            return
        end if
        
        ! Composite Simpson's rule for odd number of points
        do i = 1, n - 2, 2
            h = (x(i+2) - x(i)) / 2.0_real64
            integral = integral + (h / 3.0_real64) * (y(i) + 4.0_real64*y(i+1) + y(i+2))
        end do
        
        ! Add trapezoidal correction for even n (last panel)
        if (mod(n, 2) == 0) then
            h = x(n) - x(n-1)
            integral = integral + 0.5_real64 * (y(n-1) + y(n)) * h
        end if
        
    end function simpson_integrate


    pure function searchsorted(arr, val, n) result(idx)
        !----------------------------------------------------------------------
        ! Binary search matching numpy.searchsorted behavior.
        ! Returns insertion point index (1-based) for val in sorted array arr.
        !----------------------------------------------------------------------
        integer(int32), intent(in) :: n
        real(real64), intent(in) :: arr(n)
        real(real64), intent(in) :: val
        integer(int32) :: idx
        integer(int32) :: lo, hi, mid
        
        if (val <= arr(1)) then
            idx = 1
            return
        end if
        
        if (val >= arr(n)) then
            idx = n + 1
            return
        end if
        
        lo = 1
        hi = n + 1
        
        do while (hi - lo > 1)
            mid = (lo + hi) / 2
            if (arr(mid) <= val) then
                lo = mid
            else
                hi = mid
            end if
        end do
        
        idx = hi
        
    end function searchsorted


    pure function species_string_to_index(sp_str) result(idx)
        !----------------------------------------------------------------------
        ! Convert species string to index. Handles various formats:
        ! "SP", "S+", "S II" -> 1
        ! "S2P", "S++", "S III" -> 2
        ! etc.
        !----------------------------------------------------------------------
        character(len=*), intent(in) :: sp_str
        integer(int32) :: idx
        character(len=8) :: sp_clean
        integer(int32) :: i, j, c
        
        ! Clean string: uppercase, remove spaces
        sp_clean = ''
        j = 0
        do i = 1, min(len_trim(sp_str), 8)
            c = ichar(sp_str(i:i))
            if (c >= ichar('a') .and. c <= ichar('z')) then
                j = j + 1
                sp_clean(j:j) = char(c - 32)  ! Convert to uppercase
            else if ((c >= ichar('A') .and. c <= ichar('Z')) .or. &
                     (c >= ichar('0') .and. c <= ichar('9')) .or. &
                     c == ichar('+')) then
                j = j + 1
                sp_clean(j:j) = sp_str(i:i)
            end if
        end do
        
        select case(trim(sp_clean))
        case('SP', 'S+', 'SII', 'S2')
            idx = IDX_SP
        case('S2P', 'S++', 'SIII', 'S3')
            idx = IDX_S2P
        case('S3P', 'S+++', 'SIV', 'S4')
            idx = IDX_S3P
        case('OP', 'O+', 'OII', 'O2')
            idx = IDX_OP
        case('O2P', 'O++', 'OIII', 'O3')
            idx = IDX_O2P
        case default
            idx = 0
        end select
        
    end function species_string_to_index


    pure function get_species_key(idx) result(key)
        !----------------------------------------------------------------------
        ! Get species key string from index.
        !----------------------------------------------------------------------
        integer(int32), intent(in) :: idx
        character(len=4) :: key
        
        if (idx >= 1 .and. idx <= N_SPECIES) then
            key = SPECIES_KEYS(idx)
        else
            key = '    '
        end if
    end function get_species_key


    pure function get_display_name(idx) result(name)
        !----------------------------------------------------------------------
        ! Get display name for species (e.g., "S II", "O III").
        !----------------------------------------------------------------------
        integer(int32), intent(in) :: idx
        character(len=8) :: name
        
        select case(idx)
        case(IDX_SP);  name = 'S II'
        case(IDX_S2P); name = 'S III'
        case(IDX_S3P); name = 'S IV'
        case(IDX_OP);  name = 'O II'
        case(IDX_O2P); name = 'O III'
        case default;  name = 'Unknown'
        end select
        
    end function get_display_name

end module ipt_emission_utils


module ipt_emission_types
    !--------------------------------------------------------------------------
    ! Derived type definitions for emission data structures.
    !--------------------------------------------------------------------------
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use ipt_emission_constants
    implicit none
    
    type :: species_emission_data
        integer(int32) :: n_lines = 0
        real(real64), allocatable :: wavelengths(:)
        real(real64), allocatable :: emiss_single(:,:,:)      ! (n_lines, n_T, n_n)
        real(real64), allocatable :: emiss_double(:,:,:,:,:)  ! (n_lines, n_Tc, n_Th, n_ne, n_feh)
    end type species_emission_data
    
    type :: emission_line_result
        real(real64) :: wavelength = 0.0_real64
        real(real64) :: brightness = 0.0_real64
        integer(int32) :: species_idx = 0
    end type emission_line_result
    
end module ipt_emission_types


module ipt_emission_io
    !--------------------------------------------------------------------------
    ! HDF5 I/O routines for reading plasma models and emission tables.
    ! Requires emission table files to be preprocessed with convert_hdf5_strings.py
    ! to convert variable-length strings to fixed-length format.
    !--------------------------------------------------------------------------
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding
    use hdf5
    use ipt_emission_constants
    use ipt_emission_utils
    implicit none
    
contains

    subroutine read_1d_real_dataset(file_id, dset_name, data_out, n_out, ierr)
        !----------------------------------------------------------------------
        ! Read 1D real64 dataset from HDF5 file.
        !----------------------------------------------------------------------
        integer(hid_t), intent(in) :: file_id
        character(len=*), intent(in) :: dset_name
        real(real64), allocatable, intent(out) :: data_out(:)
        integer(int32), intent(out) :: n_out
        integer(int32), intent(out) :: ierr
        
        integer(hid_t) :: dset_id, dspace_id
        integer(hsize_t) :: dims(1), maxdims(1)
        integer :: hdferr, rank
        
        ierr = 0
        n_out = 0
        
        call h5dopen_f(file_id, trim(dset_name), dset_id, hdferr)
        if (hdferr /= 0) then
            ierr = -1
            return
        end if
        
        call h5dget_space_f(dset_id, dspace_id, hdferr)
        call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, rank)
        n_out = int(dims(1), int32)
        
        allocate(data_out(n_out))
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, dims, hdferr)
        
        call h5sclose_f(dspace_id, hdferr)
        call h5dclose_f(dset_id, hdferr)
        
    end subroutine read_1d_real_dataset


    subroutine read_species_strings_fixed(file_id, dset_name, n_strings, species_idx, ierr)
        !----------------------------------------------------------------------
        ! Read fixed-length species strings from HDF5 file.
        ! Requires preprocessing with convert_hdf5_strings.py.
        !----------------------------------------------------------------------
        integer(hid_t), intent(in) :: file_id
        character(len=*), intent(in) :: dset_name
        integer(int32), intent(in) :: n_strings
        integer(int32), intent(out) :: species_idx(:)
        integer(int32), intent(out) :: ierr
        
        integer(hid_t) :: dset_id, dspace_id, dtype_id, memtype_id
        integer :: hdferr, i
        integer(hsize_t) :: dims(1)
        integer(size_t) :: str_size
        character(len=32), allocatable :: str_buffer(:)
        
        ierr = 0
        species_idx = 0
        
        call h5dopen_f(file_id, trim(dset_name), dset_id, hdferr)
        if (hdferr /= 0) then
            ierr = -1
            return
        end if
        
        ! Get datatype and size
        call h5dget_type_f(dset_id, dtype_id, hdferr)
        call h5tget_size_f(dtype_id, str_size, hdferr)
        
        ! Get dataspace
        call h5dget_space_f(dset_id, dspace_id, hdferr)
        dims(1) = int(n_strings, hsize_t)
        
        ! Create memory type for fixed-length strings
        call h5tcopy_f(H5T_FORTRAN_S1, memtype_id, hdferr)
        call h5tset_size_f(memtype_id, 32_size_t, hdferr)
        
        ! Read strings
        allocate(str_buffer(n_strings))
        str_buffer = ''
        call h5dread_f(dset_id, memtype_id, str_buffer, dims, hdferr)
        
        if (hdferr /= 0) then
            print '(A)', '  WARNING: Failed to read species strings - check HDF5 file format'
            print '(A)', '  Run convert_hdf5_strings.py to convert variable-length strings'
            ierr = -2
        else
            ! Convert strings to species indices
            do i = 1, n_strings
                species_idx(i) = species_string_to_index(trim(adjustl(str_buffer(i))))
            end do
        end if
        
        call h5tclose_f(memtype_id, hdferr)
        call h5tclose_f(dtype_id, hdferr)
        call h5sclose_f(dspace_id, hdferr)
        call h5dclose_f(dset_id, hdferr)
        deallocate(str_buffer)
        
    end subroutine read_species_strings_fixed


    subroutine read_3d_dataset_transposed(file_id, dset_name, data_out, nx, ny, nz, ierr)
        !----------------------------------------------------------------------
        ! Read 3D dataset from HDF5 file and transpose for Fortran column-major.
        ! HDF5/Python stores as (nz, ny, nx), Fortran needs (nx, ny, nz).
        !----------------------------------------------------------------------
        integer(hid_t), intent(in) :: file_id
        character(len=*), intent(in) :: dset_name
        real(real64), allocatable, intent(out) :: data_out(:,:,:)
        integer(int32), intent(out) :: nx, ny, nz
        integer(int32), intent(out) :: ierr
        
        integer(hid_t) :: dset_id, dspace_id
        integer(hsize_t) :: dims(3), maxdims(3)
        integer :: hdferr, rank
        real(real64), allocatable :: temp(:,:,:)
        integer(int32) :: i, j, k
        
        ierr = 0
        
        call h5dopen_f(file_id, trim(dset_name), dset_id, hdferr)
        if (hdferr /= 0) then
            ierr = -1
            return
        end if
        
        call h5dget_space_f(dset_id, dspace_id, hdferr)
        call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, rank)
        
        ! HDF5 dimensions are reversed from Python: (dim3, dim2, dim1) -> (nx, ny, nz)
        nz = int(dims(1), int32)
        ny = int(dims(2), int32)
        nx = int(dims(3), int32)
        
        allocate(temp(dims(1), dims(2), dims(3)))
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, temp, dims, hdferr)
        
        call h5sclose_f(dspace_id, hdferr)
        call h5dclose_f(dset_id, hdferr)
        
        ! Transpose to Fortran order (nx, ny, nz)
        allocate(data_out(nx, ny, nz))
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    data_out(i, j, k) = temp(k, j, i)
                end do
            end do
        end do
        
        deallocate(temp)
        
    end subroutine read_3d_dataset_transposed

end module ipt_emission_io


module ipt_emission_raytracer
    !--------------------------------------------------------------------------
    ! Main raytracer module for UV/optical emission calculation.
    !--------------------------------------------------------------------------
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic
    use hdf5
    use ipt_emission_constants
    use ipt_emission_utils
    use ipt_emission_types
    use ipt_emission_io
    implicit none
    
    type :: jovian_uv_raytracer
        ! State flags
        logical :: plasma_loaded = .false.
        logical :: single_maxwellian_loaded = .false.
        logical :: double_maxwellian_loaded = .false.
        
        ! Plasma model grid
        integer(int32) :: nx = 0, ny = 0, nz = 0
        real(real64), allocatable :: x_axis(:), y_axis(:), z_axis(:)
        
        ! Plasma density fields [cm^-3]
        real(real64), allocatable :: nec(:,:,:)      ! Cold electron density
        real(real64), allocatable :: neh(:,:,:)      ! Hot electron density
        real(real64), allocatable :: ne_total(:,:,:) ! Total electron density
        real(real64), allocatable :: feh(:,:,:)      ! Hot electron fraction
        real(real64), allocatable :: n_sp(:,:,:)     ! S+ density
        real(real64), allocatable :: n_s2p(:,:,:)    ! S++ density
        real(real64), allocatable :: n_s3p(:,:,:)    ! S+++ density
        real(real64), allocatable :: n_op(:,:,:)     ! O+ density
        real(real64), allocatable :: n_o2p(:,:,:)    ! O++ density
        
        ! Temperature fields [eV]
        real(real64), allocatable :: Tec(:,:,:)      ! Cold electron temperature
        real(real64), allocatable :: Teh(:,:,:)      ! Hot electron temperature
        
        ! Single Maxwellian table grids
        integer(int32) :: n_T_single = 0, n_n_single = 0
        real(real64), allocatable :: temp_arr(:), dens_arr(:)
        real(real64), allocatable :: log_temp(:), log_dens(:)
        
        ! Double Maxwellian table grids
        integer(int32) :: n_Tc = 0, n_Th = 0, n_ne = 0, n_fh = 0
        real(real64), allocatable :: tec_arr(:), teh_arr(:), ne_arr(:), feh_arr(:)
        real(real64), allocatable :: log_tec(:), log_teh(:), log_ne(:), log_feh(:)
        
        ! Per-species emission data
        type(species_emission_data) :: species(N_SPECIES)
        
    contains
        procedure :: initialize
        procedure :: load_plasma_model
        procedure :: load_emission_tables_single
        procedure :: load_emission_tables_double
        procedure :: trace_ray
        procedure :: interp_plasma_trilinear
        procedure :: interp_emission_single_vectorized
        procedure :: interp_emission_double_vectorized
        procedure :: integrate_species_single
        procedure :: integrate_species_double
        procedure :: convolve_spectrum_erf
        procedure :: calculate_spectrum_single
        procedure :: calculate_spectrum_double
        procedure :: cleanup
    end type jovian_uv_raytracer
    
contains

    subroutine initialize(this, plasma_file, emiss_file_single, emiss_file_double, ierr)
        !----------------------------------------------------------------------
        ! Initialize raytracer with plasma model and emission tables.
        !----------------------------------------------------------------------
        class(jovian_uv_raytracer), intent(inout) :: this
        character(len=*), intent(in), optional :: plasma_file, emiss_file_single, emiss_file_double
        integer(int32), intent(out), optional :: ierr
        
        character(len=512) :: pf, ef_single, ef_double
        integer(int32) :: err
        
        ! print '(A)', 'Initializing Jovian UV/Optical Emission Raytracer...'
        
        ! Set default paths
        if (present(plasma_file)) then
            pf = plasma_file
        else
            pf = '../../3D_Torus_Model/jovian_plasma_interpolated_381x381x231.h5'
        end if
        
        if (present(emiss_file_single)) then
            ef_single = emiss_file_single
        else
            ef_single = '../../Emiss_tables/CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5'
        end if
        
        if (present(emiss_file_double)) then
            ef_double = emiss_file_double
        else
            ef_double = '../../Emiss_tables/CHIANTI_11.0.2_emiss_tables_double_maxwellian_16x8x12x10.h5'
        end if
        
        err = 0
        
        call this%load_plasma_model(trim(pf), err)
        if (err /= 0) then
            print '(A)', 'ERROR: Failed to load plasma model'
            if (present(ierr)) ierr = err
            return
        end if
        
        call this%load_emission_tables_single(trim(ef_single), err)
        if (err /= 0) then
            print '(A)', 'ERROR: Failed to load single Maxwellian emission tables'
            if (present(ierr)) ierr = err
            return
        end if
        
        call this%load_emission_tables_double(trim(ef_double), err)
        if (err /= 0) then
            print '(A)', 'NOTE: Double Maxwellian tables not available'
        end if
        
        print '(A)', 'Initialization complete.'
        if (present(ierr)) ierr = 0
        
    end subroutine initialize


    subroutine load_plasma_model(this, filename, ierr)
        !----------------------------------------------------------------------
        ! Load 3D plasma model from HDF5 file.
        !----------------------------------------------------------------------
        class(jovian_uv_raytracer), intent(inout) :: this
        character(len=*), intent(in) :: filename
        integer(int32), intent(out) :: ierr
        
        integer(hid_t) :: file_id
        integer :: hdferr
        integer(int32) :: n1, n2, n3
        logical :: file_exists
        
        ierr = 0
        ! print '(A,A)', 'Loading plasma model from ', trim(filename)
        
        inquire(file=trim(filename), exist=file_exists)
        if (.not. file_exists) then
            print '(A)', '  ERROR: Plasma model file not found'
            ierr = -1
            return
        end if
        
        call h5open_f(hdferr)
        call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, hdferr)
        if (hdferr /= 0) then
            print '(A)', '  ERROR: Cannot open plasma model file'
            ierr = -2
            call h5close_f(hdferr)
            return
        end if
        
        ! Read coordinate arrays
        call read_1d_real_dataset(file_id, 'coordinates/x', this%x_axis, this%nx, ierr)
        call read_1d_real_dataset(file_id, 'coordinates/y', this%y_axis, this%ny, ierr)
        call read_1d_real_dataset(file_id, 'coordinates/z', this%z_axis, this%nz, ierr)
        
        ! print '(A,I0,A,I0,A,I0)', '  Grid size: ', this%nx, ' x ', this%ny, ' x ', this%nz
        
        ! Read plasma fields with transposition
        call read_3d_dataset_transposed(file_id, 'data/ne_c', this%nec, n1, n2, n3, ierr)
        call read_3d_dataset_transposed(file_id, 'data/ne_h', this%neh, n1, n2, n3, ierr)
        call read_3d_dataset_transposed(file_id, 'data/nsp', this%n_sp, n1, n2, n3, ierr)
        call read_3d_dataset_transposed(file_id, 'data/ns2p', this%n_s2p, n1, n2, n3, ierr)
        call read_3d_dataset_transposed(file_id, 'data/ns3p', this%n_s3p, n1, n2, n3, ierr)
        call read_3d_dataset_transposed(file_id, 'data/nop', this%n_op, n1, n2, n3, ierr)
        call read_3d_dataset_transposed(file_id, 'data/no2p', this%n_o2p, n1, n2, n3, ierr)
        call read_3d_dataset_transposed(file_id, 'data/Te_c', this%Tec, n1, n2, n3, ierr)
        call read_3d_dataset_transposed(file_id, 'data/Te_h', this%Teh, n1, n2, n3, ierr)
        
        call h5fclose_f(file_id, hdferr)
        call h5close_f(hdferr)
        
        ! Calculate derived quantities
        allocate(this%ne_total(this%nx, this%ny, this%nz))
        allocate(this%feh(this%nx, this%ny, this%nz))
        
        this%ne_total = this%nec + this%neh
        
        where (this%ne_total > 0.0_real64)
            this%feh = this%neh / this%ne_total
        elsewhere
            this%feh = 0.0_real64
        end where
        
        ! Clean up invalid values
        where (.not. ieee_is_finite(this%nec) .or. this%nec < 0.0_real64) this%nec = 0.0_real64
        where (.not. ieee_is_finite(this%neh) .or. this%neh < 0.0_real64) this%neh = 0.0_real64
        where (.not. ieee_is_finite(this%ne_total) .or. this%ne_total < 0.0_real64) this%ne_total = 0.0_real64
        where (.not. ieee_is_finite(this%feh) .or. this%feh < 0.0_real64) this%feh = 0.0_real64
        where (.not. ieee_is_finite(this%n_sp) .or. this%n_sp < 0.0_real64) this%n_sp = 0.0_real64
        where (.not. ieee_is_finite(this%n_s2p) .or. this%n_s2p < 0.0_real64) this%n_s2p = 0.0_real64
        where (.not. ieee_is_finite(this%n_s3p) .or. this%n_s3p < 0.0_real64) this%n_s3p = 0.0_real64
        where (.not. ieee_is_finite(this%n_op) .or. this%n_op < 0.0_real64) this%n_op = 0.0_real64
        where (.not. ieee_is_finite(this%n_o2p) .or. this%n_o2p < 0.0_real64) this%n_o2p = 0.0_real64
        where (.not. ieee_is_finite(this%Tec) .or. this%Tec < 0.0_real64) this%Tec = 0.0_real64
        where (.not. ieee_is_finite(this%Teh) .or. this%Teh < 0.0_real64) this%Teh = 0.0_real64
        
        this%plasma_loaded = .true.
        
        ! print '(A,F7.2,A,F7.2,A)', '  X range: ', this%x_axis(1), ' to ', this%x_axis(this%nx), ' R_J'
        print '(A,F7.2,A,F7.2,A)', '  Y range: ', this%y_axis(1), ' to ', this%y_axis(this%ny), ' R_J'
        print '(A,F7.2,A,F7.2,A)', '  Z range: ', this%z_axis(1), ' to ', this%z_axis(this%nz), ' R_J'
        
    end subroutine load_plasma_model


    subroutine load_emission_tables_single(this, filename, ierr)
        !----------------------------------------------------------------------
        ! Load single Maxwellian emission tables from HDF5 file.
        ! Emission array stored as (n_T, n_n, n_lines) in HDF5/Python.
        ! Reorganized to per-species arrays (n_lines_sp, n_T, n_n) for vectorized interpolation.
        !----------------------------------------------------------------------
        class(jovian_uv_raytracer), intent(inout) :: this
        character(len=*), intent(in) :: filename
        integer(int32), intent(out) :: ierr
        
        integer(hid_t) :: file_id, dset_id, dspace_id
        integer :: hdferr, rank
        integer(int32) :: n_total, i, j, k, sp_idx, line_idx
        integer(hsize_t) :: dims3(3), maxdims3(3)
        logical :: file_exists
        
        real(real64), allocatable :: wavelength_all(:)
        real(real64), allocatable :: emiss_raw(:,:,:)
        integer(int32), allocatable :: species_idx(:)
        integer(int32) :: species_count(N_SPECIES), species_counter(N_SPECIES)
        
        ierr = 0
        
        inquire(file=trim(filename), exist=file_exists)
        if (.not. file_exists) then
            print '(A,A)', '  ERROR: Single Maxwellian file not found: ', trim(filename)
            ierr = -1
            return
        end if
        
        ! print '(A,A)', 'Loading single Maxwellian emission tables from ', trim(filename)
        
        call h5open_f(hdferr)
        call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, hdferr)
        if (hdferr /= 0) then
            ierr = -2
            call h5close_f(hdferr)
            return
        end if
        
        ! Read temperature and density grids
        call read_1d_real_dataset(file_id, 'T', this%temp_arr, this%n_T_single, ierr)
        call read_1d_real_dataset(file_id, 'n', this%dens_arr, this%n_n_single, ierr)
        
        ! Create log-space grids for interpolation
        allocate(this%log_temp(this%n_T_single))
        allocate(this%log_dens(this%n_n_single))
        this%log_temp = log10(this%temp_arr)
        this%log_dens = log10(this%dens_arr)
        
        ! Read wavelengths
        call read_1d_real_dataset(file_id, 'wavelength', wavelength_all, n_total, ierr)
        
        ! print '(A,I0)', '  Total emission lines: ', n_total
        ! print '(A,I0,A,I0,A)', '  Table grid: ', this%n_T_single, ' T x ', this%n_n_single, ' n'
        
        ! Read species strings (requires preprocessing with convert_hdf5_strings.py)
        allocate(species_idx(n_total))
        call read_species_strings_fixed(file_id, 'species', n_total, species_idx, ierr)
        
        if (ierr /= 0) then
            print '(A)', '  ERROR: Cannot read species strings'
            print '(A)', '  Please run: python convert_hdf5_strings.py'
            call h5fclose_f(file_id, hdferr)
            call h5close_f(hdferr)
            deallocate(wavelength_all, species_idx)
            return
        end if
        
        ! Read emission array: HDF5 stores as (n_lines, n_n, n_T) due to C-order
        call h5dopen_f(file_id, 'emiss', dset_id, hdferr)
        call h5dget_space_f(dset_id, dspace_id, hdferr)
        call h5sget_simple_extent_dims_f(dspace_id, dims3, maxdims3, rank)
        
        allocate(emiss_raw(dims3(1), dims3(2), dims3(3)))
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, emiss_raw, dims3, hdferr)
        call h5sclose_f(dspace_id, hdferr)
        call h5dclose_f(dset_id, hdferr)
        
        call h5fclose_f(file_id, hdferr)
        call h5close_f(hdferr)
        
        ! Count lines per species
        species_count = 0
        do i = 1, n_total
            sp_idx = species_idx(i)
            if (sp_idx >= 1 .and. sp_idx <= N_SPECIES) then
                species_count(sp_idx) = species_count(sp_idx) + 1
            end if
        end do
        
        ! print '(A)', '  Lines per species:'
        do j = 1, N_SPECIES
            if (species_count(j) > 0) then
                print '(A,A,A,I6)', '    ', trim(get_species_key(j)), ': ', species_count(j)
            end if
        end do
        
        ! Allocate per-species arrays
        do j = 1, N_SPECIES
            if (species_count(j) == 0) cycle
            
            this%species(j)%n_lines = species_count(j)
            allocate(this%species(j)%wavelengths(species_count(j)))
            allocate(this%species(j)%emiss_single(species_count(j), this%n_T_single, this%n_n_single))
            this%species(j)%emiss_single = 0.0_real64
        end do
        
        ! Reorganize data by species
        ! HDF5 dims: (n_lines, n_n, n_T) in C-order -> Fortran reads as reversed
        ! emiss_raw indices: (dim1=n_lines, dim2=n_n, dim3=n_T) but actually (T, n, lines) in memory
        species_counter = 0
        do i = 1, n_total
            sp_idx = species_idx(i)
            if (sp_idx < 1 .or. sp_idx > N_SPECIES) cycle
            
            species_counter(sp_idx) = species_counter(sp_idx) + 1
            line_idx = species_counter(sp_idx)
            
            this%species(sp_idx)%wavelengths(line_idx) = wavelength_all(i)
            
            ! Correct transposition: emiss_raw(i, k, j) where i=line, k=n, j=T
            ! Target: emiss_single(line, T, n)
            do j = 1, this%n_T_single
                do k = 1, this%n_n_single
                    this%species(sp_idx)%emiss_single(line_idx, j, k) = emiss_raw(i, k, j)
                end do
            end do
        end do
        
        deallocate(wavelength_all, emiss_raw, species_idx)
        this%single_maxwellian_loaded = .true.
        
        ! print '(A)', 'Single Maxwellian tables loaded successfully:'
        ! print '(A,F8.3,A,F8.1,A)', '  Temperature range: ', this%temp_arr(1), ' - ', &
            ! this%temp_arr(this%n_T_single), ' eV'
        ! print '(A,ES8.1,A,ES8.1,A)', '  Density range: ', this%dens_arr(1), ' - ', &
            ! this%dens_arr(this%n_n_single), ' cm^-3'
        
    end subroutine load_emission_tables_single


    subroutine load_emission_tables_double(this, filename, ierr)
        !----------------------------------------------------------------------
        ! Load double Maxwellian emission tables from HDF5 file.
        ! Emission array stored as (n_Tc, n_Th, n_n, n_feh, n_lines) in HDF5/Python.
        ! Reorganized to per-species arrays (n_lines_sp, n_Tc, n_Th, n_ne, n_feh).
        !----------------------------------------------------------------------
        class(jovian_uv_raytracer), intent(inout) :: this
        character(len=*), intent(in) :: filename
        integer(int32), intent(out) :: ierr
        
        integer(hid_t) :: file_id, dset_id, dspace_id
        integer :: hdferr, rank
        integer(int32) :: n_total, i, j, sp_idx, line_idx
        integer(int32) :: i_Tc, i_Th, i_n, i_f
        integer(hsize_t) :: dims5(5), maxdims5(5)
        logical :: file_exists
        
        real(real64), allocatable :: wavelength_all(:)
        real(real64), allocatable :: emiss_raw(:,:,:,:,:)
        integer(int32), allocatable :: species_idx(:)
        integer(int32) :: species_count(N_SPECIES), species_counter(N_SPECIES)
        
        ierr = 0
        
        inquire(file=trim(filename), exist=file_exists)
        if (.not. file_exists) then
            ierr = -1
            return
        end if
        
        ! print '(A,A)', 'Loading double Maxwellian emission tables from ', trim(filename)
        
        call h5open_f(hdferr)
        call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, hdferr)
        if (hdferr /= 0) then
            ierr = -2
            call h5close_f(hdferr)
            return
        end if
        
        ! Read parameter grids
        call read_1d_real_dataset(file_id, 'T_cold', this%tec_arr, this%n_Tc, ierr)
        call read_1d_real_dataset(file_id, 'T_hot', this%teh_arr, this%n_Th, ierr)
        call read_1d_real_dataset(file_id, 'n', this%ne_arr, this%n_ne, ierr)
        call read_1d_real_dataset(file_id, 'feh', this%feh_arr, this%n_fh, ierr)
        
        ! Create log-space grids
        allocate(this%log_tec(this%n_Tc), this%log_teh(this%n_Th))
        allocate(this%log_ne(this%n_ne), this%log_feh(this%n_fh))
        this%log_tec = log10(this%tec_arr)
        this%log_teh = log10(this%teh_arr)
        this%log_ne = log10(this%ne_arr)
        this%log_feh = log10(this%feh_arr)
        
        ! Read wavelengths and species
        call read_1d_real_dataset(file_id, 'wavelength', wavelength_all, n_total, ierr)
        allocate(species_idx(n_total))
        call read_species_strings_fixed(file_id, 'species', n_total, species_idx, ierr)
        
        if (ierr /= 0) then
            print '(A)', '  ERROR: Cannot read species strings for double Maxwellian'
            call h5fclose_f(file_id, hdferr)
            call h5close_f(hdferr)
            deallocate(wavelength_all, species_idx)
            return
        end if
        
        ! Read emission array
        call h5dopen_f(file_id, 'emiss', dset_id, hdferr)
        call h5dget_space_f(dset_id, dspace_id, hdferr)
        call h5sget_simple_extent_dims_f(dspace_id, dims5, maxdims5, rank)
        
        allocate(emiss_raw(dims5(1), dims5(2), dims5(3), dims5(4), dims5(5)))
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, emiss_raw, dims5, hdferr)
        call h5sclose_f(dspace_id, hdferr)
        call h5dclose_f(dset_id, hdferr)
        
        call h5fclose_f(file_id, hdferr)
        call h5close_f(hdferr)
        
        ! Count lines per species
        species_count = 0
        do i = 1, n_total
            sp_idx = species_idx(i)
            if (sp_idx >= 1 .and. sp_idx <= N_SPECIES) then
                species_count(sp_idx) = species_count(sp_idx) + 1
            end if
        end do
        
        ! Allocate per-species arrays
        do j = 1, N_SPECIES
            if (species_count(j) == 0) cycle
            allocate(this%species(j)%emiss_double(species_count(j), &
                this%n_Tc, this%n_Th, this%n_ne, this%n_fh))
            this%species(j)%emiss_double = 0.0_real64
        end do
        
        ! Reorganize data by species with proper transposition
        ! HDF5/Python: (n_Tc, n_Th, n_n, n_feh, n_lines) -> reversed in Fortran memory
        ! emiss_raw(i, i_f, i_n, i_Th, i_Tc) where i=line index
        species_counter = 0
        do i = 1, n_total
            sp_idx = species_idx(i)
            if (sp_idx < 1 .or. sp_idx > N_SPECIES) cycle
            
            species_counter(sp_idx) = species_counter(sp_idx) + 1
            line_idx = species_counter(sp_idx)
            
            do i_Tc = 1, this%n_Tc
                do i_Th = 1, this%n_Th
                    do i_n = 1, this%n_ne
                        do i_f = 1, this%n_fh
                            this%species(sp_idx)%emiss_double(line_idx, i_Tc, i_Th, i_n, i_f) = &
                                emiss_raw(i, i_f, i_n, i_Th, i_Tc)
                        end do
                    end do
                end do
            end do
        end do
        
        deallocate(wavelength_all, emiss_raw, species_idx)
        this%double_maxwellian_loaded = .true.
        
        ! print '(A)', 'Double Maxwellian tables loaded successfully:'
        ! print '(A,F8.3,A,F8.1,A)', '  Core T range:    ', this%tec_arr(1), ' - ', &
            ! this%tec_arr(this%n_Tc), ' eV'
        ! print '(A,F8.1,A,F8.0,A)', '  Hot T range:     ', this%teh_arr(1), ' - ', &
            ! this%teh_arr(this%n_Th), ' eV'
        ! print '(A,I0,A,I0,A,I0,A,I0,A)', '  Grid: ', this%n_Tc, ' Tc x ', this%n_Th, &
            ! ' Th x ', this%n_ne, ' n x ', this%n_fh, ' feh'
        
    end subroutine load_emission_tables_double


    subroutine trace_ray(this, start_pos, direction, ds, max_dist, s_arr, positions, n_pts)
        !----------------------------------------------------------------------
        ! Trace ray through plasma model with uniform step size.
        !----------------------------------------------------------------------
        class(jovian_uv_raytracer), intent(in) :: this
        real(real64), intent(in) :: start_pos(3), direction(3), ds, max_dist
        real(real64), allocatable, intent(out) :: s_arr(:), positions(:,:)
        integer(int32), intent(out) :: n_pts
        
        real(real64) :: dir_hat(3), norm
        integer(int32) :: i
        
        ! Normalize direction vector
        norm = sqrt(direction(1)**2 + direction(2)**2 + direction(3)**2)
        if (norm > 0.0_real64) then
            dir_hat = direction / norm
        else
            dir_hat = [0.0_real64, 1.0_real64, 0.0_real64]
        end if
        
        ! Create ray points matching Python's np.linspace behavior
        n_pts = int(max_dist / ds) + 1
        allocate(s_arr(n_pts), positions(n_pts, 3))
        
        do i = 1, n_pts
            s_arr(i) = real(i - 1, real64) * max_dist / real(n_pts - 1, real64)
            positions(i, 1) = start_pos(1) + s_arr(i) * dir_hat(1)
            positions(i, 2) = start_pos(2) + s_arr(i) * dir_hat(2)
            positions(i, 3) = start_pos(3) + s_arr(i) * dir_hat(3)
        end do
        
    end subroutine trace_ray


    pure function interp_plasma_trilinear(this, x, y, z, field) result(val)
        !----------------------------------------------------------------------
        ! Trilinear interpolation in 3D plasma field.
        ! Returns 0 for points outside grid (fill_value=0 matches Python).
        !----------------------------------------------------------------------
        class(jovian_uv_raytracer), intent(in) :: this
        real(real64), intent(in) :: x, y, z
        real(real64), intent(in) :: field(:,:,:)
        real(real64) :: val
        
        integer(int32) :: ix, iy, iz
        real(real64) :: wx, wy, wz, w1x, w1y, w1z
        real(real64) :: c000, c001, c010, c011, c100, c101, c110, c111
        
        val = 0.0_real64
        
        ! Return 0 if outside bounds
        if (x < this%x_axis(1) .or. x > this%x_axis(this%nx) .or. &
            y < this%y_axis(1) .or. y > this%y_axis(this%ny) .or. &
            z < this%z_axis(1) .or. z > this%z_axis(this%nz)) return
        
        ! Find bracketing indices
        ix = searchsorted(this%x_axis, x, this%nx)
        iy = searchsorted(this%y_axis, y, this%ny)
        iz = searchsorted(this%z_axis, z, this%nz)
        
        ix = max(2, min(ix, this%nx))
        iy = max(2, min(iy, this%ny))
        iz = max(2, min(iz, this%nz))
        
        ! Compute interpolation weights
        if (this%x_axis(ix) /= this%x_axis(ix-1)) then
            wx = (x - this%x_axis(ix-1)) / (this%x_axis(ix) - this%x_axis(ix-1))
        else
            wx = 0.0_real64
        end if
        
        if (this%y_axis(iy) /= this%y_axis(iy-1)) then
            wy = (y - this%y_axis(iy-1)) / (this%y_axis(iy) - this%y_axis(iy-1))
        else
            wy = 0.0_real64
        end if
        
        if (this%z_axis(iz) /= this%z_axis(iz-1)) then
            wz = (z - this%z_axis(iz-1)) / (this%z_axis(iz) - this%z_axis(iz-1))
        else
            wz = 0.0_real64
        end if
        
        w1x = 1.0_real64 - wx
        w1y = 1.0_real64 - wy
        w1z = 1.0_real64 - wz
        
        ! Get corner values
        c000 = field(ix-1, iy-1, iz-1)
        c001 = field(ix-1, iy-1, iz)
        c010 = field(ix-1, iy,   iz-1)
        c011 = field(ix-1, iy,   iz)
        c100 = field(ix,   iy-1, iz-1)
        c101 = field(ix,   iy-1, iz)
        c110 = field(ix,   iy,   iz-1)
        c111 = field(ix,   iy,   iz)
        
        ! Trilinear interpolation
        val = c000*w1x*w1y*w1z + c001*w1x*w1y*wz + c010*w1x*wy*w1z + c011*w1x*wy*wz + &
              c100*wx*w1y*w1z  + c101*wx*w1y*wz  + c110*wx*wy*w1z  + c111*wx*wy*wz
        
    end function interp_plasma_trilinear


    subroutine interp_emission_single_vectorized(this, Te_los, ne_los, n_los, sp_idx, &
            min_wav, max_wav, wav_out, emiss_out, n_lines)
        !----------------------------------------------------------------------
        ! Vectorized bilinear interpolation of emission rates in log10(T)-log10(n) space.
        ! Interpolates all emission lines for a species at all LOS points simultaneously.
        !----------------------------------------------------------------------
        class(jovian_uv_raytracer), intent(in) :: this
        real(real64), intent(in) :: Te_los(:), ne_los(:)
        integer(int32), intent(in) :: n_los, sp_idx
        real(real64), intent(in) :: min_wav, max_wav
        real(real64), allocatable, intent(out) :: wav_out(:)
        real(real64), allocatable, intent(out) :: emiss_out(:,:)
        integer(int32), intent(out) :: n_lines
        
        integer(int32) :: i, j, k, n_all, n_valid
        integer(int32), allocatable :: line_map(:), valid_idx(:)
        integer(int32), allocatable :: i_T0(:), i_T1(:), i_n0(:), i_n1(:)
        real(real64), allocatable :: log_T(:), log_n(:)
        real(real64), allocatable :: w_T(:), w_n(:), w00(:), w01(:), w10(:), w11(:)
        real(real64) :: dT, dn
        logical, allocatable :: in_bounds(:)
        
        n_lines = 0
        
        if (sp_idx < 1 .or. sp_idx > N_SPECIES) return
        if (this%species(sp_idx)%n_lines == 0) return
        
        n_all = this%species(sp_idx)%n_lines
        allocate(line_map(n_all))
        
        ! Filter wavelengths to range
        do j = 1, n_all
            if (this%species(sp_idx)%wavelengths(j) >= min_wav .and. &
                this%species(sp_idx)%wavelengths(j) <= max_wav) then
                n_lines = n_lines + 1
                line_map(n_lines) = j
            end if
        end do
        
        if (n_lines == 0) then
            deallocate(line_map)
            return
        end if
        
        allocate(wav_out(n_lines))
        allocate(emiss_out(n_los, n_lines))
        emiss_out = 0.0_real64
        
        do k = 1, n_lines
            wav_out(k) = this%species(sp_idx)%wavelengths(line_map(k))
        end do
        
        ! Find valid points (positive, finite, within bounds)
        allocate(in_bounds(n_los))
        in_bounds = (Te_los > 0.0_real64) .and. (ne_los > 0.0_real64) .and. &
                    (Te_los >= this%temp_arr(1)) .and. (Te_los <= this%temp_arr(this%n_T_single)) .and. &
                    (ne_los >= this%dens_arr(1)) .and. (ne_los <= this%dens_arr(this%n_n_single))
        
        n_valid = count(in_bounds)
        if (n_valid == 0) then
            deallocate(line_map, in_bounds)
            return
        end if
        
        allocate(valid_idx(n_valid))
        j = 0
        do i = 1, n_los
            if (in_bounds(i)) then
                j = j + 1
                valid_idx(j) = i
            end if
        end do
        
        ! Convert to log space
        allocate(log_T(n_valid), log_n(n_valid))
        do j = 1, n_valid
            log_T(j) = log10(Te_los(valid_idx(j)))
            log_n(j) = log10(ne_los(valid_idx(j)))
        end do
        
        ! Find bracketing indices
        allocate(i_T0(n_valid), i_T1(n_valid), i_n0(n_valid), i_n1(n_valid))
        do j = 1, n_valid
            i_T1(j) = searchsorted(this%log_temp, log_T(j), this%n_T_single)
            i_T1(j) = max(2, min(i_T1(j), this%n_T_single))
            i_T0(j) = i_T1(j) - 1
            
            i_n1(j) = searchsorted(this%log_dens, log_n(j), this%n_n_single)
            i_n1(j) = max(2, min(i_n1(j), this%n_n_single))
            i_n0(j) = i_n1(j) - 1
        end do
        
        ! Compute interpolation weights
        allocate(w_T(n_valid), w_n(n_valid))
        allocate(w00(n_valid), w01(n_valid), w10(n_valid), w11(n_valid))
        
        do j = 1, n_valid
            dT = this%log_temp(i_T1(j)) - this%log_temp(i_T0(j))
            dn = this%log_dens(i_n1(j)) - this%log_dens(i_n0(j))
            
            if (dT /= 0.0_real64) then
                w_T(j) = (log_T(j) - this%log_temp(i_T0(j))) / dT
            else
                w_T(j) = 0.0_real64
            end if
            
            if (dn /= 0.0_real64) then
                w_n(j) = (log_n(j) - this%log_dens(i_n0(j))) / dn
            else
                w_n(j) = 0.0_real64
            end if
            
            w00(j) = (1.0_real64 - w_T(j)) * (1.0_real64 - w_n(j))
            w01(j) = (1.0_real64 - w_T(j)) * w_n(j)
            w10(j) = w_T(j) * (1.0_real64 - w_n(j))
            w11(j) = w_T(j) * w_n(j)
        end do
        
        ! Vectorized bilinear interpolation for all lines and valid points
        do k = 1, n_lines
            do j = 1, n_valid
                i = valid_idx(j)
                emiss_out(i, k) = &
                    w00(j) * this%species(sp_idx)%emiss_single(line_map(k), i_T0(j), i_n0(j)) + &
                    w01(j) * this%species(sp_idx)%emiss_single(line_map(k), i_T0(j), i_n1(j)) + &
                    w10(j) * this%species(sp_idx)%emiss_single(line_map(k), i_T1(j), i_n0(j)) + &
                    w11(j) * this%species(sp_idx)%emiss_single(line_map(k), i_T1(j), i_n1(j))
            end do
        end do
        
        deallocate(line_map, in_bounds, valid_idx)
        deallocate(log_T, log_n, i_T0, i_T1, i_n0, i_n1)
        deallocate(w_T, w_n, w00, w01, w10, w11)
        
    end subroutine interp_emission_single_vectorized


    subroutine interp_emission_double_vectorized(this, Tec_los, Teh_los, ne_los, feh_los, &
            n_los, sp_idx, min_wav, max_wav, wav_out, emiss_out, n_lines)
        !----------------------------------------------------------------------
        ! Vectorized quadrilinear interpolation over 16 hypercube corners.
        ! Falls back to single Maxwellian when hot fraction or temperature is too low.
        !----------------------------------------------------------------------
        class(jovian_uv_raytracer), intent(in) :: this
        real(real64), intent(in) :: Tec_los(:), Teh_los(:), ne_los(:), feh_los(:)
        integer(int32), intent(in) :: n_los, sp_idx
        real(real64), intent(in) :: min_wav, max_wav
        real(real64), allocatable, intent(out) :: wav_out(:)
        real(real64), allocatable, intent(out) :: emiss_out(:,:)
        integer(int32), intent(out) :: n_lines
        
        integer(int32) :: i, j, k, m, n_all, n_valid
        integer(int32) :: b_Tc, b_Th, b_n, b_f
        integer(int32) :: idx_Tc, idx_Th, idx_n, idx_f
        integer(int32), allocatable :: line_map(:), valid_idx(:)
        integer(int32), allocatable :: i0_Tc(:), i1_Tc(:), i0_Th(:), i1_Th(:)
        integer(int32), allocatable :: i0_n(:), i1_n(:), i0_f(:), i1_f(:)
        real(real64), allocatable :: log_Tc(:), log_Th(:), log_n(:), log_f(:)
        real(real64), allocatable :: w_Tc(:), w_Th(:), w_n(:), w_f(:)
        real(real64) :: dTc, dTh, dn, df, weight, wt_Tc, wt_Th, wt_n, wt_f
        logical, allocatable :: use_double(:), in_bounds(:)
        real(real64), allocatable :: wav_single(:), emiss_single(:,:)
        integer(int32) :: n_lines_single
        
        n_lines = 0
        
        ! Fall back to single Maxwellian if double not loaded
        if (.not. this%double_maxwellian_loaded) then
            call this%interp_emission_single_vectorized(Tec_los, ne_los, n_los, sp_idx, &
                min_wav, max_wav, wav_out, emiss_out, n_lines)
            return
        end if
        
        if (sp_idx < 1 .or. sp_idx > N_SPECIES) return
        if (.not. allocated(this%species(sp_idx)%emiss_double)) return
        
        n_all = this%species(sp_idx)%n_lines
        allocate(line_map(n_all))
        
        ! Filter wavelengths to range
        do j = 1, n_all
            if (this%species(sp_idx)%wavelengths(j) >= min_wav .and. &
                this%species(sp_idx)%wavelengths(j) <= max_wav) then
                n_lines = n_lines + 1
                line_map(n_lines) = j
            end if
        end do
        
        if (n_lines == 0) then
            deallocate(line_map)
            return
        end if
        
        allocate(wav_out(n_lines))
        allocate(emiss_out(n_los, n_lines))
        emiss_out = 0.0_real64
        
        do k = 1, n_lines
            wav_out(k) = this%species(sp_idx)%wavelengths(line_map(k))
        end do
        
        ! Determine which points use single vs double Maxwellian
        allocate(use_double(n_los))
        use_double = (feh_los >= this%feh_arr(1)) .and. (Teh_los >= this%teh_arr(1))
        
        ! Handle single Maxwellian fallback points
        if (any(.not. use_double)) then
            call this%interp_emission_single_vectorized(Tec_los, ne_los, n_los, sp_idx, &
                min_wav, max_wav, wav_single, emiss_single, n_lines_single)
            if (n_lines_single > 0) then
                do i = 1, n_los
                    if (.not. use_double(i)) then
                        emiss_out(i, :) = emiss_single(i, :)
                    end if
                end do
            end if
            if (allocated(wav_single)) deallocate(wav_single)
            if (allocated(emiss_single)) deallocate(emiss_single)
        end if
        
        ! Find valid double Maxwellian points
        allocate(in_bounds(n_los))
        in_bounds = use_double .and. &
                    (Tec_los > 0.0_real64) .and. (Teh_los > 0.0_real64) .and. &
                    (ne_los > 0.0_real64) .and. (feh_los >= 0.0_real64) .and. &
                    (Tec_los >= this%tec_arr(1)) .and. (Tec_los <= this%tec_arr(this%n_Tc)) .and. &
                    (Teh_los >= this%teh_arr(1)) .and. (Teh_los <= this%teh_arr(this%n_Th)) .and. &
                    (ne_los >= this%ne_arr(1)) .and. (ne_los <= this%ne_arr(this%n_ne)) .and. &
                    (feh_los >= this%feh_arr(1)) .and. (feh_los <= this%feh_arr(this%n_fh))
        
        n_valid = count(in_bounds)
        if (n_valid == 0) then
            deallocate(line_map, use_double, in_bounds)
            return
        end if
        
        allocate(valid_idx(n_valid))
        j = 0
        do i = 1, n_los
            if (in_bounds(i)) then
                j = j + 1
                valid_idx(j) = i
            end if
        end do
        
        ! Convert to log space
        allocate(log_Tc(n_valid), log_Th(n_valid), log_n(n_valid), log_f(n_valid))
        do j = 1, n_valid
            log_Tc(j) = log10(Tec_los(valid_idx(j)))
            log_Th(j) = log10(Teh_los(valid_idx(j)))
            log_n(j) = log10(ne_los(valid_idx(j)))
            log_f(j) = log10(feh_los(valid_idx(j)))
        end do
        
        ! Find bracketing indices for all 4 dimensions
        allocate(i0_Tc(n_valid), i1_Tc(n_valid), i0_Th(n_valid), i1_Th(n_valid))
        allocate(i0_n(n_valid), i1_n(n_valid), i0_f(n_valid), i1_f(n_valid))
        allocate(w_Tc(n_valid), w_Th(n_valid), w_n(n_valid), w_f(n_valid))
        
        do j = 1, n_valid
            i1_Tc(j) = searchsorted(this%log_tec, log_Tc(j), this%n_Tc)
            i1_Tc(j) = max(2, min(i1_Tc(j), this%n_Tc))
            i0_Tc(j) = i1_Tc(j) - 1
            
            i1_Th(j) = searchsorted(this%log_teh, log_Th(j), this%n_Th)
            i1_Th(j) = max(2, min(i1_Th(j), this%n_Th))
            i0_Th(j) = i1_Th(j) - 1
            
            i1_n(j) = searchsorted(this%log_ne, log_n(j), this%n_ne)
            i1_n(j) = max(2, min(i1_n(j), this%n_ne))
            i0_n(j) = i1_n(j) - 1
            
            i1_f(j) = searchsorted(this%log_feh, log_f(j), this%n_fh)
            i1_f(j) = max(2, min(i1_f(j), this%n_fh))
            i0_f(j) = i1_f(j) - 1
            
            ! Compute weights
            dTc = this%log_tec(i1_Tc(j)) - this%log_tec(i0_Tc(j))
            dTh = this%log_teh(i1_Th(j)) - this%log_teh(i0_Th(j))
            dn = this%log_ne(i1_n(j)) - this%log_ne(i0_n(j))
            df = this%log_feh(i1_f(j)) - this%log_feh(i0_f(j))
            
            if (dTc /= 0.0_real64) then
                w_Tc(j) = (log_Tc(j) - this%log_tec(i0_Tc(j))) / dTc
            else
                w_Tc(j) = 0.0_real64
            end if
            
            if (dTh /= 0.0_real64) then
                w_Th(j) = (log_Th(j) - this%log_teh(i0_Th(j))) / dTh
            else
                w_Th(j) = 0.0_real64
            end if
            
            if (dn /= 0.0_real64) then
                w_n(j) = (log_n(j) - this%log_ne(i0_n(j))) / dn
            else
                w_n(j) = 0.0_real64
            end if
            
            if (df /= 0.0_real64) then
                w_f(j) = (log_f(j) - this%log_feh(i0_f(j))) / df
            else
                w_f(j) = 0.0_real64
            end if
        end do
        
        ! Quadrilinear interpolation: loop over 16 corners of 4D hypercube
        do m = 1, n_lines
            k = line_map(m)
            
            do j = 1, n_valid
                i = valid_idx(j)
                
                do b_Tc = 0, 1
                    do b_Th = 0, 1
                        do b_n = 0, 1
                            do b_f = 0, 1
                                if (b_Tc == 0) then
                                    idx_Tc = i0_Tc(j)
                                    wt_Tc = 1.0_real64 - w_Tc(j)
                                else
                                    idx_Tc = i1_Tc(j)
                                    wt_Tc = w_Tc(j)
                                end if
                                
                                if (b_Th == 0) then
                                    idx_Th = i0_Th(j)
                                    wt_Th = 1.0_real64 - w_Th(j)
                                else
                                    idx_Th = i1_Th(j)
                                    wt_Th = w_Th(j)
                                end if
                                
                                if (b_n == 0) then
                                    idx_n = i0_n(j)
                                    wt_n = 1.0_real64 - w_n(j)
                                else
                                    idx_n = i1_n(j)
                                    wt_n = w_n(j)
                                end if
                                
                                if (b_f == 0) then
                                    idx_f = i0_f(j)
                                    wt_f = 1.0_real64 - w_f(j)
                                else
                                    idx_f = i1_f(j)
                                    wt_f = w_f(j)
                                end if
                                
                                weight = wt_Tc * wt_Th * wt_n * wt_f
                                emiss_out(i, m) = emiss_out(i, m) + &
                                    weight * this%species(sp_idx)%emiss_double(k, idx_Tc, idx_Th, idx_n, idx_f)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        
        deallocate(line_map, use_double, in_bounds, valid_idx)
        deallocate(log_Tc, log_Th, log_n, log_f)
        deallocate(i0_Tc, i1_Tc, i0_Th, i1_Th, i0_n, i1_n, i0_f, i1_f)
        deallocate(w_Tc, w_Th, w_n, w_f)
        
    end subroutine interp_emission_double_vectorized


    subroutine integrate_species_single(this, s_arr, positions, n_pts, sp_idx, &
            Te_los, ne_los, nion_los, min_wav, max_wav, wav_out, bright_out, n_lines)
        !----------------------------------------------------------------------
        ! Integrate emission for one species along line of sight (single Maxwellian).
        ! Uses Simpson's rule integration matching scipy.integrate.simpson.
        !----------------------------------------------------------------------
        class(jovian_uv_raytracer), intent(in) :: this
        real(real64), intent(in) :: s_arr(:), positions(:,:)
        integer(int32), intent(in) :: n_pts, sp_idx
        real(real64), intent(in) :: Te_los(:), ne_los(:), nion_los(:)
        real(real64), intent(in) :: min_wav, max_wav
        real(real64), allocatable, intent(out) :: wav_out(:), bright_out(:)
        integer(int32), intent(out) :: n_lines
        
        real(real64), allocatable :: emiss_rates(:,:), integrand(:)
        integer(int32) :: k
        
        n_lines = 0
        
        ! Get emission rates at all LOS points
        call this%interp_emission_single_vectorized(Te_los, ne_los, n_pts, sp_idx, &
            min_wav, max_wav, wav_out, emiss_rates, n_lines)
        
        if (n_lines == 0) return
        
        allocate(bright_out(n_lines))
        bright_out = 0.0_real64
        
        allocate(integrand(n_pts))
        
        ! Integrate each line along LOS
        do k = 1, n_lines
            ! Integrand = emission_rate * n_ion
            integrand = emiss_rates(:, k) * nion_los
            bright_out(k) = RAYLEIGH_FACTOR * simpson_integrate(integrand, s_arr, n_pts)
        end do
        
        deallocate(emiss_rates, integrand)
        
    end subroutine integrate_species_single


    subroutine integrate_species_double(this, s_arr, positions, n_pts, sp_idx, &
            Tec_los, Teh_los, ne_los, feh_los, nion_los, min_wav, max_wav, &
            wav_out, bright_out, n_lines)
        !----------------------------------------------------------------------
        ! Integrate emission for one species along line of sight (double Maxwellian).
        !----------------------------------------------------------------------
        class(jovian_uv_raytracer), intent(in) :: this
        real(real64), intent(in) :: s_arr(:), positions(:,:)
        integer(int32), intent(in) :: n_pts, sp_idx
        real(real64), intent(in) :: Tec_los(:), Teh_los(:), ne_los(:), feh_los(:), nion_los(:)
        real(real64), intent(in) :: min_wav, max_wav
        real(real64), allocatable, intent(out) :: wav_out(:), bright_out(:)
        integer(int32), intent(out) :: n_lines
        
        real(real64), allocatable :: emiss_rates(:,:), integrand(:)
        integer(int32) :: k
        
        n_lines = 0
        
        ! Get emission rates at all LOS points
        call this%interp_emission_double_vectorized(Tec_los, Teh_los, ne_los, feh_los, &
            n_pts, sp_idx, min_wav, max_wav, wav_out, emiss_rates, n_lines)
        
        if (n_lines == 0) return
        
        allocate(bright_out(n_lines))
        bright_out = 0.0_real64
        
        allocate(integrand(n_pts))
        
        ! Integrate each line along LOS
        do k = 1, n_lines
            integrand = emiss_rates(:, k) * nion_los
            bright_out(k) = RAYLEIGH_FACTOR * simpson_integrate(integrand, s_arr, n_pts)
        end do
        
        deallocate(emiss_rates, integrand)
        
    end subroutine integrate_species_double


    subroutine convolve_spectrum_erf(this, wave_grid, n_bins, bin_width, &
            line_wav, line_bright, n_lines, fwhm, spectrum)
        !----------------------------------------------------------------------
        ! ERF-based Gaussian convolution for instrument response.
        ! Exactly integrates Gaussian profile over each wavelength bin.
        !----------------------------------------------------------------------
        class(jovian_uv_raytracer), intent(in) :: this
        real(real64), intent(in) :: wave_grid(:), line_wav(:), line_bright(:)
        integer(int32), intent(in) :: n_bins, n_lines
        real(real64), intent(in) :: bin_width, fwhm
        real(real64), intent(out) :: spectrum(:)
        
        real(real64) :: sigma, sigma_sqrt2
        real(real64) :: bin_lo, bin_hi, erf_lo, erf_hi, contrib
        integer(int32) :: i, j
        
        spectrum = 0.0_real64
        
        if (n_lines == 0) return
        
        ! Convert FWHM to sigma: sigma = FWHM / (2 * sqrt(2 * ln(2)))
        sigma = fwhm * FWHM_TO_SIGMA
        sigma_sqrt2 = sigma * SQRT2
        
        ! Convolve each line with instrument response
        do j = 1, n_lines
            if (line_bright(j) <= 0.0_real64) cycle
            
            do i = 1, n_bins
                bin_lo = wave_grid(i) - bin_width / 2.0_real64
                bin_hi = wave_grid(i) + bin_width / 2.0_real64
                
                ! ERF arguments
                erf_lo = (bin_lo - line_wav(j)) / sigma_sqrt2
                erf_hi = (bin_hi - line_wav(j)) / sigma_sqrt2
                
                ! Contribution: 0.5 * (erf(hi) - erf(lo)) * brightness
                contrib = 0.5_real64 * (erf(erf_hi) - erf(erf_lo)) * line_bright(j)
                spectrum(i) = spectrum(i) + contrib
            end do
        end do
        
        ! Convert to R/Angstrom
        spectrum = spectrum / bin_width
        
    end subroutine convolve_spectrum_erf


    subroutine calculate_spectrum_single(this, slit_pos, norm_vec, wave_range, &
            bin_width, fwhm, ds, wave_out, spectrum, lines, n_bins, n_lines_total)
        !----------------------------------------------------------------------
        ! Calculate complete LOS-integrated spectrum for single Maxwellian.
        !----------------------------------------------------------------------
        class(jovian_uv_raytracer), intent(in) :: this
        real(real64), intent(in) :: slit_pos(3), norm_vec(3), wave_range(2)
        real(real64), intent(in) :: bin_width, fwhm, ds
        real(real64), allocatable, intent(out) :: wave_out(:), spectrum(:)
        type(emission_line_result), allocatable, intent(out) :: lines(:)
        integer(int32), intent(out) :: n_bins, n_lines_total
        
        real(real64) :: min_wav, max_wav, max_dist
        integer(int32) :: n_pts, i, j, sp_idx, n_lines_sp
        
        real(real64), allocatable :: s_arr(:), positions(:,:)
        real(real64), allocatable :: ne_los(:), Te_los(:), nion_los(:)
        real(real64), allocatable :: wav_sp(:), bright_sp(:)
        real(real64), allocatable :: all_wav(:), all_bright(:)
        integer(int32), allocatable :: all_species(:)
        integer(int32) :: n_all
        
        ! print '(A)', 'Calculating single Maxwellian spectrum:'
        ! print '(A,3F8.1,A)', '  LOS start: [', slit_pos, '] R_J'
        ! print '(A,3F6.2,A)', '  Direction: [', norm_vec, ']'
        
        min_wav = wave_range(1)
        max_wav = wave_range(2)
        max_dist = 40.0_real64
        
        ! Create wavelength output grid
        n_bins = int((max_wav - min_wav) / bin_width) + 1
        allocate(wave_out(n_bins), spectrum(n_bins))
        do i = 1, n_bins
            wave_out(i) = min_wav + real(i - 1, real64) * (max_wav - min_wav) / real(n_bins - 1, real64)
        end do
        
        ! Trace ray through plasma model
        call this%trace_ray(slit_pos, norm_vec, ds, max_dist, s_arr, positions, n_pts)
        
        ! Interpolate plasma parameters along LOS
        allocate(ne_los(n_pts), Te_los(n_pts), nion_los(n_pts))
        
        do i = 1, n_pts
            ne_los(i) = this%interp_plasma_trilinear(positions(i,1), positions(i,2), positions(i,3), this%nec)
            Te_los(i) = this%interp_plasma_trilinear(positions(i,1), positions(i,2), positions(i,3), this%Tec)
        end do
        
        ! Clean up invalid values
        where (.not. ieee_is_finite(ne_los) .or. ne_los < 0.0_real64) ne_los = 0.0_real64
        where (.not. ieee_is_finite(Te_los) .or. Te_los < 0.0_real64) Te_los = 0.0_real64
        
        ! Collect emission from all species
        allocate(all_wav(MAX_LINES), all_bright(MAX_LINES), all_species(MAX_LINES))
        n_all = 0
        
        do sp_idx = 1, N_SPECIES
            if (this%species(sp_idx)%n_lines == 0) cycle
            
            ! Get ion density for this species
            select case(sp_idx)
            case(IDX_SP)
                do i = 1, n_pts
                    nion_los(i) = this%interp_plasma_trilinear(positions(i,1), positions(i,2), positions(i,3), this%n_sp)
                end do
            case(IDX_S2P)
                do i = 1, n_pts
                    nion_los(i) = this%interp_plasma_trilinear(positions(i,1), positions(i,2), positions(i,3), this%n_s2p)
                end do
            case(IDX_S3P)
                do i = 1, n_pts
                    nion_los(i) = this%interp_plasma_trilinear(positions(i,1), positions(i,2), positions(i,3), this%n_s3p)
                end do
            case(IDX_OP)
                do i = 1, n_pts
                    nion_los(i) = this%interp_plasma_trilinear(positions(i,1), positions(i,2), positions(i,3), this%n_op)
                end do
            case(IDX_O2P)
                do i = 1, n_pts
                    nion_los(i) = this%interp_plasma_trilinear(positions(i,1), positions(i,2), positions(i,3), this%n_o2p)
                end do
            end select
            
            where (.not. ieee_is_finite(nion_los) .or. nion_los < 0.0_real64) nion_los = 0.0_real64
            
            ! Skip if no valid plasma
            if (.not. any(ne_los > 0.0_real64 .and. Te_los > 0.0_real64 .and. nion_los > 0.0_real64)) cycle
            
            ! Integrate emission for this species
            call this%integrate_species_single(s_arr, positions, n_pts, sp_idx, &
                Te_los, ne_los, nion_los, min_wav, max_wav, wav_sp, bright_sp, n_lines_sp)
            
            if (n_lines_sp > 0) then
                do j = 1, n_lines_sp
                    n_all = n_all + 1
                    all_wav(n_all) = wav_sp(j)
                    all_bright(n_all) = bright_sp(j)
                    all_species(n_all) = sp_idx
                end do
                deallocate(wav_sp, bright_sp)
            end if
        end do
        
        ! Count lines with non-zero brightness
        n_lines_total = 0
        do i = 1, n_all
            if (all_bright(i) > 1.0e-10_real64) n_lines_total = n_lines_total + 1
        end do
        
        ! Create output line list
        allocate(lines(max(1, n_lines_total)))
        j = 0
        do i = 1, n_all
            if (all_bright(i) > 1.0e-10_real64) then
                j = j + 1
                lines(j)%wavelength = all_wav(i)
                lines(j)%brightness = all_bright(i)
                lines(j)%species_idx = all_species(i)
            end if
        end do
        
        ! print '(A,I0,A,I0,A)', '  Processed ', n_all, ' lines, ', n_lines_total, ' with non-zero brightness'
        
        ! Convolve with instrument response
        call this%convolve_spectrum_erf(wave_out, n_bins, bin_width, &
            all_wav(1:n_all), all_bright(1:n_all), n_all, fwhm, spectrum)
        
        deallocate(s_arr, positions, ne_los, Te_los, nion_los)
        deallocate(all_wav, all_bright, all_species)
        
    end subroutine calculate_spectrum_single


    subroutine calculate_spectrum_double(this, slit_pos, norm_vec, wave_range, &
            bin_width, fwhm, ds, wave_out, spectrum, lines, n_bins, n_lines_total)
        !----------------------------------------------------------------------
        ! Calculate complete LOS-integrated spectrum for double Maxwellian.
        !----------------------------------------------------------------------
        class(jovian_uv_raytracer), intent(in) :: this
        real(real64), intent(in) :: slit_pos(3), norm_vec(3), wave_range(2)
        real(real64), intent(in) :: bin_width, fwhm, ds
        real(real64), allocatable, intent(out) :: wave_out(:), spectrum(:)
        type(emission_line_result), allocatable, intent(out) :: lines(:)
        integer(int32), intent(out) :: n_bins, n_lines_total
        
        real(real64) :: min_wav, max_wav, max_dist
        integer(int32) :: n_pts, i, j, sp_idx, n_lines_sp
        
        real(real64), allocatable :: s_arr(:), positions(:,:)
        real(real64), allocatable :: ne_los(:), Tec_los(:), Teh_los(:), feh_los(:), nion_los(:)
        real(real64), allocatable :: wav_sp(:), bright_sp(:)
        real(real64), allocatable :: all_wav(:), all_bright(:)
        integer(int32), allocatable :: all_species(:)
        integer(int32) :: n_all
        
        if (.not. this%double_maxwellian_loaded) then
            print '(A)', 'ERROR: Double Maxwellian tables not loaded'
            n_bins = 0
            n_lines_total = 0
            return
        end if
        
        ! print '(A)', 'Calculating double Maxwellian spectrum:'
        ! print '(A,3F8.1,A)', '  LOS start: [', slit_pos, '] R_J'
        ! print '(A,3F6.2,A)', '  Direction: [', norm_vec, ']'
        
        min_wav = wave_range(1)
        max_wav = wave_range(2)
        max_dist = 40.0_real64
        
        ! Create wavelength output grid
        n_bins = int((max_wav - min_wav) / bin_width) + 1
        allocate(wave_out(n_bins), spectrum(n_bins))
        do i = 1, n_bins
            wave_out(i) = min_wav + real(i - 1, real64) * (max_wav - min_wav) / real(n_bins - 1, real64)
        end do
        
        ! Trace ray
        call this%trace_ray(slit_pos, norm_vec, ds, max_dist, s_arr, positions, n_pts)
        
        ! Interpolate plasma parameters
        allocate(ne_los(n_pts), Tec_los(n_pts), Teh_los(n_pts), feh_los(n_pts), nion_los(n_pts))
        
        do i = 1, n_pts
            ne_los(i) = this%interp_plasma_trilinear(positions(i,1), positions(i,2), positions(i,3), this%ne_total)
            Tec_los(i) = this%interp_plasma_trilinear(positions(i,1), positions(i,2), positions(i,3), this%Tec)
            Teh_los(i) = this%interp_plasma_trilinear(positions(i,1), positions(i,2), positions(i,3), this%Teh)
            feh_los(i) = this%interp_plasma_trilinear(positions(i,1), positions(i,2), positions(i,3), this%feh)
        end do
        
        ! Clean up invalid values
        where (.not. ieee_is_finite(ne_los) .or. ne_los < 0.0_real64) ne_los = 0.0_real64
        where (.not. ieee_is_finite(Tec_los) .or. Tec_los < 0.0_real64) Tec_los = 0.0_real64
        where (.not. ieee_is_finite(Teh_los) .or. Teh_los < 0.0_real64) Teh_los = 0.0_real64
        where (.not. ieee_is_finite(feh_los) .or. feh_los < 0.0_real64) feh_los = 0.0_real64
        
        ! Collect emission from all species
        allocate(all_wav(MAX_LINES), all_bright(MAX_LINES), all_species(MAX_LINES))
        n_all = 0
        
        do sp_idx = 1, N_SPECIES
            if (this%species(sp_idx)%n_lines == 0) cycle
            if (.not. allocated(this%species(sp_idx)%emiss_double)) cycle
            
            ! Get ion density for this species
            select case(sp_idx)
            case(IDX_SP)
                do i = 1, n_pts
                    nion_los(i) = this%interp_plasma_trilinear(positions(i,1), positions(i,2), positions(i,3), this%n_sp)
                end do
            case(IDX_S2P)
                do i = 1, n_pts
                    nion_los(i) = this%interp_plasma_trilinear(positions(i,1), positions(i,2), positions(i,3), this%n_s2p)
                end do
            case(IDX_S3P)
                do i = 1, n_pts
                    nion_los(i) = this%interp_plasma_trilinear(positions(i,1), positions(i,2), positions(i,3), this%n_s3p)
                end do
            case(IDX_OP)
                do i = 1, n_pts
                    nion_los(i) = this%interp_plasma_trilinear(positions(i,1), positions(i,2), positions(i,3), this%n_op)
                end do
            case(IDX_O2P)
                do i = 1, n_pts
                    nion_los(i) = this%interp_plasma_trilinear(positions(i,1), positions(i,2), positions(i,3), this%n_o2p)
                end do
            end select
            
            where (.not. ieee_is_finite(nion_los) .or. nion_los < 0.0_real64) nion_los = 0.0_real64
            
            if (.not. any(ne_los > 0.0_real64 .and. Tec_los > 0.0_real64 .and. nion_los > 0.0_real64)) cycle
            
            ! Integrate emission for this species
            call this%integrate_species_double(s_arr, positions, n_pts, sp_idx, &
                Tec_los, Teh_los, ne_los, feh_los, nion_los, min_wav, max_wav, &
                wav_sp, bright_sp, n_lines_sp)
            
            if (n_lines_sp > 0) then
                do j = 1, n_lines_sp
                    n_all = n_all + 1
                    all_wav(n_all) = wav_sp(j)
                    all_bright(n_all) = bright_sp(j)
                    all_species(n_all) = sp_idx
                end do
                deallocate(wav_sp, bright_sp)
            end if
        end do
        
        ! Count lines with non-zero brightness
        n_lines_total = 0
        do i = 1, n_all
            if (all_bright(i) > 1.0e-10_real64) n_lines_total = n_lines_total + 1
        end do
        
        ! Create output line list
        allocate(lines(max(1, n_lines_total)))
        j = 0
        do i = 1, n_all
            if (all_bright(i) > 1.0e-10_real64) then
                j = j + 1
                lines(j)%wavelength = all_wav(i)
                lines(j)%brightness = all_bright(i)
                lines(j)%species_idx = all_species(i)
            end if
        end do
        
        ! print '(A,I0,A,I0,A)', '  Processed ', n_all, ' lines, ', n_lines_total, ' with non-zero brightness'
        
        ! Convolve with instrument response
        call this%convolve_spectrum_erf(wave_out, n_bins, bin_width, &
            all_wav(1:n_all), all_bright(1:n_all), n_all, fwhm, spectrum)
        
        deallocate(s_arr, positions, ne_los, Tec_los, Teh_los, feh_los, nion_los)
        deallocate(all_wav, all_bright, all_species)
        
    end subroutine calculate_spectrum_double


    subroutine cleanup(this)
        !----------------------------------------------------------------------
        ! Deallocate all arrays and reset state.
        !----------------------------------------------------------------------
        class(jovian_uv_raytracer), intent(inout) :: this
        integer(int32) :: j
        
        if (allocated(this%x_axis)) deallocate(this%x_axis)
        if (allocated(this%y_axis)) deallocate(this%y_axis)
        if (allocated(this%z_axis)) deallocate(this%z_axis)
        if (allocated(this%nec)) deallocate(this%nec)
        if (allocated(this%neh)) deallocate(this%neh)
        if (allocated(this%ne_total)) deallocate(this%ne_total)
        if (allocated(this%feh)) deallocate(this%feh)
        if (allocated(this%n_sp)) deallocate(this%n_sp)
        if (allocated(this%n_s2p)) deallocate(this%n_s2p)
        if (allocated(this%n_s3p)) deallocate(this%n_s3p)
        if (allocated(this%n_op)) deallocate(this%n_op)
        if (allocated(this%n_o2p)) deallocate(this%n_o2p)
        if (allocated(this%Tec)) deallocate(this%Tec)
        if (allocated(this%Teh)) deallocate(this%Teh)
        if (allocated(this%temp_arr)) deallocate(this%temp_arr)
        if (allocated(this%dens_arr)) deallocate(this%dens_arr)
        if (allocated(this%log_temp)) deallocate(this%log_temp)
        if (allocated(this%log_dens)) deallocate(this%log_dens)
        if (allocated(this%tec_arr)) deallocate(this%tec_arr)
        if (allocated(this%teh_arr)) deallocate(this%teh_arr)
        if (allocated(this%ne_arr)) deallocate(this%ne_arr)
        if (allocated(this%feh_arr)) deallocate(this%feh_arr)
        if (allocated(this%log_tec)) deallocate(this%log_tec)
        if (allocated(this%log_teh)) deallocate(this%log_teh)
        if (allocated(this%log_ne)) deallocate(this%log_ne)
        if (allocated(this%log_feh)) deallocate(this%log_feh)
        
        do j = 1, N_SPECIES
            if (allocated(this%species(j)%wavelengths)) deallocate(this%species(j)%wavelengths)
            if (allocated(this%species(j)%emiss_single)) deallocate(this%species(j)%emiss_single)
            if (allocated(this%species(j)%emiss_double)) deallocate(this%species(j)%emiss_double)
        end do
        
        this%plasma_loaded = .false.
        this%single_maxwellian_loaded = .false.
        this%double_maxwellian_loaded = .false.
        
    end subroutine cleanup

end module ipt_emission_raytracer
