!==============================================================================
! IPT_emiss_MOP_community_code.f90
!
! Io Plasma Torus (IPT) UV/Optical Emission Line-of-Sight Integration Module
! ==========================================================================
!
! This module provides comprehensive functionality for calculating UV and optical 
! emission spectra from the Io Plasma Torus (IPT) using ray tracing through a 3D 
! plasma model with CHIANTI atomic emission tables.
!
! AUTHOR: Edward (Eddie) G. Nerney
! INSTITUTION: Laboratory for Atmospheric and Space Physics, University of Colorado Boulder
! LICENSE: Open source for academic and research use
! VERSION: 2.0
! DATE: November 2025
!==============================================================================

module ipt_emission_constants
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none
    
    real(real64), parameter :: R_J_KM = 71492.0_real64
    real(real64), parameter :: R_J_CM = 7.1492e9_real64
    real(real64), parameter :: RAYLEIGH_FACTOR = 1.0e-6_real64 * R_J_CM
    real(real64), parameter :: PI = 3.14159265358979323846_real64
    real(real64), parameter :: SQRT2 = 1.41421356237309504880_real64
    real(real64), parameter :: LN2 = 0.69314718055994530942_real64
    
    integer(int32), parameter :: MAX_LINES = 50000
    integer(int32), parameter :: MAX_LOS_POINTS = 5000
    integer(int32), parameter :: MAX_WAVE_BINS = 5000
    integer(int32), parameter :: N_SPECIES = 5
    
    integer(int32), parameter :: IDX_SP = 1
    integer(int32), parameter :: IDX_S2P = 2
    integer(int32), parameter :: IDX_S3P = 3
    integer(int32), parameter :: IDX_OP = 4
    integer(int32), parameter :: IDX_O2P = 5
    
end module ipt_emission_constants

module ipt_emission_utils
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use ipt_emission_constants
    implicit none
contains
    
    function simpson_integrate(y, x) result(integral)
        real(real64), intent(in) :: y(:), x(:)
        real(real64) :: integral
        integer(int32) :: n, i
        real(real64) :: h
        
        n = size(y)
        if (n < 2) then
            integral = 0.0_real64
            return
        end if
        if (n == 2) then
            integral = 0.5_real64 * (y(1) + y(2)) * (x(2) - x(1))
            return
        end if
        
        integral = 0.0_real64
        do i = 1, n-2, 2
            h = (x(i+2) - x(i)) / 2.0_real64
            integral = integral + h/3.0_real64 * (y(i) + 4.0_real64*y(i+1) + y(i+2))
        end do
        if (mod(n, 2) == 0) then
            h = x(n) - x(n-1)
            integral = integral + 0.5_real64 * (y(n-1) + y(n)) * h
        end if
    end function simpson_integrate
    
    function binary_search(arr, val) result(idx)
        real(real64), intent(in) :: arr(:)
        real(real64), intent(in) :: val
        integer(int32) :: idx, lo, hi, mid, n
        
        n = size(arr)
        lo = 1; hi = n
        if (val <= arr(1)) then; idx = 1; return; end if
        if (val >= arr(n)) then; idx = n - 1; return; end if
        
        do while (hi - lo > 1)
            mid = (lo + hi) / 2
            if (arr(mid) <= val) then; lo = mid; else; hi = mid; end if
        end do
        idx = lo
    end function binary_search
    
    function trilinear_interp(xp, yp, zp, x_axis, y_axis, z_axis, field) result(val)
        real(real64), intent(in) :: xp, yp, zp
        real(real64), intent(in) :: x_axis(:), y_axis(:), z_axis(:)
        real(real64), intent(in) :: field(:,:,:)
        real(real64) :: val
        integer(int32) :: ix, iy, iz, nx, ny, nz
        real(real64) :: wx, wy, wz, w1x, w1y, w1z
        real(real64) :: c000, c001, c010, c011, c100, c101, c110, c111
        
        nx = size(x_axis); ny = size(y_axis); nz = size(z_axis)
        if (xp < x_axis(1) .or. xp > x_axis(nx) .or. &
            yp < y_axis(1) .or. yp > y_axis(ny) .or. &
            zp < z_axis(1) .or. zp > z_axis(nz)) then
            val = 0.0_real64; return
        end if
        
        ix = binary_search(x_axis, xp); iy = binary_search(y_axis, yp); iz = binary_search(z_axis, zp)
        ix = max(1, min(ix, nx-1)); iy = max(1, min(iy, ny-1)); iz = max(1, min(iz, nz-1))
        
        if (x_axis(ix+1) /= x_axis(ix)) then
            wx = (xp - x_axis(ix)) / (x_axis(ix+1) - x_axis(ix))
        else; wx = 0.0_real64; end if
        if (y_axis(iy+1) /= y_axis(iy)) then
            wy = (yp - y_axis(iy)) / (y_axis(iy+1) - y_axis(iy))
        else; wy = 0.0_real64; end if
        if (z_axis(iz+1) /= z_axis(iz)) then
            wz = (zp - z_axis(iz)) / (z_axis(iz+1) - z_axis(iz))
        else; wz = 0.0_real64; end if
        
        w1x = 1.0_real64 - wx; w1y = 1.0_real64 - wy; w1z = 1.0_real64 - wz
        
        c000 = field(ix, iy, iz);     c001 = field(ix, iy, iz+1)
        c010 = field(ix, iy+1, iz);   c011 = field(ix, iy+1, iz+1)
        c100 = field(ix+1, iy, iz);   c101 = field(ix+1, iy, iz+1)
        c110 = field(ix+1, iy+1, iz); c111 = field(ix+1, iy+1, iz+1)
        
        val = c000*w1x*w1y*w1z + c001*w1x*w1y*wz + c010*w1x*wy*w1z + c011*w1x*wy*wz + &
              c100*wx*w1y*w1z + c101*wx*w1y*wz + c110*wx*wy*w1z + c111*wx*wy*wz
    end function trilinear_interp
    
    function vec_norm(v) result(n)
        real(real64), intent(in) :: v(3)
        real(real64) :: n
        n = sqrt(v(1)**2 + v(2)**2 + v(3)**2)
    end function vec_norm
    
    function get_species_name(idx) result(name)
        integer(int32), intent(in) :: idx
        character(len=8) :: name
        select case(idx)
        case(IDX_SP);  name = 'SP'
        case(IDX_S2P); name = 'S2P'
        case(IDX_S3P); name = 'S3P'
        case(IDX_OP);  name = 'OP'
        case(IDX_O2P); name = 'O2P'
        case default;  name = 'UNKNOWN'
        end select
    end function get_species_name
    
    function get_display_name(idx) result(name)
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
    
    function species_string_to_index(sp_str) result(idx)
        character(len=*), intent(in) :: sp_str
        integer(int32) :: idx
        character(len=8) :: sp_clean
        integer :: i
        
        sp_clean = ''
        do i = 1, min(len_trim(sp_str), 8)
            if (ichar(sp_str(i:i)) >= 32 .and. ichar(sp_str(i:i)) <= 126) then
                sp_clean = trim(sp_clean) // sp_str(i:i)
            end if
        end do
        sp_clean = adjustl(sp_clean)
        
        select case(trim(sp_clean))
        case('SP', 'sp', 'Sp'); idx = IDX_SP
        case('S2P', 's2p', 'S2p'); idx = IDX_S2P
        case('S3P', 's3p', 'S3p'); idx = IDX_S3P
        case('OP', 'op', 'Op'); idx = IDX_OP
        case('O2P', 'o2p', 'O2p'); idx = IDX_O2P
        case default; idx = 0
        end select
    end function species_string_to_index
    
end module ipt_emission_utils

module ipt_emission_types
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use ipt_emission_constants
    implicit none
    
    type :: species_emission_data
        character(len=8) :: name = ''
        integer(int32) :: n_lines = 0
        real(real64), allocatable :: wavelengths(:)
        real(real64), allocatable :: emissivities_single(:,:,:)
        real(real64), allocatable :: emissivities_double(:,:,:,:,:)
    end type species_emission_data
    
    type :: emission_line
        real(real64) :: wavelength = 0.0_real64
        real(real64) :: brightness = 0.0_real64
        integer(int32) :: species_idx = 0
    end type emission_line
    
end module ipt_emission_types

module ipt_emission_raytracer
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic
    use hdf5
    use ipt_emission_constants
    use ipt_emission_utils
    use ipt_emission_types
    implicit none
    
    type :: jovian_uv_raytracer
        logical :: plasma_loaded = .false.
        logical :: single_maxwellian_loaded = .false.
        logical :: double_maxwellian_loaded = .false.
        
        integer(int32) :: nx = 0, ny = 0, nz = 0
        real(real64), allocatable :: x_axis(:), y_axis(:), z_axis(:)
        real(real64), allocatable :: nec(:,:,:), neh(:,:,:), ne_total(:,:,:), feh(:,:,:)
        real(real64), allocatable :: nsp(:,:,:), ns2p(:,:,:), ns3p(:,:,:)
        real(real64), allocatable :: nop(:,:,:), no2p(:,:,:)
        real(real64), allocatable :: Tec(:,:,:), Teh(:,:,:)
        
        integer(int32) :: n_temp_single = 0, n_dens_single = 0
        real(real64), allocatable :: temp_arr(:), dens_arr(:), log_temp(:), log_dens(:)
        
        integer(int32) :: n_tec = 0, n_teh = 0, n_ne = 0, n_feh = 0
        real(real64), allocatable :: tec_arr(:), teh_arr(:), ne_arr(:), feh_arr(:)
        real(real64), allocatable :: log_tec(:), log_teh(:), log_ne(:), log_feh(:)
        
        type(species_emission_data) :: species(N_SPECIES)
        
    contains
        procedure :: initialize
        procedure :: load_plasma_model
        procedure :: load_emission_tables_single
        procedure :: load_emission_tables_double
        procedure :: trace_ray
        procedure :: interpolate_plasma_trilinear
        procedure :: interpolate_emission_single
        procedure :: interpolate_emission_double
        procedure :: integrate_species_single
        procedure :: integrate_species_double
        procedure :: convolve_spectrum_erf
        procedure :: calculate_spectrum_single
        procedure :: calculate_spectrum_double
        procedure :: cleanup
    end type jovian_uv_raytracer
    
contains
    
    subroutine initialize(this, plasma_file, emission_file_single, emission_file_double)
        class(jovian_uv_raytracer), intent(inout) :: this
        character(len=*), intent(in) :: plasma_file, emission_file_single
        character(len=*), intent(in), optional :: emission_file_double
        integer(int32) :: ierr
        
        print '(A)', 'Initializing Jovian UV/Optical Emission Raytracer...'
        call this%load_plasma_model(plasma_file, ierr)
        if (ierr /= 0) then
            print '(A,I0)', 'ERROR: Failed to load plasma model, error code: ', ierr
            return
        end if
        call this%load_emission_tables_single(emission_file_single, ierr)
        if (ierr /= 0) then
            print '(A,I0)', 'ERROR: Failed to load single Maxwellian tables, error code: ', ierr
            return
        end if
        if (present(emission_file_double)) then
            call this%load_emission_tables_double(emission_file_double, ierr)
            if (ierr /= 0) print '(A)', 'WARNING: Double Maxwellian tables not loaded'
        end if
        print '(A)', 'Initialization complete.'
    end subroutine initialize
    
    subroutine load_plasma_model(this, filename, ierr)
        class(jovian_uv_raytracer), intent(inout) :: this
        character(len=*), intent(in) :: filename
        integer(int32), intent(out) :: ierr
        
        integer(hid_t) :: file_id, dset_id, dspace_id
        integer(hsize_t) :: dims1(1), maxdims1(1), dims3(3)
        integer :: hdferr, i, j, k
        real(real64), allocatable :: temp_3d(:,:,:)
        
        ierr = 0
        print '(A,A,A)', 'Loading plasma model from ', trim(filename), '...'
        
        call h5open_f(hdferr)
        if (hdferr /= 0) then; ierr = -1; return; end if
        
        call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, hdferr)
        if (hdferr /= 0) then; ierr = -2; call h5close_f(hdferr); return; end if
        
        ! Read coordinate axes
        call h5dopen_f(file_id, 'coordinates/x', dset_id, hdferr)
        call h5dget_space_f(dset_id, dspace_id, hdferr)
        call h5sget_simple_extent_dims_f(dspace_id, dims1, maxdims1, hdferr)
        this%nx = int(dims1(1), int32)
        allocate(this%x_axis(this%nx))
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, this%x_axis, dims1, hdferr)
        call h5sclose_f(dspace_id, hdferr); call h5dclose_f(dset_id, hdferr)
        
        call h5dopen_f(file_id, 'coordinates/y', dset_id, hdferr)
        call h5dget_space_f(dset_id, dspace_id, hdferr)
        call h5sget_simple_extent_dims_f(dspace_id, dims1, maxdims1, hdferr)
        this%ny = int(dims1(1), int32)
        allocate(this%y_axis(this%ny))
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, this%y_axis, dims1, hdferr)
        call h5sclose_f(dspace_id, hdferr); call h5dclose_f(dset_id, hdferr)
        
        call h5dopen_f(file_id, 'coordinates/z', dset_id, hdferr)
        call h5dget_space_f(dset_id, dspace_id, hdferr)
        call h5sget_simple_extent_dims_f(dspace_id, dims1, maxdims1, hdferr)
        this%nz = int(dims1(1), int32)
        allocate(this%z_axis(this%nz))
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, this%z_axis, dims1, hdferr)
        call h5sclose_f(dspace_id, hdferr); call h5dclose_f(dset_id, hdferr)
        
        ! Allocate 3D fields
        allocate(this%nec(this%nx, this%ny, this%nz), this%neh(this%nx, this%ny, this%nz))
        allocate(this%ne_total(this%nx, this%ny, this%nz), this%feh(this%nx, this%ny, this%nz))
        allocate(this%nsp(this%nx, this%ny, this%nz), this%ns2p(this%nx, this%ny, this%nz))
        allocate(this%ns3p(this%nx, this%ny, this%nz), this%nop(this%nx, this%ny, this%nz))
        allocate(this%no2p(this%nx, this%ny, this%nz))
        allocate(this%Tec(this%nx, this%ny, this%nz), this%Teh(this%nx, this%ny, this%nz))
        allocate(temp_3d(this%nz, this%ny, this%nx))
        
        dims3 = [int(this%nz, hsize_t), int(this%ny, hsize_t), int(this%nx, hsize_t)]
        
        ! Read and transpose each field (HDF5 C-order to Fortran column-major)
        call h5dopen_f(file_id, 'data/ne_c', dset_id, hdferr)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, temp_3d, dims3, hdferr)
        call h5dclose_f(dset_id, hdferr)
        do k = 1, this%nz; do j = 1, this%ny; do i = 1, this%nx
            this%nec(i,j,k) = temp_3d(k,j,i)
        end do; end do; end do
        
        call h5dopen_f(file_id, 'data/ne_h', dset_id, hdferr)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, temp_3d, dims3, hdferr)
        call h5dclose_f(dset_id, hdferr)
        do k = 1, this%nz; do j = 1, this%ny; do i = 1, this%nx
            this%neh(i,j,k) = temp_3d(k,j,i)
        end do; end do; end do
        
        call h5dopen_f(file_id, 'data/nsp', dset_id, hdferr)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, temp_3d, dims3, hdferr)
        call h5dclose_f(dset_id, hdferr)
        do k = 1, this%nz; do j = 1, this%ny; do i = 1, this%nx
            this%nsp(i,j,k) = temp_3d(k,j,i)
        end do; end do; end do
        
        call h5dopen_f(file_id, 'data/ns2p', dset_id, hdferr)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, temp_3d, dims3, hdferr)
        call h5dclose_f(dset_id, hdferr)
        do k = 1, this%nz; do j = 1, this%ny; do i = 1, this%nx
            this%ns2p(i,j,k) = temp_3d(k,j,i)
        end do; end do; end do
        
        call h5dopen_f(file_id, 'data/ns3p', dset_id, hdferr)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, temp_3d, dims3, hdferr)
        call h5dclose_f(dset_id, hdferr)
        do k = 1, this%nz; do j = 1, this%ny; do i = 1, this%nx
            this%ns3p(i,j,k) = temp_3d(k,j,i)
        end do; end do; end do
        
        call h5dopen_f(file_id, 'data/nop', dset_id, hdferr)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, temp_3d, dims3, hdferr)
        call h5dclose_f(dset_id, hdferr)
        do k = 1, this%nz; do j = 1, this%ny; do i = 1, this%nx
            this%nop(i,j,k) = temp_3d(k,j,i)
        end do; end do; end do
        
        call h5dopen_f(file_id, 'data/no2p', dset_id, hdferr)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, temp_3d, dims3, hdferr)
        call h5dclose_f(dset_id, hdferr)
        do k = 1, this%nz; do j = 1, this%ny; do i = 1, this%nx
            this%no2p(i,j,k) = temp_3d(k,j,i)
        end do; end do; end do
        
        call h5dopen_f(file_id, 'data/Te_c', dset_id, hdferr)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, temp_3d, dims3, hdferr)
        call h5dclose_f(dset_id, hdferr)
        do k = 1, this%nz; do j = 1, this%ny; do i = 1, this%nx
            this%Tec(i,j,k) = temp_3d(k,j,i)
        end do; end do; end do
        
        call h5dopen_f(file_id, 'data/Te_h', dset_id, hdferr)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, temp_3d, dims3, hdferr)
        call h5dclose_f(dset_id, hdferr)
        do k = 1, this%nz; do j = 1, this%ny; do i = 1, this%nx
            this%Teh(i,j,k) = temp_3d(k,j,i)
        end do; end do; end do
        
        deallocate(temp_3d)
        call h5fclose_f(file_id, hdferr)
        call h5close_f(hdferr)
        
        ! Calculate derived quantities
        this%ne_total = this%nec + this%neh
        do k = 1, this%nz; do j = 1, this%ny; do i = 1, this%nx
            if (this%ne_total(i,j,k) > 0.0_real64) then
                this%feh(i,j,k) = this%neh(i,j,k) / this%ne_total(i,j,k)
            else
                this%feh(i,j,k) = 0.0_real64
            end if
        end do; end do; end do
        
        ! Clean invalid values
        where (.not. ieee_is_finite(this%nec) .or. this%nec < 0.0_real64) this%nec = 0.0_real64
        where (.not. ieee_is_finite(this%neh) .or. this%neh < 0.0_real64) this%neh = 0.0_real64
        where (.not. ieee_is_finite(this%nsp) .or. this%nsp < 0.0_real64) this%nsp = 0.0_real64
        where (.not. ieee_is_finite(this%ns2p) .or. this%ns2p < 0.0_real64) this%ns2p = 0.0_real64
        where (.not. ieee_is_finite(this%ns3p) .or. this%ns3p < 0.0_real64) this%ns3p = 0.0_real64
        where (.not. ieee_is_finite(this%nop) .or. this%nop < 0.0_real64) this%nop = 0.0_real64
        where (.not. ieee_is_finite(this%no2p) .or. this%no2p < 0.0_real64) this%no2p = 0.0_real64
        where (.not. ieee_is_finite(this%Tec) .or. this%Tec < 0.0_real64) this%Tec = 0.0_real64
        where (.not. ieee_is_finite(this%Teh) .or. this%Teh < 0.0_real64) this%Teh = 0.0_real64
        where (.not. ieee_is_finite(this%feh) .or. this%feh < 0.0_real64) this%feh = 0.0_real64
        where (.not. ieee_is_finite(this%ne_total) .or. this%ne_total < 0.0_real64) this%ne_total = 0.0_real64
        
        this%plasma_loaded = .true.
        
        print '(A,I0,A,I0,A,I0)', 'Plasma model loaded: grid shape ', this%nx, ' x ', this%ny, ' x ', this%nz
        print '(A,F6.1,A,F6.1,A)', '  X range: ', this%x_axis(1), ' to ', this%x_axis(this%nx), ' R_J'
        print '(A,F6.1,A,F6.1,A)', '  Y range: ', this%y_axis(1), ' to ', this%y_axis(this%ny), ' R_J'
        print '(A,F6.1,A,F6.1,A)', '  Z range: ', this%z_axis(1), ' to ', this%z_axis(this%nz), ' R_J'
        print '(A,ES10.3)', '  Max ne_c: ', maxval(this%nec)
        print '(A,ES10.3)', '  Max Te_c: ', maxval(this%Tec)
    end subroutine load_plasma_model
    
    subroutine load_emission_tables_single(this, filename, ierr)
        class(jovian_uv_raytracer), intent(inout) :: this
        character(len=*), intent(in) :: filename
        integer(int32), intent(out) :: ierr
        
        integer(hid_t) :: file_id, dset_id, dspace_id, type_id, memtype_id
        integer(hsize_t) :: dims1(1), maxdims1(1), dims3(3)
        integer(size_t) :: str_len
        integer :: hdferr, n_total, i, j, k, n_sp, sp_idx
        real(real64), allocatable :: emiss_all(:,:,:), wavelength_all(:)
        character(len=8), allocatable :: species_all(:)
        integer(int32) :: species_counts(N_SPECIES), line_counters(N_SPECIES)
        integer(int32), allocatable :: species_indices(:)
        
        ierr = 0
        print '(A,A,A)', 'Loading single Maxwellian emission tables from ', trim(filename), '...'
        
        call h5open_f(hdferr)
        if (hdferr /= 0) then; ierr = -1; return; end if
        call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, hdferr)
        if (hdferr /= 0) then; ierr = -2; call h5close_f(hdferr); return; end if
        
        ! Read temperature grid
        call h5dopen_f(file_id, 'T', dset_id, hdferr)
        call h5dget_space_f(dset_id, dspace_id, hdferr)
        call h5sget_simple_extent_dims_f(dspace_id, dims1, maxdims1, hdferr)
        this%n_temp_single = int(dims1(1), int32)
        allocate(this%temp_arr(this%n_temp_single), this%log_temp(this%n_temp_single))
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, this%temp_arr, dims1, hdferr)
        call h5sclose_f(dspace_id, hdferr); call h5dclose_f(dset_id, hdferr)
        this%log_temp = log10(this%temp_arr)
        
        ! Read density grid
        call h5dopen_f(file_id, 'n', dset_id, hdferr)
        call h5dget_space_f(dset_id, dspace_id, hdferr)
        call h5sget_simple_extent_dims_f(dspace_id, dims1, maxdims1, hdferr)
        this%n_dens_single = int(dims1(1), int32)
        allocate(this%dens_arr(this%n_dens_single), this%log_dens(this%n_dens_single))
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, this%dens_arr, dims1, hdferr)
        call h5sclose_f(dspace_id, hdferr); call h5dclose_f(dset_id, hdferr)
        this%log_dens = log10(this%dens_arr)
        
        ! Read wavelengths
        call h5dopen_f(file_id, 'wavelength', dset_id, hdferr)
        call h5dget_space_f(dset_id, dspace_id, hdferr)
        call h5sget_simple_extent_dims_f(dspace_id, dims1, maxdims1, hdferr)
        n_total = int(dims1(1), int32)
        allocate(wavelength_all(n_total))
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, wavelength_all, dims1, hdferr)
        call h5sclose_f(dspace_id, hdferr); call h5dclose_f(dset_id, hdferr)
        
        ! Read species labels
        allocate(species_all(n_total), species_indices(n_total))
        species_all = ''
        call h5dopen_f(file_id, 'species', dset_id, hdferr)
        call h5dget_type_f(dset_id, type_id, hdferr)
        call h5tget_size_f(type_id, str_len, hdferr)
        call h5tcopy_f(H5T_FORTRAN_S1, memtype_id, hdferr)
        call h5tset_size_f(memtype_id, int(8, size_t), hdferr)
        call h5dread_f(dset_id, memtype_id, species_all, dims1, hdferr)
        call h5tclose_f(memtype_id, hdferr); call h5tclose_f(type_id, hdferr)
        call h5dclose_f(dset_id, hdferr)
        
        ! Map species strings to indices
        species_counts = 0
        do i = 1, n_total
            sp_idx = species_string_to_index(species_all(i))
            species_indices(i) = sp_idx
            if (sp_idx > 0 .and. sp_idx <= N_SPECIES) species_counts(sp_idx) = species_counts(sp_idx) + 1
        end do
        
        print '(A,I0)', '  Total emission lines: ', n_total
        do j = 1, N_SPECIES
            if (species_counts(j) > 0) print '(A,A,A,I0)', '    ', trim(get_species_name(j)), ': ', species_counts(j)
        end do
        
        ! Read emissivity array - Python shape (n_T, n_ne, n_lines) -> Fortran reads as (n_lines, n_ne, n_T)
        dims3 = [int(n_total, hsize_t), int(this%n_dens_single, hsize_t), int(this%n_temp_single, hsize_t)]
        allocate(emiss_all(n_total, this%n_dens_single, this%n_temp_single))
        call h5dopen_f(file_id, 'emiss', dset_id, hdferr)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, emiss_all, dims3, hdferr)
        call h5dclose_f(dset_id, hdferr)
        
        call h5fclose_f(file_id, hdferr)
        call h5close_f(hdferr)
        
        ! Organize by species - store as (n_lines, n_T, n_ne)
        do j = 1, N_SPECIES
            n_sp = species_counts(j)
            if (n_sp == 0) cycle
            this%species(j)%name = get_species_name(j)
            this%species(j)%n_lines = n_sp
            allocate(this%species(j)%wavelengths(n_sp))
            allocate(this%species(j)%emissivities_single(n_sp, this%n_temp_single, this%n_dens_single))
            this%species(j)%emissivities_single = 0.0_real64
        end do
        
        ! Fill species arrays with proper index reordering
        line_counters = 0
        do i = 1, n_total
            sp_idx = species_indices(i)
            if (sp_idx < 1 .or. sp_idx > N_SPECIES) cycle
            line_counters(sp_idx) = line_counters(sp_idx) + 1
            k = line_counters(sp_idx)
            this%species(sp_idx)%wavelengths(k) = wavelength_all(i)
            ! emiss_all is (line, ne, T), store as (line, T, ne)
            do j = 1, this%n_temp_single
                do n_sp = 1, this%n_dens_single
                    this%species(sp_idx)%emissivities_single(k, j, n_sp) = emiss_all(i, n_sp, j)
                end do
            end do
        end do
        
        deallocate(emiss_all, wavelength_all, species_all, species_indices)
        this%single_maxwellian_loaded = .true.
        
        print '(A)', 'Single Maxwellian tables loaded successfully:'
        print '(A,F8.3,A,F8.1,A)', '  Temperature range: ', this%temp_arr(1), ' - ', this%temp_arr(this%n_temp_single), ' eV'
        print '(A,F8.1,A,F8.0,A)', '  Density range: ', this%dens_arr(1), ' - ', this%dens_arr(this%n_dens_single), ' cm^-3'
        print '(A,I0,A,I0)', '  Grid size: ', this%n_temp_single, ' x ', this%n_dens_single
        write(*, '(A)', advance='no') '  Species: '
        do j = 1, N_SPECIES
            if (this%species(j)%n_lines > 0) write(*, '(A,A,I0,A)', advance='no') trim(this%species(j)%name), '(', this%species(j)%n_lines, ') '
        end do
        print '(A)', ''
    end subroutine load_emission_tables_single
    
    subroutine load_emission_tables_double(this, filename, ierr)
        class(jovian_uv_raytracer), intent(inout) :: this
        character(len=*), intent(in) :: filename
        integer(int32), intent(out) :: ierr
        
        integer(hid_t) :: file_id, dset_id, dspace_id, type_id, memtype_id
        integer(hsize_t) :: dims1(1), maxdims1(1), dims5(5)
        integer(size_t) :: str_len
        integer :: hdferr, n_total, i, j, k, n_sp, sp_idx, i_tc, i_th, i_ne, i_feh
        real(real64), allocatable :: emiss_all(:,:,:,:,:), wavelength_all(:)
        character(len=8), allocatable :: species_all(:)
        integer(int32) :: species_counts(N_SPECIES), line_counters(N_SPECIES)
        integer(int32), allocatable :: species_indices(:)
        logical :: file_exists
        
        ierr = 0
        inquire(file=trim(filename), exist=file_exists)
        if (.not. file_exists) then; ierr = -1; return; end if
        
        print '(A,A,A)', 'Loading double Maxwellian emission tables from ', trim(filename), '...'
        
        call h5open_f(hdferr)
        if (hdferr /= 0) then; ierr = -2; return; end if
        call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, hdferr)
        if (hdferr /= 0) then; ierr = -3; call h5close_f(hdferr); return; end if
        
        ! Read parameter grids
        call h5dopen_f(file_id, 'T_cold', dset_id, hdferr)
        call h5dget_space_f(dset_id, dspace_id, hdferr)
        call h5sget_simple_extent_dims_f(dspace_id, dims1, maxdims1, hdferr)
        this%n_tec = int(dims1(1), int32)
        allocate(this%tec_arr(this%n_tec), this%log_tec(this%n_tec))
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, this%tec_arr, dims1, hdferr)
        call h5sclose_f(dspace_id, hdferr); call h5dclose_f(dset_id, hdferr)
        this%log_tec = log10(this%tec_arr)
        
        call h5dopen_f(file_id, 'T_hot', dset_id, hdferr)
        call h5dget_space_f(dset_id, dspace_id, hdferr)
        call h5sget_simple_extent_dims_f(dspace_id, dims1, maxdims1, hdferr)
        this%n_teh = int(dims1(1), int32)
        allocate(this%teh_arr(this%n_teh), this%log_teh(this%n_teh))
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, this%teh_arr, dims1, hdferr)
        call h5sclose_f(dspace_id, hdferr); call h5dclose_f(dset_id, hdferr)
        this%log_teh = log10(this%teh_arr)
        
        call h5dopen_f(file_id, 'n', dset_id, hdferr)
        call h5dget_space_f(dset_id, dspace_id, hdferr)
        call h5sget_simple_extent_dims_f(dspace_id, dims1, maxdims1, hdferr)
        this%n_ne = int(dims1(1), int32)
        allocate(this%ne_arr(this%n_ne), this%log_ne(this%n_ne))
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, this%ne_arr, dims1, hdferr)
        call h5sclose_f(dspace_id, hdferr); call h5dclose_f(dset_id, hdferr)
        this%log_ne = log10(this%ne_arr)
        
        call h5dopen_f(file_id, 'feh', dset_id, hdferr)
        call h5dget_space_f(dset_id, dspace_id, hdferr)
        call h5sget_simple_extent_dims_f(dspace_id, dims1, maxdims1, hdferr)
        this%n_feh = int(dims1(1), int32)
        allocate(this%feh_arr(this%n_feh), this%log_feh(this%n_feh))
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, this%feh_arr, dims1, hdferr)
        call h5sclose_f(dspace_id, hdferr); call h5dclose_f(dset_id, hdferr)
        this%log_feh = log10(this%feh_arr)
        
        ! Read wavelengths
        call h5dopen_f(file_id, 'wavelength', dset_id, hdferr)
        call h5dget_space_f(dset_id, dspace_id, hdferr)
        call h5sget_simple_extent_dims_f(dspace_id, dims1, maxdims1, hdferr)
        n_total = int(dims1(1), int32)
        allocate(wavelength_all(n_total))
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, wavelength_all, dims1, hdferr)
        call h5sclose_f(dspace_id, hdferr); call h5dclose_f(dset_id, hdferr)
        
        ! Read species labels
        allocate(species_all(n_total), species_indices(n_total))
        species_all = ''
        call h5dopen_f(file_id, 'species', dset_id, hdferr)
        call h5dget_type_f(dset_id, type_id, hdferr)
        call h5tget_size_f(type_id, str_len, hdferr)
        call h5tcopy_f(H5T_FORTRAN_S1, memtype_id, hdferr)
        call h5tset_size_f(memtype_id, int(8, size_t), hdferr)
        call h5dread_f(dset_id, memtype_id, species_all, dims1, hdferr)
        call h5tclose_f(memtype_id, hdferr); call h5tclose_f(type_id, hdferr)
        call h5dclose_f(dset_id, hdferr)
        
        species_counts = 0
        do i = 1, n_total
            sp_idx = species_string_to_index(species_all(i))
            species_indices(i) = sp_idx
            if (sp_idx > 0 .and. sp_idx <= N_SPECIES) species_counts(sp_idx) = species_counts(sp_idx) + 1
        end do
        
        ! Read emissivity - Python (n_Tc, n_Th, n_ne, n_feh, n_lines) -> Fortran reads (n_lines, n_feh, n_ne, n_Th, n_Tc)
        dims5 = [int(n_total, hsize_t), int(this%n_feh, hsize_t), int(this%n_ne, hsize_t), &
                 int(this%n_teh, hsize_t), int(this%n_tec, hsize_t)]
        allocate(emiss_all(n_total, this%n_feh, this%n_ne, this%n_teh, this%n_tec))
        call h5dopen_f(file_id, 'emiss', dset_id, hdferr)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, emiss_all, dims5, hdferr)
        call h5dclose_f(dset_id, hdferr)
        
        call h5fclose_f(file_id, hdferr)
        call h5close_f(hdferr)
        
        ! Allocate species arrays
        do j = 1, N_SPECIES
            n_sp = species_counts(j)
            if (n_sp == 0) cycle
            if (.not. allocated(this%species(j)%emissivities_double)) then
                allocate(this%species(j)%emissivities_double(n_sp, this%n_tec, this%n_teh, this%n_ne, this%n_feh))
                this%species(j)%emissivities_double = 0.0_real64
            end if
        end do
        
        ! Fill with reordering: from (line, feh, ne, Th, Tc) to (line, Tc, Th, ne, feh)
        line_counters = 0
        do i = 1, n_total
            sp_idx = species_indices(i)
            if (sp_idx < 1 .or. sp_idx > N_SPECIES) cycle
            line_counters(sp_idx) = line_counters(sp_idx) + 1
            k = line_counters(sp_idx)
            do i_tc = 1, this%n_tec
                do i_th = 1, this%n_teh
                    do i_ne = 1, this%n_ne
                        do i_feh = 1, this%n_feh
                            this%species(sp_idx)%emissivities_double(k, i_tc, i_th, i_ne, i_feh) = &
                                emiss_all(i, i_feh, i_ne, i_th, i_tc)
                        end do
                    end do
                end do
            end do
        end do
        
        deallocate(emiss_all, wavelength_all, species_all, species_indices)
        this%double_maxwellian_loaded = .true.
        
        print '(A)', 'Double Maxwellian tables loaded successfully:'
        print '(A,F8.3,A,F8.1,A)', '  Core temperature range: ', this%tec_arr(1), ' - ', this%tec_arr(this%n_tec), ' eV'
        print '(A,F8.1,A,F8.0,A)', '  Hot temperature range: ', this%teh_arr(1), ' - ', this%teh_arr(this%n_teh), ' eV'
        print '(A,F8.1,A,F8.0,A)', '  Density range: ', this%ne_arr(1), ' - ', this%ne_arr(this%n_ne), ' cm^-3'
        print '(A,ES10.3,A,ES10.3)', '  Hot fraction range: ', this%feh_arr(1), ' - ', this%feh_arr(this%n_feh)
        print '(A,I0,A,I0,A,I0,A,I0)', '  Grid size: ', this%n_tec, ' x ', this%n_teh, ' x ', this%n_ne, ' x ', this%n_feh
        write(*, '(A)', advance='no') '  Species: '
        do j = 1, N_SPECIES
            if (allocated(this%species(j)%emissivities_double)) &
                write(*, '(A,A,I0,A)', advance='no') trim(this%species(j)%name), '(', this%species(j)%n_lines, ') '
        end do
        print '(A)', ''
    end subroutine load_emission_tables_double
    
    subroutine trace_ray(this, start_pos, direction, ds, max_distance, s_values, positions, n_points)
        class(jovian_uv_raytracer), intent(in) :: this
        real(real64), intent(in) :: start_pos(3), direction(3), ds, max_distance
        real(real64), intent(out) :: s_values(:), positions(:,:)
        integer(int32), intent(out) :: n_points
        
        real(real64) :: dir_norm(3), norm_val
        integer(int32) :: i
        
        norm_val = vec_norm(direction)
        if (norm_val > 0.0_real64) then
            dir_norm = direction / norm_val
        else
            dir_norm = [0.0_real64, 1.0_real64, 0.0_real64]
        end if
        
        n_points = int(max_distance / ds) + 1
        n_points = min(n_points, size(s_values))
        
        do i = 1, n_points
            s_values(i) = real(i-1, real64) * ds
            positions(i, 1) = start_pos(1) + s_values(i) * dir_norm(1)
            positions(i, 2) = start_pos(2) + s_values(i) * dir_norm(2)
            positions(i, 3) = start_pos(3) + s_values(i) * dir_norm(3)
        end do
    end subroutine trace_ray
    
    subroutine interpolate_plasma_trilinear(this, positions, n_pos, ne_los, Te_los, &
                                             nsp_los, ns2p_los, ns3p_los, nop_los, no2p_los, feh_los, Teh_los)
        class(jovian_uv_raytracer), intent(in) :: this
        real(real64), intent(in) :: positions(:,:)
        integer(int32), intent(in) :: n_pos
        real(real64), intent(out) :: ne_los(:), Te_los(:), nsp_los(:), ns2p_los(:), ns3p_los(:)
        real(real64), intent(out) :: nop_los(:), no2p_los(:), feh_los(:), Teh_los(:)
        
        integer(int32) :: i
        real(real64) :: xp, yp, zp, val
        
        !$omp parallel do private(i, xp, yp, zp, val)
        do i = 1, n_pos
            xp = positions(i, 1); yp = positions(i, 2); zp = positions(i, 3)
            
            val = trilinear_interp(xp, yp, zp, this%x_axis, this%y_axis, this%z_axis, this%nec)
            if (.not. ieee_is_finite(val) .or. val < 0.0_real64) val = 0.0_real64
            ne_los(i) = val
            
            val = trilinear_interp(xp, yp, zp, this%x_axis, this%y_axis, this%z_axis, this%Tec)
            if (.not. ieee_is_finite(val) .or. val < 0.0_real64) val = 0.0_real64
            Te_los(i) = val
            
            val = trilinear_interp(xp, yp, zp, this%x_axis, this%y_axis, this%z_axis, this%nsp)
            if (.not. ieee_is_finite(val) .or. val < 0.0_real64) val = 0.0_real64
            nsp_los(i) = val
            
            val = trilinear_interp(xp, yp, zp, this%x_axis, this%y_axis, this%z_axis, this%ns2p)
            if (.not. ieee_is_finite(val) .or. val < 0.0_real64) val = 0.0_real64
            ns2p_los(i) = val
            
            val = trilinear_interp(xp, yp, zp, this%x_axis, this%y_axis, this%z_axis, this%ns3p)
            if (.not. ieee_is_finite(val) .or. val < 0.0_real64) val = 0.0_real64
            ns3p_los(i) = val
            
            val = trilinear_interp(xp, yp, zp, this%x_axis, this%y_axis, this%z_axis, this%nop)
            if (.not. ieee_is_finite(val) .or. val < 0.0_real64) val = 0.0_real64
            nop_los(i) = val
            
            val = trilinear_interp(xp, yp, zp, this%x_axis, this%y_axis, this%z_axis, this%no2p)
            if (.not. ieee_is_finite(val) .or. val < 0.0_real64) val = 0.0_real64
            no2p_los(i) = val
            
            val = trilinear_interp(xp, yp, zp, this%x_axis, this%y_axis, this%z_axis, this%feh)
            if (.not. ieee_is_finite(val) .or. val < 0.0_real64) val = 0.0_real64
            feh_los(i) = val
            
            val = trilinear_interp(xp, yp, zp, this%x_axis, this%y_axis, this%z_axis, this%Teh)
            if (.not. ieee_is_finite(val) .or. val < 0.0_real64) val = 0.0_real64
            Teh_los(i) = val
        end do
        !$omp end parallel do
    end subroutine interpolate_plasma_trilinear
    
    subroutine interpolate_emission_single(this, Te_arr, ne_arr, n_los, species_idx, &
                                            min_wav, max_wav, wavelengths, emission_rates, n_lines_out)
        class(jovian_uv_raytracer), intent(in) :: this
        real(real64), intent(in) :: Te_arr(:), ne_arr(:)
        integer(int32), intent(in) :: n_los, species_idx
        real(real64), intent(in) :: min_wav, max_wav
        real(real64), intent(out) :: wavelengths(:), emission_rates(:,:)
        integer(int32), intent(out) :: n_lines_out
        
        integer(int32) :: i, j, k, n_lines_all, i_T0, i_T1, i_n0, i_n1
        real(real64) :: log_T, log_n, dT, dn, w_T, w_n, w00, w01, w10, w11
        logical :: in_bounds
        integer(int32), allocatable :: line_indices(:)
        
        n_lines_out = 0
        emission_rates = 0.0_real64
        if (species_idx < 1 .or. species_idx > N_SPECIES) return
        if (this%species(species_idx)%n_lines == 0) return
        
        n_lines_all = this%species(species_idx)%n_lines
        allocate(line_indices(n_lines_all))
        
        do j = 1, n_lines_all
            if (this%species(species_idx)%wavelengths(j) >= min_wav .and. &
                this%species(species_idx)%wavelengths(j) <= max_wav) then
                n_lines_out = n_lines_out + 1
                wavelengths(n_lines_out) = this%species(species_idx)%wavelengths(j)
                line_indices(n_lines_out) = j
            end if
        end do
        
        if (n_lines_out == 0) then; deallocate(line_indices); return; end if
        
        do i = 1, n_los
            if (Te_arr(i) <= 0.0_real64 .or. ne_arr(i) <= 0.0_real64) cycle
            if (.not. ieee_is_finite(Te_arr(i)) .or. .not. ieee_is_finite(ne_arr(i))) cycle
            
            in_bounds = (Te_arr(i) >= this%temp_arr(1)) .and. (Te_arr(i) <= this%temp_arr(this%n_temp_single)) .and. &
                       (ne_arr(i) >= this%dens_arr(1)) .and. (ne_arr(i) <= this%dens_arr(this%n_dens_single))
            if (.not. in_bounds) cycle
            
            log_T = log10(Te_arr(i)); log_n = log10(ne_arr(i))
            i_T0 = max(1, min(binary_search(this%log_temp, log_T), this%n_temp_single - 1)); i_T1 = i_T0 + 1
            i_n0 = max(1, min(binary_search(this%log_dens, log_n), this%n_dens_single - 1)); i_n1 = i_n0 + 1
            
            dT = this%log_temp(i_T1) - this%log_temp(i_T0)
            dn = this%log_dens(i_n1) - this%log_dens(i_n0)
            if (dT /= 0.0_real64) then; w_T = (log_T - this%log_temp(i_T0)) / dT; else; w_T = 0.0_real64; end if
            if (dn /= 0.0_real64) then; w_n = (log_n - this%log_dens(i_n0)) / dn; else; w_n = 0.0_real64; end if
            
            w00 = (1.0_real64 - w_T) * (1.0_real64 - w_n)
            w01 = (1.0_real64 - w_T) * w_n
            w10 = w_T * (1.0_real64 - w_n)
            w11 = w_T * w_n
            
            do k = 1, n_lines_out
                j = line_indices(k)
                emission_rates(i, k) = w00 * this%species(species_idx)%emissivities_single(j, i_T0, i_n0) + &
                                       w01 * this%species(species_idx)%emissivities_single(j, i_T0, i_n1) + &
                                       w10 * this%species(species_idx)%emissivities_single(j, i_T1, i_n0) + &
                                       w11 * this%species(species_idx)%emissivities_single(j, i_T1, i_n1)
            end do
        end do
        deallocate(line_indices)
    end subroutine interpolate_emission_single
    
    subroutine interpolate_emission_double(this, Tec_arr, Teh_arr, ne_arr, feh_local, &
                                            n_los, species_idx, min_wav, max_wav, &
                                            wavelengths, emission_rates, n_lines_out)
        class(jovian_uv_raytracer), intent(in) :: this
        real(real64), intent(in) :: Tec_arr(:), Teh_arr(:), ne_arr(:), feh_local(:)
        integer(int32), intent(in) :: n_los, species_idx
        real(real64), intent(in) :: min_wav, max_wav
        real(real64), intent(out) :: wavelengths(:), emission_rates(:,:)
        integer(int32), intent(out) :: n_lines_out
        
        integer(int32) :: i, j, k, m, n_lines_all
        integer(int32) :: i0_Tc, i1_Tc, i0_Th, i1_Th, i0_ne, i1_ne, i0_feh, i1_feh
        integer(int32) :: b_Tc, b_Th, b_ne, b_feh, idx_Tc, idx_Th, idx_ne, idx_feh
        real(real64) :: log_Tc, log_Th, log_ne_val, log_feh_val
        real(real64) :: dTc, dTh, dne, dfeh, w_Tc, w_Th, w_ne, w_feh, weight
        logical :: use_single, in_bounds
        integer(int32), allocatable :: line_indices(:)
        
        n_lines_out = 0
        emission_rates = 0.0_real64
        
        if (.not. this%double_maxwellian_loaded) then
            call this%interpolate_emission_single(Tec_arr, ne_arr, n_los, species_idx, min_wav, max_wav, &
                                                   wavelengths, emission_rates, n_lines_out)
            return
        end if
        
        if (species_idx < 1 .or. species_idx > N_SPECIES) return
        if (.not. allocated(this%species(species_idx)%emissivities_double)) return
        
        n_lines_all = this%species(species_idx)%n_lines
        allocate(line_indices(n_lines_all))
        
        do j = 1, n_lines_all
            if (this%species(species_idx)%wavelengths(j) >= min_wav .and. &
                this%species(species_idx)%wavelengths(j) <= max_wav) then
                n_lines_out = n_lines_out + 1
                wavelengths(n_lines_out) = this%species(species_idx)%wavelengths(j)
                line_indices(n_lines_out) = j
            end if
        end do
        
        if (n_lines_out == 0) then; deallocate(line_indices); return; end if
        
        do i = 1, n_los
            use_single = (feh_local(i) < this%feh_arr(1)) .or. (Teh_arr(i) < this%teh_arr(1))
            if (use_single) cycle
            if (Tec_arr(i) <= 0.0_real64 .or. Teh_arr(i) <= 0.0_real64 .or. &
                ne_arr(i) <= 0.0_real64 .or. feh_local(i) < 0.0_real64) cycle
            
            in_bounds = (Tec_arr(i) >= this%tec_arr(1)) .and. (Tec_arr(i) <= this%tec_arr(this%n_tec)) .and. &
                       (Teh_arr(i) >= this%teh_arr(1)) .and. (Teh_arr(i) <= this%teh_arr(this%n_teh)) .and. &
                       (ne_arr(i) >= this%ne_arr(1)) .and. (ne_arr(i) <= this%ne_arr(this%n_ne)) .and. &
                       (feh_local(i) >= this%feh_arr(1)) .and. (feh_local(i) <= this%feh_arr(this%n_feh))
            if (.not. in_bounds) cycle
            
            log_Tc = log10(Tec_arr(i)); log_Th = log10(Teh_arr(i))
            log_ne_val = log10(ne_arr(i)); log_feh_val = log10(feh_local(i))
            
            i0_Tc = max(1, min(binary_search(this%log_tec, log_Tc), this%n_tec - 1)); i1_Tc = i0_Tc + 1
            i0_Th = max(1, min(binary_search(this%log_teh, log_Th), this%n_teh - 1)); i1_Th = i0_Th + 1
            i0_ne = max(1, min(binary_search(this%log_ne, log_ne_val), this%n_ne - 1)); i1_ne = i0_ne + 1
            i0_feh = max(1, min(binary_search(this%log_feh, log_feh_val), this%n_feh - 1)); i1_feh = i0_feh + 1
            
            dTc = this%log_tec(i1_Tc) - this%log_tec(i0_Tc)
            dTh = this%log_teh(i1_Th) - this%log_teh(i0_Th)
            dne = this%log_ne(i1_ne) - this%log_ne(i0_ne)
            dfeh = this%log_feh(i1_feh) - this%log_feh(i0_feh)
            
            if (dTc /= 0.0_real64) then; w_Tc = (log_Tc - this%log_tec(i0_Tc)) / dTc; else; w_Tc = 0.0_real64; end if
            if (dTh /= 0.0_real64) then; w_Th = (log_Th - this%log_teh(i0_Th)) / dTh; else; w_Th = 0.0_real64; end if
            if (dne /= 0.0_real64) then; w_ne = (log_ne_val - this%log_ne(i0_ne)) / dne; else; w_ne = 0.0_real64; end if
            if (dfeh /= 0.0_real64) then; w_feh = (log_feh_val - this%log_feh(i0_feh)) / dfeh; else; w_feh = 0.0_real64; end if
            
            do m = 1, n_lines_out
                k = line_indices(m)
                do b_Tc = 0, 1
                    do b_Th = 0, 1
                        do b_ne = 0, 1
                            do b_feh = 0, 1
                                if (b_Tc == 0) then; idx_Tc = i0_Tc; weight = 1.0_real64 - w_Tc
                                else; idx_Tc = i1_Tc; weight = w_Tc; end if
                                if (b_Th == 0) then; idx_Th = i0_Th; weight = weight * (1.0_real64 - w_Th)
                                else; idx_Th = i1_Th; weight = weight * w_Th; end if
                                if (b_ne == 0) then; idx_ne = i0_ne; weight = weight * (1.0_real64 - w_ne)
                                else; idx_ne = i1_ne; weight = weight * w_ne; end if
                                if (b_feh == 0) then; idx_feh = i0_feh; weight = weight * (1.0_real64 - w_feh)
                                else; idx_feh = i1_feh; weight = weight * w_feh; end if
                                
                                emission_rates(i, m) = emission_rates(i, m) + &
                                    weight * this%species(species_idx)%emissivities_double(k, idx_Tc, idx_Th, idx_ne, idx_feh)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        deallocate(line_indices)
    end subroutine interpolate_emission_double
    
    subroutine integrate_species_single(this, s_values, n_pos, species_idx, Te_los, ne_los, n_ion_los, &
                                         min_wav, max_wav, wavelengths, brightnesses, n_lines)
        class(jovian_uv_raytracer), intent(in) :: this
        real(real64), intent(in) :: s_values(:)
        integer(int32), intent(in) :: n_pos, species_idx
        real(real64), intent(in) :: Te_los(:), ne_los(:), n_ion_los(:)
        real(real64), intent(in) :: min_wav, max_wav
        real(real64), intent(out) :: wavelengths(:), brightnesses(:)
        integer(int32), intent(out) :: n_lines
        
        real(real64), allocatable :: emission_rates(:,:), integrand(:)
        integer(int32) :: j
        
        allocate(emission_rates(n_pos, MAX_LINES))
        emission_rates = 0.0_real64
        
        call this%interpolate_emission_single(Te_los, ne_los, n_pos, species_idx, min_wav, max_wav, &
                                               wavelengths, emission_rates, n_lines)
        
        if (n_lines == 0) then; deallocate(emission_rates); return; end if
        
        brightnesses = 0.0_real64
        allocate(integrand(n_pos))
        
        do j = 1, n_lines
            integrand = emission_rates(1:n_pos, j) * n_ion_los(1:n_pos)
            brightnesses(j) = RAYLEIGH_FACTOR * simpson_integrate(integrand, s_values(1:n_pos))
        end do
        
        deallocate(emission_rates, integrand)
    end subroutine integrate_species_single
    
    subroutine integrate_species_double(this, s_values, n_pos, species_idx, Tec_los, Teh_los, ne_los, feh_los, n_ion_los, &
                                         min_wav, max_wav, wavelengths, brightnesses, n_lines)
        class(jovian_uv_raytracer), intent(in) :: this
        real(real64), intent(in) :: s_values(:)
        integer(int32), intent(in) :: n_pos, species_idx
        real(real64), intent(in) :: Tec_los(:), Teh_los(:), ne_los(:), feh_los(:), n_ion_los(:)
        real(real64), intent(in) :: min_wav, max_wav
        real(real64), intent(out) :: wavelengths(:), brightnesses(:)
        integer(int32), intent(out) :: n_lines
        
        real(real64), allocatable :: emission_rates(:,:), integrand(:)
        integer(int32) :: j
        
        allocate(emission_rates(n_pos, MAX_LINES))
        emission_rates = 0.0_real64
        
        call this%interpolate_emission_double(Tec_los, Teh_los, ne_los, feh_los, n_pos, species_idx, min_wav, max_wav, &
                                               wavelengths, emission_rates, n_lines)
        
        if (n_lines == 0) then; deallocate(emission_rates); return; end if
        
        brightnesses = 0.0_real64
        allocate(integrand(n_pos))
        
        do j = 1, n_lines
            integrand = emission_rates(1:n_pos, j) * n_ion_los(1:n_pos)
            brightnesses(j) = RAYLEIGH_FACTOR * simpson_integrate(integrand, s_values(1:n_pos))
        end do
        
        deallocate(emission_rates, integrand)
    end subroutine integrate_species_double
    
    subroutine convolve_spectrum_erf(this, wavelength_grid, n_bins, bin_width, line_wavelengths, &
                                      line_brightnesses, n_lines, fwhm, spectrum)
        class(jovian_uv_raytracer), intent(in) :: this
        real(real64), intent(in) :: wavelength_grid(:), line_wavelengths(:), line_brightnesses(:)
        integer(int32), intent(in) :: n_bins, n_lines
        real(real64), intent(in) :: bin_width, fwhm
        real(real64), intent(out) :: spectrum(:)
        
        real(real64) :: sigma, sigma_sqrt2, bin_lo, bin_hi, erf_lo, erf_hi, contrib
        integer(int32) :: i, j
        
        spectrum = 0.0_real64
        if (n_lines == 0) return
        
        sigma = fwhm / (2.0_real64 * sqrt(2.0_real64 * LN2))
        sigma_sqrt2 = sigma * SQRT2
        
        !$omp parallel do private(i, j, bin_lo, bin_hi, erf_lo, erf_hi, contrib)
        do i = 1, n_bins
            bin_lo = wavelength_grid(i) - bin_width / 2.0_real64
            bin_hi = wavelength_grid(i) + bin_width / 2.0_real64
            
            do j = 1, n_lines
                if (line_brightnesses(j) <= 0.0_real64) cycle
                erf_lo = (bin_lo - line_wavelengths(j)) / sigma_sqrt2
                erf_hi = (bin_hi - line_wavelengths(j)) / sigma_sqrt2
                contrib = 0.5_real64 * (erf(erf_hi) - erf(erf_lo)) * line_brightnesses(j)
                spectrum(i) = spectrum(i) + contrib
            end do
            spectrum(i) = spectrum(i) / bin_width
        end do
        !$omp end parallel do
    end subroutine convolve_spectrum_erf
    
    subroutine calculate_spectrum_single(this, slit_pos_vec, norm_vec, wavelength_range, bin_width, fwhm, ds, &
                                          wave_bins, spectrum, line_list, n_bins, n_lines_total)
        class(jovian_uv_raytracer), intent(in) :: this
        real(real64), intent(in) :: slit_pos_vec(3), norm_vec(3), wavelength_range(2), bin_width, fwhm, ds
        real(real64), intent(out) :: wave_bins(:), spectrum(:)
        type(emission_line), intent(out) :: line_list(:)
        integer(int32), intent(out) :: n_bins, n_lines_total
        
        real(real64) :: min_wav, max_wav, max_distance
        integer(int32) :: n_wav, n_pos, i, j, species_idx, n_lines_sp, n_all
        real(real64), allocatable :: s_values(:), positions(:,:)
        real(real64), allocatable :: ne_los(:), Te_los(:), nsp_los(:), ns2p_los(:), ns3p_los(:)
        real(real64), allocatable :: nop_los(:), no2p_los(:), feh_los(:), Teh_los(:), n_ion_los(:)
        real(real64), allocatable :: wavelengths_sp(:), brightnesses_sp(:)
        real(real64), allocatable :: all_wavelengths(:), all_brightnesses(:)
        integer(int32), allocatable :: all_species(:)
        
        print '(A)', 'Calculating single Maxwellian spectrum:'
        print '(A,F6.1,A,F6.1,A,F6.1,A)', '  LOS start: [', slit_pos_vec(1), ', ', slit_pos_vec(2), ', ', slit_pos_vec(3), '] R_J'
        print '(A,F5.2,A,F5.2,A,F5.2,A)', '  Direction: [', norm_vec(1), ', ', norm_vec(2), ', ', norm_vec(3), ']'
        
        min_wav = wavelength_range(1); max_wav = wavelength_range(2); max_distance = 40.0_real64
        
        n_wav = int((max_wav - min_wav) / bin_width) + 1
        n_bins = n_wav
        do i = 1, n_wav
            wave_bins(i) = min_wav + real(i-1, real64) * bin_width
        end do
        
        allocate(s_values(MAX_LOS_POINTS), positions(MAX_LOS_POINTS, 3))
        allocate(ne_los(MAX_LOS_POINTS), Te_los(MAX_LOS_POINTS))
        allocate(nsp_los(MAX_LOS_POINTS), ns2p_los(MAX_LOS_POINTS), ns3p_los(MAX_LOS_POINTS))
        allocate(nop_los(MAX_LOS_POINTS), no2p_los(MAX_LOS_POINTS))
        allocate(feh_los(MAX_LOS_POINTS), Teh_los(MAX_LOS_POINTS), n_ion_los(MAX_LOS_POINTS))
        allocate(wavelengths_sp(MAX_LINES), brightnesses_sp(MAX_LINES))
        allocate(all_wavelengths(MAX_LINES), all_brightnesses(MAX_LINES), all_species(MAX_LINES))
        
        call this%trace_ray(slit_pos_vec, norm_vec, ds, max_distance, s_values, positions, n_pos)
        call this%interpolate_plasma_trilinear(positions, n_pos, ne_los, Te_los, nsp_los, ns2p_los, ns3p_los, &
                                                nop_los, no2p_los, feh_los, Teh_los)
        
        n_all = 0
        do species_idx = 1, N_SPECIES
            if (this%species(species_idx)%n_lines == 0) cycle
            
            select case(species_idx)
            case(IDX_SP);  n_ion_los(1:n_pos) = nsp_los(1:n_pos)
            case(IDX_S2P); n_ion_los(1:n_pos) = ns2p_los(1:n_pos)
            case(IDX_S3P); n_ion_los(1:n_pos) = ns3p_los(1:n_pos)
            case(IDX_OP);  n_ion_los(1:n_pos) = nop_los(1:n_pos)
            case(IDX_O2P); n_ion_los(1:n_pos) = no2p_los(1:n_pos)
            end select
            
            if (.not. any(ne_los(1:n_pos) > 0.0_real64 .and. Te_los(1:n_pos) > 0.0_real64 .and. &
                         n_ion_los(1:n_pos) > 0.0_real64)) cycle
            
            call this%integrate_species_single(s_values, n_pos, species_idx, Te_los, ne_los, n_ion_los, &
                                                min_wav, max_wav, wavelengths_sp, brightnesses_sp, n_lines_sp)
            
            do j = 1, n_lines_sp
                n_all = n_all + 1
                all_wavelengths(n_all) = wavelengths_sp(j)
                all_brightnesses(n_all) = brightnesses_sp(j)
                all_species(n_all) = species_idx
            end do
        end do
        
        n_lines_total = 0
        do i = 1, n_all
            if (all_brightnesses(i) > 1.0e-10_real64) then
                n_lines_total = n_lines_total + 1
                line_list(n_lines_total)%wavelength = all_wavelengths(i)
                line_list(n_lines_total)%brightness = all_brightnesses(i)
                line_list(n_lines_total)%species_idx = all_species(i)
            end if
        end do
        
        print '(A,I0,A,I0,A)', '  Processed ', n_all, ' lines, ', n_lines_total, ' with non-zero brightness'
        
        call this%convolve_spectrum_erf(wave_bins, n_bins, bin_width, all_wavelengths, all_brightnesses, n_all, fwhm, spectrum)
        
        deallocate(s_values, positions, ne_los, Te_los, nsp_los, ns2p_los, ns3p_los, nop_los, no2p_los)
        deallocate(feh_los, Teh_los, n_ion_los, wavelengths_sp, brightnesses_sp)
        deallocate(all_wavelengths, all_brightnesses, all_species)
    end subroutine calculate_spectrum_single
    
    subroutine calculate_spectrum_double(this, slit_pos_vec, norm_vec, wavelength_range, bin_width, fwhm, ds, &
                                          wave_bins, spectrum, line_list, n_bins, n_lines_total)
        class(jovian_uv_raytracer), intent(in) :: this
        real(real64), intent(in) :: slit_pos_vec(3), norm_vec(3), wavelength_range(2), bin_width, fwhm, ds
        real(real64), intent(out) :: wave_bins(:), spectrum(:)
        type(emission_line), intent(out) :: line_list(:)
        integer(int32), intent(out) :: n_bins, n_lines_total
        
        real(real64) :: min_wav, max_wav, max_distance
        integer(int32) :: n_wav, n_pos, i, j, species_idx, n_lines_sp, n_all
        real(real64), allocatable :: s_values(:), positions(:,:)
        real(real64), allocatable :: ne_los(:), Te_los(:), nsp_los(:), ns2p_los(:), ns3p_los(:)
        real(real64), allocatable :: nop_los(:), no2p_los(:), feh_los(:), Teh_los(:), n_ion_los(:)
        real(real64), allocatable :: wavelengths_sp(:), brightnesses_sp(:)
        real(real64), allocatable :: all_wavelengths(:), all_brightnesses(:)
        integer(int32), allocatable :: all_species(:)
        
        if (.not. this%double_maxwellian_loaded) then
            print '(A)', 'ERROR: Double Maxwellian tables not loaded'
            n_bins = 0; n_lines_total = 0; return
        end if
        
        print '(A)', 'Calculating double Maxwellian spectrum:'
        print '(A,F6.1,A,F6.1,A,F6.1,A)', '  LOS start: [', slit_pos_vec(1), ', ', slit_pos_vec(2), ', ', slit_pos_vec(3), '] R_J'
        print '(A,F5.2,A,F5.2,A,F5.2,A)', '  Direction: [', norm_vec(1), ', ', norm_vec(2), ', ', norm_vec(3), ']'
        
        min_wav = wavelength_range(1); max_wav = wavelength_range(2); max_distance = 40.0_real64
        
        n_wav = int((max_wav - min_wav) / bin_width) + 1
        n_bins = n_wav
        do i = 1, n_wav
            wave_bins(i) = min_wav + real(i-1, real64) * bin_width
        end do
        
        allocate(s_values(MAX_LOS_POINTS), positions(MAX_LOS_POINTS, 3))
        allocate(ne_los(MAX_LOS_POINTS), Te_los(MAX_LOS_POINTS))
        allocate(nsp_los(MAX_LOS_POINTS), ns2p_los(MAX_LOS_POINTS), ns3p_los(MAX_LOS_POINTS))
        allocate(nop_los(MAX_LOS_POINTS), no2p_los(MAX_LOS_POINTS))
        allocate(feh_los(MAX_LOS_POINTS), Teh_los(MAX_LOS_POINTS), n_ion_los(MAX_LOS_POINTS))
        allocate(wavelengths_sp(MAX_LINES), brightnesses_sp(MAX_LINES))
        allocate(all_wavelengths(MAX_LINES), all_brightnesses(MAX_LINES), all_species(MAX_LINES))
        
        call this%trace_ray(slit_pos_vec, norm_vec, ds, max_distance, s_values, positions, n_pos)
        call this%interpolate_plasma_trilinear(positions, n_pos, ne_los, Te_los, nsp_los, ns2p_los, ns3p_los, &
                                                nop_los, no2p_los, feh_los, Teh_los)
        
        ! For double Maxwellian, use total ne
        do i = 1, n_pos
            ne_los(i) = trilinear_interp(positions(i,1), positions(i,2), positions(i,3), &
                                          this%x_axis, this%y_axis, this%z_axis, this%ne_total)
        end do
        
        n_all = 0
        do species_idx = 1, N_SPECIES
            if (this%species(species_idx)%n_lines == 0) cycle
            if (.not. allocated(this%species(species_idx)%emissivities_double)) cycle
            
            select case(species_idx)
            case(IDX_SP);  n_ion_los(1:n_pos) = nsp_los(1:n_pos)
            case(IDX_S2P); n_ion_los(1:n_pos) = ns2p_los(1:n_pos)
            case(IDX_S3P); n_ion_los(1:n_pos) = ns3p_los(1:n_pos)
            case(IDX_OP);  n_ion_los(1:n_pos) = nop_los(1:n_pos)
            case(IDX_O2P); n_ion_los(1:n_pos) = no2p_los(1:n_pos)
            end select
            
            if (.not. any(ne_los(1:n_pos) > 0.0_real64 .and. Te_los(1:n_pos) > 0.0_real64 .and. &
                         n_ion_los(1:n_pos) > 0.0_real64)) cycle
            
            call this%integrate_species_double(s_values, n_pos, species_idx, Te_los, Teh_los, ne_los, feh_los, n_ion_los, &
                                                min_wav, max_wav, wavelengths_sp, brightnesses_sp, n_lines_sp)
            
            do j = 1, n_lines_sp
                n_all = n_all + 1
                all_wavelengths(n_all) = wavelengths_sp(j)
                all_brightnesses(n_all) = brightnesses_sp(j)
                all_species(n_all) = species_idx
            end do
        end do
        
        n_lines_total = 0
        do i = 1, n_all
            if (all_brightnesses(i) > 1.0e-10_real64) then
                n_lines_total = n_lines_total + 1
                line_list(n_lines_total)%wavelength = all_wavelengths(i)
                line_list(n_lines_total)%brightness = all_brightnesses(i)
                line_list(n_lines_total)%species_idx = all_species(i)
            end if
        end do
        
        print '(A,I0,A,I0,A)', '  Processed ', n_all, ' lines, ', n_lines_total, ' with non-zero brightness'
        
        call this%convolve_spectrum_erf(wave_bins, n_bins, bin_width, all_wavelengths, all_brightnesses, n_all, fwhm, spectrum)
        
        deallocate(s_values, positions, ne_los, Te_los, nsp_los, ns2p_los, ns3p_los, nop_los, no2p_los)
        deallocate(feh_los, Teh_los, n_ion_los, wavelengths_sp, brightnesses_sp)
        deallocate(all_wavelengths, all_brightnesses, all_species)
    end subroutine calculate_spectrum_double
    
    subroutine cleanup(this)
        class(jovian_uv_raytracer), intent(inout) :: this
        integer(int32) :: j
        
        if (allocated(this%x_axis)) deallocate(this%x_axis)
        if (allocated(this%y_axis)) deallocate(this%y_axis)
        if (allocated(this%z_axis)) deallocate(this%z_axis)
        if (allocated(this%nec)) deallocate(this%nec)
        if (allocated(this%neh)) deallocate(this%neh)
        if (allocated(this%ne_total)) deallocate(this%ne_total)
        if (allocated(this%feh)) deallocate(this%feh)
        if (allocated(this%nsp)) deallocate(this%nsp)
        if (allocated(this%ns2p)) deallocate(this%ns2p)
        if (allocated(this%ns3p)) deallocate(this%ns3p)
        if (allocated(this%nop)) deallocate(this%nop)
        if (allocated(this%no2p)) deallocate(this%no2p)
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
            if (allocated(this%species(j)%emissivities_single)) deallocate(this%species(j)%emissivities_single)
            if (allocated(this%species(j)%emissivities_double)) deallocate(this%species(j)%emissivities_double)
        end do
        
        this%plasma_loaded = .false.
        this%single_maxwellian_loaded = .false.
        this%double_maxwellian_loaded = .false.
    end subroutine cleanup
    
end module ipt_emission_raytracer

module ipt_emission_plotting
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use ipt_emission_constants
    use ipt_emission_types
    use ipt_emission_utils
    implicit none
contains
    
    subroutine plot_spectrum_gnuplot(wave_bins, spectrum, n_bins, line_list, n_lines, &
                                      title, color_index, output_file)
        real(real64), intent(in) :: wave_bins(:), spectrum(:)
        integer(int32), intent(in) :: n_bins
        type(emission_line), intent(in) :: line_list(:)
        integer(int32), intent(in) :: n_lines
        character(len=*), intent(in) :: title, output_file
        integer(int32), intent(in) :: color_index
        
        integer :: unit_data, unit_gp, i, ios
        character(len=256) :: data_file, gp_file
        character(len=32) :: color_str
        real(real64) :: total_brightness, y_max
        
        if (color_index == 0) then; color_str = '#1f77b4'; else; color_str = '#ff7f0e'; end if
        
        data_file = trim(output_file) // '.dat'
        gp_file = trim(output_file) // '.gp'
        
        open(newunit=unit_data, file=trim(data_file), status='replace', iostat=ios)
        if (ios /= 0) then; print *, 'Error: Could not open data file'; return; end if
        write(unit_data, '(A)') '# Wavelength [A]    Brightness [R/A]'
        do i = 1, n_bins
            write(unit_data, '(F12.4, ES16.6)') wave_bins(i), spectrum(i)
        end do
        close(unit_data)
        
        total_brightness = simpson_integrate(spectrum(1:n_bins), wave_bins(1:n_bins))
        y_max = maxval(spectrum(1:n_bins)) * 1.1_real64
        if (y_max < 1.0e-10_real64) y_max = 1.0_real64
        
        open(newunit=unit_gp, file=trim(gp_file), status='replace', iostat=ios)
        if (ios /= 0) then; print *, 'Error: Could not open gnuplot script'; return; end if
        
        write(unit_gp, '(A)') 'set terminal pngcairo size 1200,600 enhanced font "Arial,12"'
        write(unit_gp, '(A)') 'set output "' // trim(output_file) // '.png"'
        write(unit_gp, '(A)') 'set title "' // trim(title) // '" font "Arial,14"'
        write(unit_gp, '(A)') 'set xlabel "Wavelength [Angstrom]" font "Arial,12"'
        write(unit_gp, '(A)') 'set ylabel "Brightness [Rayleighs/Angstrom]" font "Arial,12"'
        write(unit_gp, '(A,F8.1,A,F8.1,A)') 'set xrange [', wave_bins(1), ':', wave_bins(n_bins), ']'
        write(unit_gp, '(A,ES10.3,A)') 'set yrange [0:', y_max, ']'
        write(unit_gp, '(A)') 'set grid'
        write(unit_gp, '(A)') 'set key off'
        write(unit_gp, '(A,F8.1,A)') 'set label 1 "Total: ', total_brightness, ' R" at graph 0.02,0.95 left front'
        write(unit_gp, '(A)') 'plot "' // trim(data_file) // '" using 1:2 with lines lc rgb "' // &
            trim(color_str) // '" lw 1.5 notitle'
        close(unit_gp)
        
        call execute_command_line('gnuplot ' // trim(gp_file), wait=.true.)
        print '(A,A)', 'Saved: ', trim(output_file) // '.png'
    end subroutine plot_spectrum_gnuplot
    
    subroutine plot_spectrum_comparison_gnuplot(wave_bins, spectrum1, spectrum2, n_bins, &
                                                 label1, label2, total1, total2, title, output_file)
        real(real64), intent(in) :: wave_bins(:), spectrum1(:), spectrum2(:)
        integer(int32), intent(in) :: n_bins
        character(len=*), intent(in) :: label1, label2, title, output_file
        real(real64), intent(in) :: total1, total2
        
        integer :: unit_data, unit_gp, i, ios
        character(len=256) :: data_file, gp_file
        real(real64) :: y_max
        
        data_file = trim(output_file) // '.dat'
        gp_file = trim(output_file) // '.gp'
        
        open(newunit=unit_data, file=trim(data_file), status='replace', iostat=ios)
        if (ios /= 0) return
        write(unit_data, '(A)') '# Wavelength    Spectrum1    Spectrum2'
        do i = 1, n_bins
            write(unit_data, '(F12.4, 2ES16.6)') wave_bins(i), spectrum1(i), spectrum2(i)
        end do
        close(unit_data)
        
        y_max = max(maxval(spectrum1(1:n_bins)), maxval(spectrum2(1:n_bins))) * 1.1_real64
        if (y_max < 1.0e-10_real64) y_max = 1.0_real64
        
        open(newunit=unit_gp, file=trim(gp_file), status='replace', iostat=ios)
        if (ios /= 0) return
        
        write(unit_gp, '(A)') 'set terminal pngcairo size 1200,600 enhanced font "Arial,12"'
        write(unit_gp, '(A)') 'set output "' // trim(output_file) // '.png"'
        write(unit_gp, '(A)') 'set title "' // trim(title) // '" font "Arial,14"'
        write(unit_gp, '(A)') 'set xlabel "Wavelength [Angstrom]" font "Arial,12"'
        write(unit_gp, '(A)') 'set ylabel "Brightness [Rayleighs/Angstrom]" font "Arial,12"'
        write(unit_gp, '(A,F8.1,A,F8.1,A)') 'set xrange [', wave_bins(1), ':', wave_bins(n_bins), ']'
        write(unit_gp, '(A,ES10.3,A)') 'set yrange [0:', y_max, ']'
        write(unit_gp, '(A)') 'set grid'
        write(unit_gp, '(A)') 'set key top right'
        write(unit_gp, '(A,F8.1,A)') 'set label 1 "', total1, ' R" at graph 0.02,0.95 left front tc rgb "#1f77b4"'
        write(unit_gp, '(A,F8.1,A)') 'set label 2 "', total2, ' R" at graph 0.02,0.90 left front tc rgb "#ff7f0e"'
        write(unit_gp, '(A)') 'plot "' // trim(data_file) // '" using 1:2 with lines lc rgb "#1f77b4" lw 1.5 title "' // &
            trim(label1) // '", \'
        write(unit_gp, '(A)') '     "" using 1:3 with lines lc rgb "#ff7f0e" lw 1.5 title "' // trim(label2) // '"'
        close(unit_gp)
        
        call execute_command_line('gnuplot ' // trim(gp_file), wait=.true.)
        print '(A,A)', 'Saved: ', trim(output_file) // '.png'
    end subroutine plot_spectrum_comparison_gnuplot
    
    subroutine plot_single_vs_double_gnuplot(wave_bins, spectrum_single, spectrum_double, n_bins, &
                                              total_single, total_double, enhancement, title, output_file)
        real(real64), intent(in) :: wave_bins(:), spectrum_single(:), spectrum_double(:)
        integer(int32), intent(in) :: n_bins
        real(real64), intent(in) :: total_single, total_double, enhancement
        character(len=*), intent(in) :: title, output_file
        
        integer :: unit_data, unit_gp, i, ios
        character(len=256) :: data_file, gp_file
        real(real64) :: y_max
        
        data_file = trim(output_file) // '.dat'
        gp_file = trim(output_file) // '.gp'
        
        open(newunit=unit_data, file=trim(data_file), status='replace', iostat=ios)
        if (ios /= 0) return
        write(unit_data, '(A)') '# Wavelength    Single    Double'
        do i = 1, n_bins
            write(unit_data, '(F12.4, 2ES16.6)') wave_bins(i), spectrum_single(i), spectrum_double(i)
        end do
        close(unit_data)
        
        y_max = max(maxval(spectrum_single(1:n_bins)), maxval(spectrum_double(1:n_bins))) * 1.1_real64
        if (y_max < 1.0e-10_real64) y_max = 1.0_real64
        
        open(newunit=unit_gp, file=trim(gp_file), status='replace', iostat=ios)
        if (ios /= 0) return
        
        write(unit_gp, '(A)') 'set terminal pngcairo size 1200,600 enhanced font "Arial,12"'
        write(unit_gp, '(A)') 'set output "' // trim(output_file) // '.png"'
        write(unit_gp, '(A)') 'set title "' // trim(title) // '" font "Arial,14"'
        write(unit_gp, '(A)') 'set xlabel "Wavelength [Angstrom]" font "Arial,12"'
        write(unit_gp, '(A)') 'set ylabel "Brightness [Rayleighs/Angstrom]" font "Arial,12"'
        write(unit_gp, '(A,F8.1,A,F8.1,A)') 'set xrange [', wave_bins(1), ':', wave_bins(n_bins), ']'
        write(unit_gp, '(A,ES10.3,A)') 'set yrange [0:', y_max, ']'
        write(unit_gp, '(A)') 'set grid'
        write(unit_gp, '(A)') 'set key top right'
        write(unit_gp, '(A,F8.1,A)') 'set label 1 "Single: ', total_single, ' R" at graph 0.02,0.95 left front'
        write(unit_gp, '(A,F8.1,A)') 'set label 2 "Double: ', total_double, ' R" at graph 0.02,0.90 left front'
        write(unit_gp, '(A,F6.3,A)') 'set label 3 "Enhancement: ', enhancement, 'x" at graph 0.02,0.85 left front'
        write(unit_gp, '(A)') 'plot "' // trim(data_file) // '" using 1:2 with lines lc rgb "#1f77b4" lw 1.5 title "Single Maxwellian", \'
        write(unit_gp, '(A)') '     "" using 1:3 with lines lc rgb "#ff7f0e" lw 1.5 title "Double Maxwellian"'
        close(unit_gp)
        
        call execute_command_line('gnuplot ' // trim(gp_file), wait=.true.)
        print '(A,A)', 'Saved: ', trim(output_file) // '.png'
    end subroutine plot_single_vs_double_gnuplot
    
end module ipt_emission_plotting
