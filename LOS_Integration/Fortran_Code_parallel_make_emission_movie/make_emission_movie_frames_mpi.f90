!==============================================================================
! make_emission_movie_frames_mpi.f90
!
! MPI-Parallelized Hybrid Optical+UV Emission Movie Frame Generator
! ==================================================================
!
! This program generates emission data for movie frames showing three specific
! emission lines from the Jovian plasma torus in a stacked 3-panel view:
!
!     Top panel:    S+ 6731 Å (optical)
!     Middle panel: S++ 680 Å (UV)
!     Bottom panel: O+ 833 Å (UV)
!
! PARALLELIZATION STRATEGY:
!   - MPI parallelization over lines of sight within each frame
!   - Each MPI rank loads a shared copy of the plasma model and emission tables
!   - Work is distributed across ranks using block decomposition
!   - Rank 0 collects results and writes binary output files
!   - Frame generation is sequential; LOS calculations are parallel
!
! OUTPUT FORMAT:
!   Binary files in output_data/ directory:
!     frame_XXXX_s1p_6731.bin  - S+ 6731 Å emission data
!     frame_XXXX_s2p_680.bin   - S++ 680 Å emission data
!     frame_XXXX_op_833.bin    - O+ 833 Å emission data
!     grid_info.bin            - Grid coordinates and metadata
!
!   A separate Python script reads these files to generate PNG frames.
!
! MOVIE PARAMETERS:
!   - Frame rate: 15 fps
!   - Duration: 5 seconds
!   - Total frames: 75 (covering 360° rotation)
!   - CML step: 4.8° per frame
!
! SPECTRAL PARAMETERS:
!   S+ 6731 Å (Optical):
!     - FWHM: 1.2 Å (narrow optical spectrometer)
!     - Bin width: 0.6 Å
!     - Window: 6721 - 6741 Å
!
!   S++ 680 Å (UV):
!     - FWHM: 6.0 Å (broader UV)
!     - Bin width: 1.0 Å
!     - Window: 670 - 690 Å
!
!   O+ 833 Å (UV):
!     - FWHM: 6.0 Å
!     - Bin width: 1.0 Å
!     - Window: 823 - 843 Å
!
! VIEWING GEOMETRY:
!   For each frame at CML angle θ (degrees):
!     CML = 0°:   Observer at X = -R_obs, looking in +X direction
!     CML = 90°:  Observer at Y = +R_obs, looking in -Y direction
!     CML = 180°: Observer at X = +R_obs, looking in -X direction
!     CML = 270°: Observer at Y = -R_obs, looking in +Y direction
!
! COLORBAR LIMITS:
!   - S+ 6731 Å: Determined by prescan of 4 CML angles with 6% headroom
!   - S++ 680 Å: Fixed (0 - 185 R)
!   - O+ 833 Å: Fixed (0 - 175 R)
!
! COMPILATION:
!   mpifort -O3 -march=native -o make_emission_movie_frames_mpi \
!       IPT_emiss_MOP_community_code.f90 make_emission_movie_frames_mpi.f90 \
!       -I$(HDF5_DIR)/include -L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5
!
! USAGE:
!   mpirun -np <num_procs> ./make_emission_movie_frames_mpi
!
! AUTHOR: Edward (Eddie) G. Nerney
! INSTITUTION: Laboratory for Atmospheric and Space Physics, CU Boulder
! LICENSE: Open source for academic and research use (MOP Community Code)
! DATE: November 2025
!==============================================================================

program make_emission_movie_frames_mpi
    use mpi
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use ipt_emission_constants
    use ipt_emission_utils
    use ipt_emission_types
    use ipt_emission_raytracer
    implicit none
    
    !--------------------------------------------------------------------------
    ! MPI Variables
    !--------------------------------------------------------------------------
    integer :: ierr, rank, nprocs
    integer :: status(MPI_STATUS_SIZE)
    
    !--------------------------------------------------------------------------
    ! Movie Parameters
    !--------------------------------------------------------------------------
    integer(int32), parameter :: FRAME_RATE = 30 !15
    integer(int32), parameter :: DURATION = 8 !5
    integer(int32), parameter :: N_FRAMES = FRAME_RATE * DURATION  !240 Frames ! 75 frames
    
    !--------------------------------------------------------------------------
    ! Spectral Parameters
    !--------------------------------------------------------------------------
    ! S+ 6731 Å (Optical)
    real(real64), parameter :: S1P_6731_WAV = 6731.0_real64
    real(real64), parameter :: S1P_6731_FWHM = 1.2_real64
    real(real64), parameter :: S1P_6731_BIN_WIDTH = 0.6_real64
    real(real64), parameter :: S1P_6731_RANGE_MIN = 6721.0_real64
    real(real64), parameter :: S1P_6731_RANGE_MAX = 6741.0_real64
    
    ! S++ 680 Å (UV)
    real(real64), parameter :: S2P_680_WAV = 680.0_real64
    real(real64), parameter :: S2P_680_FWHM = 6.0_real64
    real(real64), parameter :: S2P_680_BIN_WIDTH = 1.0_real64
    real(real64), parameter :: S2P_680_RANGE_MIN = 670.0_real64
    real(real64), parameter :: S2P_680_RANGE_MAX = 690.0_real64
    
    ! O+ 833 Å (UV)
    real(real64), parameter :: OP_833_WAV = 833.0_real64
    real(real64), parameter :: OP_833_FWHM = 6.0_real64
    real(real64), parameter :: OP_833_BIN_WIDTH = 1.0_real64
    real(real64), parameter :: OP_833_RANGE_MIN = 823.0_real64
    real(real64), parameter :: OP_833_RANGE_MAX = 843.0_real64
    
    !--------------------------------------------------------------------------
    ! Ray Tracing Parameters
    !--------------------------------------------------------------------------
    real(real64), parameter :: R_OBS = 20.0_real64  ! Observer distance [R_J]
    real(real64), parameter :: DS = 0.1_real64      ! Integration step [R_J]
    
    !--------------------------------------------------------------------------
    ! Manual UV Colorbar Limits
    !--------------------------------------------------------------------------
    real(real64), parameter :: S2P_680_VMAX = 185.0_real64  ! Rayleighs
    real(real64), parameter :: OP_833_VMAX = 175.0_real64   ! Rayleighs
    
    !--------------------------------------------------------------------------
    ! Prescan CML Angles for Optical Colorbar
    !--------------------------------------------------------------------------
    integer(int32), parameter :: N_PRESCAN = 4
    real(real64), parameter :: PRESCAN_ANGLES(N_PRESCAN) = &
        [0.0_real64, 90.0_real64, 180.0_real64, 270.0_real64]
    
    !--------------------------------------------------------------------------
    ! Grid Arrays
    !--------------------------------------------------------------------------
    real(real64), allocatable :: rho_grid(:), z_grid(:)
    integer(int32) :: nrho, nz, total_points
    
    !--------------------------------------------------------------------------
    ! Emission Data Arrays
    !--------------------------------------------------------------------------
    real(real64), allocatable :: s1p_6731_emission(:,:)
    real(real64), allocatable :: s2p_680_emission(:,:)
    real(real64), allocatable :: op_833_emission(:,:)
    real(real64), allocatable :: local_s1p(:), local_s2p(:), local_op(:)
    real(real64), allocatable :: global_s1p(:), global_s2p(:), global_op(:)
    
    !--------------------------------------------------------------------------
    ! Raytracer Instance
    !--------------------------------------------------------------------------
    type(jovian_uv_raytracer) :: raytracer
    
    !--------------------------------------------------------------------------
    ! Work Distribution Variables
    !--------------------------------------------------------------------------
    integer(int32), allocatable :: work_indices(:,:)  ! (2, total_points) - (i,j) pairs
    integer(int32) :: local_start, local_end, local_count
    integer(int32), allocatable :: counts(:), displs(:)
    
    !--------------------------------------------------------------------------
    ! Loop Variables
    !--------------------------------------------------------------------------
    integer(int32) :: frame_idx, i, j, k, idx, load_err
    real(real64) :: cml_deg, cml_rad
    real(real64) :: rho_sky, z_sky
    real(real64) :: slit_pos(3), norm_vec(3)
    real(real64) :: s1p_val, s2p_val, op_val
    real(real64) :: optical_vmax, prescan_max
    real(real64), allocatable :: cml_values(:)
    
    !--------------------------------------------------------------------------
    ! Timing Variables
    !--------------------------------------------------------------------------
    real(real64) :: start_time, frame_start, frame_time
    real(real64) :: elapsed, remaining, total_time
    
    !--------------------------------------------------------------------------
    ! File I/O
    !--------------------------------------------------------------------------
    character(len=512) :: output_dir, filename
    character(len=64) :: frame_str
    integer :: file_unit
    logical :: dir_exists
    
    !--------------------------------------------------------------------------
    ! Initialize MPI
    !--------------------------------------------------------------------------
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
    
    !--------------------------------------------------------------------------
    ! Print Header (Rank 0 Only)
    !--------------------------------------------------------------------------
    if (rank == 0) then
        print '(A)', '======================================================================'
        print '(A)', 'Io Plasma Torus Hybrid Emission Movie Frame Generator'
        print '(A)', 'MPI-Parallelized Fortran Implementation'
        print '(A)', 'Optical (S+ 6731 A) + UV (S++ 680 A, O+ 833 A)'
        print '(A)', '======================================================================'
        print '(A)', ''
        print '(A,I0,A)', 'Running with ', nprocs, ' MPI processes'
        print '(A)', ''
    end if
    
    !--------------------------------------------------------------------------
    ! Setup Output Directory (Rank 0 Only)
    !--------------------------------------------------------------------------
    output_dir = 'output_data'
    if (rank == 0) then
        print '(A)', 'Setting up output directory...'
        call execute_command_line('mkdir -p ' // trim(output_dir), wait=.true.)
        print '(A,A)', '  Created: ', trim(output_dir)
        print '(A)', ''
    end if
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    
    !--------------------------------------------------------------------------
    ! Create Non-Uniform Grid
    !--------------------------------------------------------------------------
    call create_nonuniform_grids(rho_grid, z_grid, nrho, nz)
    total_points = nrho * nz
    
    if (rank == 0) then
        print '(A)', 'Grid configuration:'
        print '(A,I0,A,I0,A,I0,A)', '  Grid: ', nrho, ' x ', nz, ' = ', total_points, ' points per frame'
        print '(A,F6.1,A,F6.1,A)', '  rho range: ', rho_grid(1), ' to ', rho_grid(nrho), ' R_J'
        print '(A,F6.1,A,F6.1,A)', '  z range: ', z_grid(1), ' to ', z_grid(nz), ' R_J'
        print '(A,I0,A)', '  Total frames: ', N_FRAMES, ' (covering 360 deg rotation)'
        print '(A,F5.2,A)', '  CML step: ', 360.0_real64 / real(N_FRAMES, real64), ' deg per frame'
        print '(A)', ''
    end if
    
    !--------------------------------------------------------------------------
    ! Save Grid Information (Rank 0 Only)
    !--------------------------------------------------------------------------
    if (rank == 0) then
        filename = trim(output_dir) // '/grid_info.bin'
        open(newunit=file_unit, file=filename, form='unformatted', access='stream', status='replace')
        write(file_unit) nrho, nz
        write(file_unit) rho_grid
        write(file_unit) z_grid
        write(file_unit) N_FRAMES
        write(file_unit) S2P_680_VMAX, OP_833_VMAX
        close(file_unit)
        print '(A,A)', 'Saved grid info: ', trim(filename)
        print '(A)', ''
    end if
    
    !--------------------------------------------------------------------------
    ! Create Work Index Array
    !--------------------------------------------------------------------------
    allocate(work_indices(2, total_points))
    idx = 0
    do i = 1, nrho
        do j = 1, nz
            idx = idx + 1
            work_indices(1, idx) = i
            work_indices(2, idx) = j
        end do
    end do
    
    !--------------------------------------------------------------------------
    ! Compute Work Distribution
    !--------------------------------------------------------------------------
    allocate(counts(nprocs), displs(nprocs))
    call compute_work_distribution(total_points, nprocs, counts, displs)
    
    local_start = displs(rank + 1) + 1
    local_end = displs(rank + 1) + counts(rank + 1)
    local_count = counts(rank + 1)
    
    if (rank == 0) then
        print '(A)', 'Work distribution:'
        do i = 1, nprocs
            print '(A,I0,A,I0,A)', '  Rank ', i-1, ': ', counts(i), ' points'
        end do
        print '(A)', ''
    end if
    
    !--------------------------------------------------------------------------
    ! Allocate Local and Global Arrays
    !--------------------------------------------------------------------------
    allocate(local_s1p(local_count))
    allocate(local_s2p(local_count))
    allocate(local_op(local_count))
    
    if (rank == 0) then
        allocate(global_s1p(total_points))
        allocate(global_s2p(total_points))
        allocate(global_op(total_points))
        allocate(s1p_6731_emission(nz, nrho))
        allocate(s2p_680_emission(nz, nrho))
        allocate(op_833_emission(nz, nrho))
    end if
    
    !--------------------------------------------------------------------------
    ! Initialize Raytracer (All Ranks)
    !--------------------------------------------------------------------------
    if (rank == 0) then
        print '(A)', 'Loading plasma model and emission tables...'
    end if
    
    call raytracer%initialize( &
        plasma_file='../../3D_Torus_Model/jovian_plasma_interpolated_381x381x231.h5', &
        emiss_file_single='../../Emiss_tables/CHIANTI_11.0.2_emiss_tables_single_maxwellian_50x50.h5', &
        ierr=load_err)
    
    if (load_err /= 0) then
        if (rank == 0) then
            print '(A)', 'ERROR: Failed to load plasma model or emission tables'
        end if
        call MPI_Finalize(ierr)
        stop 1
    end if
    
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    
    if (rank == 0) then
        print '(A)', 'All ranks initialized successfully.'
        print '(A)', ''
    end if
    
    !--------------------------------------------------------------------------
    ! Prescan for Optical Colorbar Limit
    !--------------------------------------------------------------------------
    if (rank == 0) then
        print '(A)', '======================================================================'
        print '(A,I0,A)', 'Pre-scanning ', N_PRESCAN, ' CML angles for optical colorbar limit...'
        print '(A)', '======================================================================'
    end if
    
    optical_vmax = 0.0_real64
    
    do i = 1, N_PRESCAN
        cml_deg = PRESCAN_ANGLES(i)
        
        if (rank == 0) then
            print '(A,F5.1,A)', '  Scanning CML = ', cml_deg, ' deg...'
        end if
        
        call compute_frame_emission(cml_deg, rho_grid, z_grid, nrho, nz, &
            work_indices, local_start, local_end, local_count, &
            counts, displs, raytracer, &
            local_s1p, local_s2p, local_op, &
            global_s1p, global_s2p, global_op, &
            s1p_6731_emission, s2p_680_emission, op_833_emission, &
            rank, nprocs)
        
        if (rank == 0) then
            prescan_max = maxval(s1p_6731_emission)
            if (prescan_max > optical_vmax) optical_vmax = prescan_max
            print '(A,F8.1,A)', '    S+ 6731 max = ', prescan_max, ' R'
        end if
    end do
    
    ! Apply 6% headroom
    optical_vmax = 1.06_real64 * optical_vmax
    
    ! Broadcast optical_vmax to all ranks
    call MPI_Bcast(optical_vmax, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    
    if (rank == 0) then
        print '(A)', ''
        print '(A)', 'Final colorbar limits:'
        print '(A,F8.1,A)', '  S+ 6731 A:  0.0 - ', optical_vmax, ' R (prescanned with 6% headroom)'
        print '(A,F8.1,A)', '  S++ 680 A:  0.0 - ', S2P_680_VMAX, ' R (manual)'
        print '(A,F8.1,A)', '  O+ 833 A:   0.0 - ', OP_833_VMAX, ' R (manual)'
        print '(A)', ''
        
        ! Save colorbar limits
        filename = trim(output_dir) // '/colorbar_limits.bin'
        open(newunit=file_unit, file=filename, form='unformatted', access='stream', status='replace')
        write(file_unit) optical_vmax, S2P_680_VMAX, OP_833_VMAX
        close(file_unit)
    end if
    
    !--------------------------------------------------------------------------
    ! Generate All Frames
    !--------------------------------------------------------------------------
    if (rank == 0) then
        print '(A)', '======================================================================'
        print '(A,I0,A)', 'Generating ', N_FRAMES, ' frames...'
        print '(A)', '======================================================================'
        print '(A)', ''
        call cpu_time(start_time)
    end if
    
    ! Create CML values array
    allocate(cml_values(N_FRAMES))
    do frame_idx = 1, N_FRAMES
        cml_values(frame_idx) = real(frame_idx - 1, real64) * 360.0_real64 / real(N_FRAMES, real64)
    end do
    
    ! Generate each frame
    do frame_idx = 1, N_FRAMES
        cml_deg = cml_values(frame_idx)
        
        if (rank == 0) then
            call cpu_time(frame_start)
            print '(A,I0,A,I0,A,F6.1,A)', 'Frame ', frame_idx, '/', N_FRAMES, ' | CML = ', cml_deg, ' deg'
        end if
        
        ! Compute emission for this frame
        call compute_frame_emission(cml_deg, rho_grid, z_grid, nrho, nz, &
            work_indices, local_start, local_end, local_count, &
            counts, displs, raytracer, &
            local_s1p, local_s2p, local_op, &
            global_s1p, global_s2p, global_op, &
            s1p_6731_emission, s2p_680_emission, op_833_emission, &
            rank, nprocs)
        
        ! Save frame data (Rank 0 only)
        if (rank == 0) then
            write(frame_str, '(I4.4)') frame_idx - 1
            
            ! Save S+ 6731 data
            filename = trim(output_dir) // '/frame_' // trim(frame_str) // '_s1p_6731.bin'
            open(newunit=file_unit, file=filename, form='unformatted', access='stream', status='replace')
            write(file_unit) s1p_6731_emission
            close(file_unit)
            
            ! Save S++ 680 data
            filename = trim(output_dir) // '/frame_' // trim(frame_str) // '_s2p_680.bin'
            open(newunit=file_unit, file=filename, form='unformatted', access='stream', status='replace')
            write(file_unit) s2p_680_emission
            close(file_unit)
            
            ! Save O+ 833 data
            filename = trim(output_dir) // '/frame_' // trim(frame_str) // '_op_833.bin'
            open(newunit=file_unit, file=filename, form='unformatted', access='stream', status='replace')
            write(file_unit) op_833_emission
            close(file_unit)
            
            ! Print timing info
            call cpu_time(frame_time)
            frame_time = frame_time - frame_start
            call cpu_time(elapsed)
            elapsed = elapsed - start_time
            remaining = real(N_FRAMES - frame_idx, real64) * (elapsed / real(frame_idx, real64))
            
            print '(A,F6.1,A,F6.1,A,F6.1,A)', '  Time: ', frame_time, 's | Elapsed: ', &
                elapsed/60.0_real64, 'min | Remaining: ~', remaining/60.0_real64, 'min'
            print '(A,F8.1,A,F8.1,A,F8.1,A)', '  S+ 6731: max=', maxval(s1p_6731_emission), &
                ' R | S++ 680: max=', maxval(s2p_680_emission), &
                ' R | O+ 833: max=', maxval(op_833_emission), ' R'
        end if
        
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end do
    
    !--------------------------------------------------------------------------
    ! Summary
    !--------------------------------------------------------------------------
    if (rank == 0) then
        call cpu_time(total_time)
        total_time = total_time - start_time
        
        print '(A)', ''
        print '(A)', '======================================================================'
        print '(A)', 'Frame Generation Complete!'
        print '(A)', '======================================================================'
        print '(A)', ''
        print '(A,F8.1,A,F6.2,A)', 'Total computation time: ', total_time/60.0_real64, &
            ' minutes (', total_time/3600.0_real64, ' hours)'
        print '(A,F6.1,A)', 'Average time per frame: ', total_time/real(N_FRAMES, real64), ' seconds'
        print '(A)', ''
        print '(A)', 'Output files saved to: ' // trim(output_dir)
        print '(A,I0,A)', '  ', N_FRAMES * 3, ' binary data files'
        print '(A)', '  grid_info.bin'
        print '(A)', '  colorbar_limits.bin'
        print '(A)', ''
        print '(A)', 'Next step: Run generate_frames_from_binary.py to create PNG frames.'
    end if
    
    !--------------------------------------------------------------------------
    ! Cleanup
    !--------------------------------------------------------------------------
    call raytracer%cleanup()
    
    deallocate(rho_grid, z_grid)
    deallocate(work_indices)
    deallocate(counts, displs)
    deallocate(local_s1p, local_s2p, local_op)
    deallocate(cml_values)
    
    if (rank == 0) then
        deallocate(global_s1p, global_s2p, global_op)
        deallocate(s1p_6731_emission, s2p_680_emission, op_833_emission)
    end if
    
    call MPI_Finalize(ierr)

contains

    !--------------------------------------------------------------------------
    ! Create Non-Uniform Grids
    !--------------------------------------------------------------------------
    subroutine create_nonuniform_grids(rho_grid, z_grid, nrho, nz)
        real(real64), allocatable, intent(out) :: rho_grid(:), z_grid(:)
        integer(int32), intent(out) :: nrho, nz
        
        real(real64), allocatable :: rho_parts(:), z_parts(:)
        integer(int32) :: n_rho_parts, n_z_parts
        integer(int32) :: i, n_seg
        real(real64) :: start_val, stop_val, step_val
        
        ! Radial grid segments: (start, stop, step)
        ! Total ~221 points from -8 to +8 R_J
        ! Array shape: (3, 7) where first index is (start, stop, step), second is segment number
        real(real64), parameter :: RHO_SEGS(3, 7) = reshape([ &
            -8.0_real64, -6.0_real64, 0.1_real64, &
            -6.0_real64, -4.5_real64, 0.025_real64, &
            -4.5_real64, -3.0_real64, 0.1_real64, &
            -3.0_real64, 3.0_real64, 0.2_real64, &
            3.0_real64, 4.5_real64, 0.1_real64, &
            4.5_real64, 6.0_real64, 0.025_real64, &
            6.0_real64, 8.0_real64, 0.1_real64], [3, 7])
        
        ! Vertical grid segments: (start, stop, step)
        ! Total ~131 points from -2 to +2 R_J
        ! Array shape: (3, 3) where first index is (start, stop, step), second is segment number
        real(real64), parameter :: Z_SEGS(3, 3) = reshape([ &
            -2.0_real64, -0.5_real64, 0.1_real64, &
            -0.5_real64, 0.5_real64, 0.01_real64, &
            0.5_real64, 2.0_real64, 0.1_real64], [3, 3])
        
        ! Build rho grid
        call build_grid_from_segments(RHO_SEGS, 7, rho_grid, nrho)
        
        ! Build z grid
        call build_grid_from_segments(Z_SEGS, 3, z_grid, nz)
        
    end subroutine create_nonuniform_grids
    
    
    !--------------------------------------------------------------------------
    ! Build Grid from Segment Specifications
    !--------------------------------------------------------------------------
    subroutine build_grid_from_segments(segs, n_segs, grid, n_points)
        real(real64), intent(in) :: segs(3, *)
        integer(int32), intent(in) :: n_segs
        real(real64), allocatable, intent(out) :: grid(:)
        integer(int32), intent(out) :: n_points
        
        real(real64), allocatable :: temp_grid(:)
        integer(int32) :: i, j, n_in_seg, total_est, idx
        real(real64) :: start_val, stop_val, step_val, val
        
        ! Estimate total points
        total_est = 0
        do i = 1, n_segs
            start_val = segs(1, i)
            stop_val = segs(2, i)
            step_val = segs(3, i)
            total_est = total_est + ceiling((stop_val - start_val) / step_val) + 1
        end do
        
        allocate(temp_grid(total_est))
        
        idx = 0
        do i = 1, n_segs
            start_val = segs(1, i)
            stop_val = segs(2, i)
            step_val = segs(3, i)
            
            val = start_val
            do while (val < stop_val - step_val * 0.5_real64)
                idx = idx + 1
                temp_grid(idx) = val
                val = val + step_val
            end do
            
            ! Include endpoint for last segment
            if (i == n_segs) then
                idx = idx + 1
                temp_grid(idx) = stop_val
            end if
        end do
        
        n_points = idx
        allocate(grid(n_points))
        grid = temp_grid(1:n_points)
        deallocate(temp_grid)
        
    end subroutine build_grid_from_segments
    
    
    !--------------------------------------------------------------------------
    ! Compute Work Distribution
    !--------------------------------------------------------------------------
    subroutine compute_work_distribution(total_work, nprocs, counts, displs)
        integer(int32), intent(in) :: total_work, nprocs
        integer(int32), intent(out) :: counts(:), displs(:)
        
        integer(int32) :: base_count, remainder, i
        
        base_count = total_work / nprocs
        remainder = mod(total_work, nprocs)
        
        displs(1) = 0
        do i = 1, nprocs
            if (i <= remainder) then
                counts(i) = base_count + 1
            else
                counts(i) = base_count
            end if
            if (i < nprocs) displs(i + 1) = displs(i) + counts(i)
        end do
        
    end subroutine compute_work_distribution
    
    
    !--------------------------------------------------------------------------
    ! Compute Ray Geometry
    !--------------------------------------------------------------------------
    subroutine compute_ray_geometry(rho_sky, z_sky, cml_deg, R_obs, slit_pos, norm_vec)
        real(real64), intent(in) :: rho_sky, z_sky, cml_deg, R_obs
        real(real64), intent(out) :: slit_pos(3), norm_vec(3)
        
        real(real64) :: theta
        real(real64), parameter :: DEG_TO_RAD = 3.14159265358979323846_real64 / 180.0_real64
        
        theta = cml_deg * DEG_TO_RAD
        
        ! Viewing direction (toward Jupiter center)
        norm_vec(1) = cos(theta)
        norm_vec(2) = -sin(theta)
        norm_vec(3) = 0.0_real64
        
        ! Ray starting position
        slit_pos(1) = -R_obs * cos(theta) - rho_sky * sin(theta)
        slit_pos(2) =  R_obs * sin(theta) - rho_sky * cos(theta)
        slit_pos(3) = z_sky
        
    end subroutine compute_ray_geometry
    
    
    !--------------------------------------------------------------------------
    ! Calculate Point Emission
    !--------------------------------------------------------------------------
    subroutine calculate_point_emission(raytracer, slit_pos, norm_vec, &
            s1p_val, s2p_val, op_val)
        type(jovian_uv_raytracer), intent(in) :: raytracer
        real(real64), intent(in) :: slit_pos(3), norm_vec(3)
        real(real64), intent(out) :: s1p_val, s2p_val, op_val
        
        real(real64), allocatable :: wave_bins(:), spectrum(:)
        type(emission_line_result), allocatable :: lines(:)
        integer(int32) :: n_bins, n_lines
        real(real64) :: wave_range(2)
        
        ! S+ 6731 Å (Optical)
        wave_range = [S1P_6731_RANGE_MIN, S1P_6731_RANGE_MAX]
        call raytracer%calculate_spectrum_single(slit_pos, norm_vec, wave_range, &
            S1P_6731_BIN_WIDTH, S1P_6731_FWHM, DS, wave_bins, spectrum, lines, n_bins, n_lines)
        
        if (n_bins > 0) then
            s1p_val = simpson_integrate(spectrum, wave_bins, n_bins)
        else
            s1p_val = 0.0_real64
        end if
        
        if (allocated(wave_bins)) deallocate(wave_bins)
        if (allocated(spectrum)) deallocate(spectrum)
        if (allocated(lines)) deallocate(lines)
        
        ! S++ 680 Å (UV)
        wave_range = [S2P_680_RANGE_MIN, S2P_680_RANGE_MAX]
        call raytracer%calculate_spectrum_single(slit_pos, norm_vec, wave_range, &
            S2P_680_BIN_WIDTH, S2P_680_FWHM, DS, wave_bins, spectrum, lines, n_bins, n_lines)
        
        if (n_bins > 0) then
            s2p_val = simpson_integrate(spectrum, wave_bins, n_bins)
        else
            s2p_val = 0.0_real64
        end if
        
        if (allocated(wave_bins)) deallocate(wave_bins)
        if (allocated(spectrum)) deallocate(spectrum)
        if (allocated(lines)) deallocate(lines)
        
        ! O+ 833 Å (UV)
        wave_range = [OP_833_RANGE_MIN, OP_833_RANGE_MAX]
        call raytracer%calculate_spectrum_single(slit_pos, norm_vec, wave_range, &
            OP_833_BIN_WIDTH, OP_833_FWHM, DS, wave_bins, spectrum, lines, n_bins, n_lines)
        
        if (n_bins > 0) then
            op_val = simpson_integrate(spectrum, wave_bins, n_bins)
        else
            op_val = 0.0_real64
        end if
        
        if (allocated(wave_bins)) deallocate(wave_bins)
        if (allocated(spectrum)) deallocate(spectrum)
        if (allocated(lines)) deallocate(lines)
        
    end subroutine calculate_point_emission
    
    
    !--------------------------------------------------------------------------
    ! Compute Frame Emission (MPI Parallel)
    !--------------------------------------------------------------------------
    subroutine compute_frame_emission(cml_deg, rho_grid, z_grid, nrho, nz, &
            work_indices, local_start, local_end, local_count, &
            counts, displs, raytracer, &
            local_s1p, local_s2p, local_op, &
            global_s1p, global_s2p, global_op, &
            s1p_6731_emission, s2p_680_emission, op_833_emission, &
            rank, nprocs)
        
        real(real64), intent(in) :: cml_deg
        real(real64), intent(in) :: rho_grid(:), z_grid(:)
        integer(int32), intent(in) :: nrho, nz
        integer(int32), intent(in) :: work_indices(:,:)
        integer(int32), intent(in) :: local_start, local_end, local_count
        integer(int32), intent(in) :: counts(:), displs(:)
        type(jovian_uv_raytracer), intent(in) :: raytracer
        real(real64), intent(out) :: local_s1p(:), local_s2p(:), local_op(:)
        real(real64), intent(out) :: global_s1p(:), global_s2p(:), global_op(:)
        real(real64), intent(out) :: s1p_6731_emission(:,:), s2p_680_emission(:,:), op_833_emission(:,:)
        integer, intent(in) :: rank, nprocs
        
        integer(int32) :: k, idx, i, j, local_idx, ierr
        real(real64) :: rho_sky, z_sky
        real(real64) :: slit_pos(3), norm_vec(3)
        real(real64) :: s1p_val, s2p_val, op_val
        
        ! Each rank computes its assigned points
        local_idx = 0
        do k = local_start, local_end
            local_idx = local_idx + 1
            i = work_indices(1, k)
            j = work_indices(2, k)
            
            rho_sky = rho_grid(i)
            z_sky = z_grid(j)
            
            call compute_ray_geometry(rho_sky, z_sky, cml_deg, R_OBS, slit_pos, norm_vec)
            call calculate_point_emission(raytracer, slit_pos, norm_vec, s1p_val, s2p_val, op_val)
            
            local_s1p(local_idx) = s1p_val
            local_s2p(local_idx) = s2p_val
            local_op(local_idx) = op_val
        end do
        
        ! Gather results to rank 0
        call MPI_Gatherv(local_s1p, local_count, MPI_DOUBLE_PRECISION, &
            global_s1p, counts, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Gatherv(local_s2p, local_count, MPI_DOUBLE_PRECISION, &
            global_s2p, counts, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Gatherv(local_op, local_count, MPI_DOUBLE_PRECISION, &
            global_op, counts, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        
        ! Rank 0 reshapes into 2D arrays
        if (rank == 0) then
            do k = 1, nrho * nz
                i = work_indices(1, k)
                j = work_indices(2, k)
                s1p_6731_emission(j, i) = global_s1p(k)
                s2p_680_emission(j, i) = global_s2p(k)
                op_833_emission(j, i) = global_op(k)
            end do
        end if
        
    end subroutine compute_frame_emission

end program make_emission_movie_frames_mpi
