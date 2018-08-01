module MOM_energy

use MOM_coms, only : num_PEs, root_PE, PE_here
use MOM_EOS, only : calculate_density
use MOM_grid, only : ocean_grid_type
use MOM_error_handler, only : is_root_pe, MOM_error, FATAL
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

! use raw MPI for its easier collectives (particularly reduce and exscan)
include "mpif.h"

type, public :: rpe_elem
  real :: density !< Density of a cell (sort key)
  real :: volume  !< Volume of a cell (needs to be in sorted order to calculate position)
end type rpe_elem

type, public :: rpe_CS
  integer :: rpe_type !< MPI type for contiguous reals, corresponding to rpe_elem
  integer :: np !< num_PEs
end type rpe_CS

public calculate_RPE_init, calculate_RPE

contains

  !> Partition a list around its last value
  function qksrt_partition(n, list, start, end) result(top)
    integer,        intent(in)    :: n
    type(rpe_elem), intent(inout) :: list(n)
    integer,        intent(in)    :: start, end

    ! Local variables
    type(rpe_elem) :: pivot
    integer        :: bottom, top
    logical        :: done

    pivot = list(end)                                 ! Partition around the last value
    bottom = start-1                                  ! Start outside the area to be partitioned
    top = end                                         ! Ditto

    done = .false.
    do while (.not. done)                             ! Until all elements are partitioned...

      do while (.not. done)                           ! Until we find an out of place element...
        bottom = bottom+1                             ! ... move the bottom up.

        if(bottom == top) then                        ! If we hit the top...
          done = .true.                               ! ... we are done.
          exit
        endif

        if(list(bottom)%density > pivot%density) then ! Is the bottom out of place?
          list(top) = list(bottom)                    ! Then put it at the top...
          exit                                        ! ... and start searching from the top.
        endif
      enddo

      do while (.not. done)                           ! Until we find an out of place element...
        top = top-1                                   ! ... move the top down.

        if(top == bottom) then                        ! If we hit the bottom...
          done = .true.                               ! ... we are done.
          exit
        endif

        if(list(top)%density < pivot%density) then    ! Is the top out of place?
          list(bottom) = list(top)                    ! Then put it at the bottom...
          exit                                        ! ...and start searching from the bottom.
        endif
      enddo
    enddo

    list(top) = pivot                              ! Put the pivot in its place.
    ! Return the split point
  end function qksrt_partition

  !> Recursively sort the subarray of list from start to end (inclusive)
  recursive subroutine quicksort(n, list, start, end)
    implicit none
    integer,        intent(in)    :: n
    type(rpe_elem), intent(inout) :: list(n)
    integer,        intent(in)    :: start, end

    integer :: split ! Split point between the sublists (from partition)

    if(start < end) then                            ! If there are two or more elements...
      split = qksrt_partition(n, list, start, end)  ! ... partition the sublist...
      call quicksort(n, list,  start, split-1)      ! ... and sort both halves.
      call quicksort(n, list, split+1, end)
    endif
  end subroutine quicksort

  !> Wrapper function to quicksort an array
  subroutine sort(A, n)
    integer,        intent(in)    :: n
    type(rpe_elem), intent(inout) :: A(n)

    call quicksort(n, A, 1, n)
  end subroutine sort

  subroutine merge(n, size, A, counts, B)
    integer, intent(in) :: n, size
    type(rpe_elem), intent(in) :: A(size)
    type(rpe_elem), intent(inout) :: B(size)
    integer, intent(in) :: counts(n)

    ! pointer into each bucket
    integer, dimension(n) :: pos, max_pos
    ! loop indices
    integer :: i, j, smallest

    ! calculate initial pointers
    pos(1) = 1
    do j = 2, n
      ! index of each bucket into A
      pos(j) = pos(j-1) + counts(j-1)
      ! previous bucket can't overrun this one
      max_pos(j-1) = pos(j)
    enddo
    max_pos(n) = size + 1

    ! merge loop over entire bucket
    do i = 1, size
      ! initially no bucket is a candidate
      smallest = -1

      do j = 1, n
        ! for each bucket, check that we haven't hit the end
        if (pos(j) >= max_pos(j)) cycle

        ! set smallest if it's not set now
        if (smallest == -1) smallest = j

        ! whether we're smaller than the minimum
        if (A(pos(j))%density < A(pos(smallest))%density) smallest = j
      enddo

      ! unable to merge?
      if (smallest == -1) call MOM_error(FATAL, "merge: couldn't find smallest element")

      ! assign element
      B(i) = A(pos(smallest))
      ! increase pointer for bucket we used
      pos(smallest) = pos(smallest) + 1
    enddo
  end subroutine merge

  subroutine calculate_RPE_init(CS)
    type(rpe_CS), pointer :: CS

    if (associated(CS)) call MOM_error(FATAL, "calculate_RPE_init: CS already associated!")

    allocate(CS)
    CS%rpe_type = MPI_2DOUBLE_PRECISION
    CS%np = num_PEs()
  end subroutine calculate_RPE_init

  subroutine calculate_RPE(CS, G, GV, h, tv, RPE)
    type(rpe_CS),                             intent(in)  :: CS  !< RPE control structure
    type(ocean_grid_type),                    intent(in)  :: G   !< Ocean (horizontal) grid
    type(verticalGrid_type),                  intent(in)  :: GV  !< Ocean (vertical) grid
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)  :: h   !< Thicknesses
    type(thermo_var_ptrs),                    intent(in)  :: tv  !< Thermodynamics
    real,                                     intent(out) :: RPE !< Globally integrated RPE

    ! Local variables
    type(rpe_elem), target :: cells_full(G%isc:G%iec,G%jsc:G%jec,G%ke) ! RPE elements only in computational domain
    type(rpe_elem), pointer :: cells(:)
    type(rpe_elem) :: pivots(CS%np), all_pivots(CS%np*CS%np) ! pivot elements for bucket sort
    type(rpe_elem), allocatable :: bucket(:), bucket_sorted(:)

    integer :: counts(CS%np), recv_counts(CS%np) ! local element counts for each bucket
    integer :: displs(CS%np), rdispls(CS%np) ! displacements into bucket to receive

    real :: p_ref(SZK_(G)) ! Reference pressure for every cell in a column
    real :: density(SZI_(G),SZJ_(G),SZK_(G))
    real, allocatable :: interfaces(:) ! Depth of local interfaces

    real :: local_depth, total_thickness, depth, local_RPE

    integer :: i, j, k ! loop indices
    integer :: n, np, nk ! total number of cells, PEs and layers
    integer :: bucket_size, ierr

    n = size(cells_full)
    nk = GV%ke
    np = CS%np

    p_ref(:) = tv%P_Ref ! use same reference pressure for all cells
    cells(1:n) => cells_full ! flattened version of cells

    ! Calculate density columnwise and populate cells
    do j = G%jsc,G%jec ; do i = G%isc,G%iec
      call calculate_density(tv%T(i,j,:), tv%S(i,j,:), p_ref, density(i,j,:), 1, nk, tv%eqn_of_state)

      do k = 1, nk
        cells_full(i,j,k)%density = density(i,j,k)
        cells_full(i,j,k)%volume  = h(i,j,k) * G%areaT(i,j)
      enddo
    enddo; enddo

    ! Sort local cells
    call sort(cells, n)

    ! Choose np-1 pivots
    do i = 1, np - 1
      pivots(i) = cells(i * n / np)
    enddo

    ! Gather pivots on root
    call mpi_gather(pivots, np - 1, CS%rpe_type, all_pivots, np - 1, CS%rpe_type, root_PE(), MPI_COMM_WORLD, ierr)

    ! Sort on root and choose global set of pivots
    if (is_root_pe()) then
      call sort(all_pivots, np * (np - 1))
      do i = 1, np - 1
        pivots(i) = all_pivots(i * (np - 1))
      enddo
    endif

    ! Broadcast global pivots
    call mpi_bcast(pivots, np - 1, CS%rpe_type, root_PE(), MPI_COMM_WORLD, ierr)

    ! Compute counts for each bucket
    k = 1
    j = 1
    do i = 1, np - 1
      ! advance pointer until out of range for current bucket
      do while (j < n .and. cells(j)%density <= pivots(i)%density)
        j = j + 1
      enddo

      counts(i) = j - k
      k = j
    enddo
    ! update last bucket count
    counts(np) = n - k + 1

    ! Let all PEs know what to expect from all other PEs
    call mpi_alltoall(counts, 1, MPI_INTEGER, recv_counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    ! Allocate receive buffers
    bucket_size = sum(recv_counts)
    allocate(bucket(bucket_size))
    allocate(bucket_sorted(bucket_size))
    allocate(interfaces(0:bucket_size))

    ! Calculate displacements into buckets array to receive elements
    displs(1) = 0 ; rdispls(1) = 0
    do i = 1, np - 1
      displs(i+1) = displs(i) + counts(i)
      rdispls(i+1) = rdispls(i) + recv_counts(i)
    enddo

    ! Distribute data to correct buckets
    call mpi_alltoallv(cells, counts, displs, CS%rpe_type, bucket, recv_counts, rdispls, CS%rpe_type, &
         MPI_COMM_WORLD, ierr)

    ! Merge sorted subarrays
    call merge(np, bucket_size, bucket, recv_counts, bucket_sorted)

    ! Calculate total thickness of this bucket
    total_thickness = 0.
    interfaces(0) = 0.
    do i = 1, bucket_size
      interfaces(i) = interfaces(i-1) + bucket_sorted(i)%volume / G%areaT_global
      total_thickness = total_thickness + bucket_sorted(i)%volume / G%areaT_global
    enddo

    ! Scan thickness for a running sum
    call mpi_exscan(total_thickness, local_depth, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Calculate depths of local cells and integrate into bucket RPE
    local_RPE = 0.
    do i = 1, bucket_size
      depth = local_depth + (interfaces(i-1) + interfaces(i)) / 2.
      local_RPE = local_RPE + bucket_sorted(i)%density * bucket_sorted(i)%volume * depth
    enddo
    local_RPE = local_RPE * G%g_Earth / G%areaT_global

    ! Reduce RPE to a total global sum on all PEs
    call mpi_allreduce(local_RPE, RPE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    deallocate(interfaces)
    deallocate(bucket)
    deallocate(bucket_sorted)
  end subroutine calculate_RPE
end module MOM_energy
