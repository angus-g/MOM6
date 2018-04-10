!> Regrid columns for the adaptive coordinate
module coord_adapt

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_EOS,           only : calculate_density_derivs
use MOM_error_handler, only : MOM_error, FATAL
use MOM_variables,     only : ocean_grid_type, thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type
use coord_zlike, only : init_coord_zlike, end_coord_zlike, zlike_CS

implicit none ; private

#include <MOM_memory.h>

type, public :: adapt_diag_CS
  real, dimension(:,:,:), pointer :: phys_u => null()
  real, dimension(:,:,:), pointer :: phys_v => null()

  real, dimension(:,:,:), pointer :: slope_u => null()
  real, dimension(:,:,:), pointer :: slope_v => null()

  real, dimension(:,:,:), pointer :: denom_u => null()
  real, dimension(:,:,:), pointer :: denom_v => null()

  real, dimension(:,:,:), pointer :: coord_u => null()
  real, dimension(:,:,:), pointer :: coord_v => null()

  real, dimension(:,:,:), pointer :: limiting_smooth => null()
  real, dimension(:,:,:), pointer :: limiting_dense => null()

  real, dimension(:,:,:), pointer :: dk_sig => null()
end type adapt_diag_CS

type, public :: adapt_CS
  !private

  !> Number of layers/levels
  integer :: nk

  !> Nominal near-surface resolution
  real, allocatable, dimension(:) :: coordinateResolution

  !> Density adaptivity coefficient
  real :: adaptAlphaRho = -1.0

  !> Pressure adaptivity coefficient
  real :: adaptAlphaP = -1.0

  !> Timescale for diffusivity
  real :: adaptTimescale, restoringTimescale

  !> Coordinate relaxation coefficient/timescale
  real :: adaptTau = 0.0

  real :: adaptCutoff, adaptSmooth

  logical :: mean_h = .false.
  logical :: twin_grad = .true.
  logical :: physicalSlope
  logical :: restoreMean

  type(zlike_CS), pointer :: zlike_CS => null()
end type adapt_CS

public init_coord_adapt, set_adapt_params, build_adapt_column, end_coord_adapt

contains

!> Initialise an adapt_CS with parameters
subroutine init_coord_adapt(CS, nk, coordinateResolution)
  type(adapt_CS),     pointer    :: CS !< Unassociated pointer to hold the control structure
  integer,            intent(in) :: nk
  real, dimension(:), intent(in) :: coordinateResolution

  if (associated(CS)) call MOM_error(FATAL, "init_coord_adapt: CS already associated")
  allocate(CS)
  allocate(CS%coordinateResolution(nk))

  CS%nk = nk
  CS%coordinateResolution(:) = coordinateResolution(:)

  call init_coord_zlike(CS%zlike_CS, nk, coordinateResolution)
end subroutine init_coord_adapt

!> Clean up the coordinate control structure
subroutine end_coord_adapt(CS)
  type(adapt_CS), pointer :: CS

  ! nothing to do
  if (.not. associated(CS)) return

  call end_coord_zlike(CS%zlike_CS)

  deallocate(CS%coordinateResolution)
  deallocate(CS)
end subroutine end_coord_adapt

subroutine set_adapt_params(CS, adaptAlphaRho, adaptAlphaP, adaptTimescale, adaptTau, adaptMean, adaptTwin, adaptCutoff, adaptSmooth, adaptPhysicalSlope, restoringTimescale, restoreMean)
  type(adapt_CS),    pointer    :: CS
  real, optional,    intent(in) :: adaptAlphaRho, adaptAlphaP, adaptTimescale, adaptTau, &
       adaptCutoff, adaptSmooth, restoringTimescale
  logical, optional, intent(in) :: adaptMean, adaptTwin, adaptPhysicalSlope, restoreMean

  if (.not. associated(CS)) call MOM_error(FATAL, "set_adapt_params: CS not associated")

  if (present(adaptAlphaRho)) CS%adaptAlphaRho = adaptAlphaRho
  if (present(adaptAlphaP)) CS%adaptAlphaP = adaptAlphaP

  if (present(adaptTimescale)) CS%adaptTimescale = adaptTimescale
  if (present(restoringTimescale)) CS%restoringTimescale = restoringTimescale

  if (present(adaptTau)) CS%adaptTau = adaptTau
  if (present(adaptMean)) CS%mean_h = adaptMean
  if (present(adaptTwin)) CS%twin_grad = adaptTwin
  if (present(adaptCutoff)) CS%adaptCutoff = adaptCutoff
  if (present(adaptSmooth)) CS%adaptSmooth = adaptSmooth
  if (present(adaptPhysicalSlope)) CS%physicalSlope = adaptPhysicalSlope
  if (present(restoreMean)) CS%restoreMean = restoreMean
end subroutine set_adapt_params

subroutine build_adapt_column(CS, G, GV, tv, i, j)
  type(adapt_CS),                              intent(in)    :: CS
  type(ocean_grid_type),                       intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),                     intent(in)    :: GV   !< The ocean's vertical grid structure
  type(thermo_var_ptrs),                       intent(in)    :: tv   !< A structure pointing to various thermodynamic variables
  integer,                                     intent(in)    :: i, j
end subroutine build_adapt_column

end module coord_adapt
