module user_change_diffusivity

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_error_handler, only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_variables, only : thermo_var_ptrs, vertvisc_type, p3d
use MOM_EOS, only : calculate_density
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public user_change_diff, user_change_diff_init
public user_change_diff_end

type, public :: user_change_diff_CS ; private
  real :: Kd_factor
  real :: decay_dist
  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
end type user_change_diff_CS

contains

!> This subroutine provides an interface for a user to use to modify the
!! main code to alter the diffusivities as needed.  The specific example
!! implemented here augments the diffusivity for a specified range of latitude
!! and coordinate potential density.
subroutine user_change_diff(h, tv, G, GV, CS, Kd, Kd_int, T_f, S_f, Kd_int_add)
  type(ocean_grid_type),                    intent(in)    :: G   !< The ocean's grid structure.
  type(verticalGrid_type),   intent(in)    :: GV   !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h   !< Layer thickness, in m or kg m-2.
  type(thermo_var_ptrs),                    intent(in)    :: tv  !< A structure containing pointers
                                                                 !! to any available thermodynamic
                                                                 !! fields. Absent fields have NULL ptrs.
  type(user_change_diff_CS),                pointer       :: CS  !< This module's control structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   optional, intent(inout) :: Kd !< The diapycnal diffusivity of
                                                                  !! each layer in m2 s-1.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), optional, intent(inout) :: Kd_int !< The diapycnal diffusivity
                                                                  !! at each interface in m2 s-1.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   optional, intent(in)    :: T_f !< Temperature with massless
                                                                  !! layers filled in vertically.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   optional, intent(in)    :: S_f !< Salinity with massless
                                                                  !! layers filled in vertically.
  real, dimension(:,:,:),                     optional, pointer       :: Kd_int_add !< The diapycnal
                                                                  !! diffusivity that is being added at
                                                                  !! each interface in m2 s-1.

  logical :: store_Kd_add  ! Save the added diffusivity as a diagnostic if true.
  integer :: i, j, k, is, ie, js, je, nz
  integer :: isd, ied, jsd, jed

  real :: Z_lay, Z_int
  character(len=200) :: mesg

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (.not.associated(CS)) call MOM_error(FATAL,"user_set_diffusivity: "//&
         "Module must be initialized before it is used.")

  store_Kd_add = .false.
  if (present(Kd_int_add)) store_Kd_add = associated(Kd_int_add)
  if (store_Kd_add) Kd_int_add(:,:,:) = 0.0

  do j=js,je
    if (present(Kd)) then
      do i=is,ie
        Z_lay = -G%bathyT(i,j) + 0.5*GV%H_to_m*h(i,j,nz)
        do k=nz,1,-1
          Kd(i,j,k) = Kd(i,j,k) * (1.0 + CS%Kd_factor * exp(-(Z_lay / CS%decay_dist) ** 2))
          Z_lay = Z_lay + 0.5*GV%H_to_m*(h(i,j,k) + h(i,j,k-1))
        enddo
      enddo
    endif

    if (present(Kd_int)) then
      do i=is,ie
        Z_int = -G%bathyT(i,j)
        do K=nz,1,-1
          Z_int = Z_int + GV%H_to_m*h(i,j,k)
          Kd_int(i,j,K) = Kd_int(i,j,K) * (1.0 + CS%Kd_factor * exp(-(Z_int / CS%decay_dist) ** 2))
          if (store_Kd_add) Kd_int_add(i,j,k) = CS%Kd_factor * exp(-(Z_int / CS%decay_dist) ** 2)
        enddo
      enddo
    endif
  enddo
end subroutine user_change_diff

!> Set up the module control structure.
subroutine user_change_diff_init(Time, G, param_file, diag, CS)
  type(time_type),           intent(in)    :: Time       !< The current model time.
  type(ocean_grid_type),     intent(in)    :: G          !< The ocean's grid structure.
  type(param_file_type),     intent(in)    :: param_file !< A structure indicating the
                                                         !! open file to parse for
                                                         !! model parameter values.
  type(diag_ctrl), target,   intent(inout) :: diag       !< A structure that is used to
                                                         !! regulate diagnostic output.
  type(user_change_diff_CS), pointer       :: CS         !< A pointer that is set to
                                                         !! point to the control
                                                         !! structure for this module.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "user_set_diffusivity"  ! This module's name.
  character(len=200) :: mesg

  if (associated(CS)) then
    call MOM_error(WARNING, "use_change_diff_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "USER_KD_FACTOR", CS%Kd_factor, &
       "The factor by which KD is elevated at the surface,\n"//&
       "subtracted by one.", default=0.0)
  if (CS%Kd_factor /= 0.0) then
    call get_param(param_file, mdl, "USER_KD_DECAY_DIST", CS%decay_dist, &
         "The decay distance away from the surface for elevated\n"//&
         "diffusivity.", units="m", default=1.0e2)
  endif

end subroutine user_change_diff_init

!> Clean up the module control structure.
subroutine user_change_diff_end(CS)
  type(user_change_diff_CS), pointer :: CS         !< A pointer that is set to
                                                   !! point to the control
                                                   !! structure for this module.

  if (associated(CS)) deallocate(CS)

end subroutine user_change_diff_end

!> \namespace user_change_diffusivity
!!
!!  By Robert Hallberg, May 2012
!!
!!    This file contains a subroutine that increments the diapycnal
!!  diffusivity in a specified band of latitudes and densities.
!!
!!     A small fragment of the grid is shown below:
!!
!!    j+1  x ^ x ^ x   At x:  q
!!    j+1  > o > o >   At ^:  v
!!    j    x ^ x ^ x   At >:  u
!!    j    > o > o >   At o:  h, T, S, Kd, etc.
!!    j-1  x ^ x ^ x
!!        i-1  i  i+1  At x & ^:
!!           i  i+1    At > & o:
!!
!!  The boundaries always run through q grid points (x).

end module user_change_diffusivity
