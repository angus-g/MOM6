!> This module contains I/O framework code
module MOM_io

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_domains,          only : MOM_domain_type, AGRID, BGRID_NE, CGRID_NE, get_domain_components
use MOM_domains,          only : domain1D, get_simple_array_i_ind, get_simple_array_j_ind
use MOM_dyn_horgrid,      only : dyn_horgrid_type
use MOM_error_handler,    only : MOM_error, NOTE, FATAL, WARNING
use MOM_file_parser,      only : log_version, param_file_type
use MOM_grid,             only : ocean_grid_type
use MOM_io_wrapper,       only : MOM_read_data, MOM_read_vector, MOM_write_field, read_axis_data
use MOM_io_wrapper,       only : file_exists, field_exists, read_field_chksum
use MOM_io_wrapper,       only : open_file, close_file, field_size, fieldtype, get_filename_appendix
use MOM_io_wrapper,       only : flush_file, get_file_info, get_file_atts, get_file_fields
use MOM_io_wrapper,       only : get_file_times, read_data, axistype, get_axis_data
use MOM_io_wrapper,       only : write_field, write_metadata, write_version_number, get_ensemble_id
use MOM_io_wrapper,       only : open_namelist_file, check_nml_error, io_infra_init, io_infra_end
use MOM_io_wrapper,       only : APPEND_FILE, ASCII_FILE, MULTIPLE, NETCDF_FILE, OVERWRITE_FILE
use MOM_io_wrapper,       only : READONLY_FILE, SINGLE_FILE, WRITEONLY_FILE
use MOM_io_wrapper,       only : CENTER, CORNER, NORTH_FACE, EAST_FACE
use MOM_string_functions, only : lowercase, slasher
use MOM_verticalGrid,     only : verticalGrid_type

use iso_fortran_env,      only : stdout_iso=>output_unit, stderr_iso=>error_unit
use netcdf,               only : NF90_open, NF90_inquire, NF90_inq_varids, NF90_inquire_variable
use netcdf,               only : NF90_Inquire_Dimension, NF90_max_name, NF90_max_var_dims
use netcdf,               only : NF90_STRERROR, NF90_NOWRITE, NF90_NOERR

implicit none ; private

! These interfaces are actually implemented in this file.
public :: create_file, reopen_file, num_timelevels, cmor_long_std, ensembler, MOM_io_init
public :: var_desc, modify_vardesc, query_vardesc
! The following are simple pass throughs of routines from MOM_io_wrapper or other modules
public :: close_file, field_exists, field_size, fieldtype, get_filename_appendix
public :: file_exists, flush_file, get_file_info, get_file_atts, get_file_fields
public :: get_file_times, open_file, read_axis_data, read_data, read_field_chksum
public :: MOM_read_data, MOM_read_vector, MOM_write_field, get_axis_data
public :: slasher, write_field, write_version_number
public :: open_namelist_file, check_nml_error, io_infra_init, io_infra_end
! These are encoding constants.
public :: APPEND_FILE, ASCII_FILE, MULTIPLE, NETCDF_FILE, OVERWRITE_FILE
public :: READONLY_FILE, SINGLE_FILE, WRITEONLY_FILE
public :: CENTER, CORNER, NORTH_FACE, EAST_FACE

!> Type for describing a variable, typically a tracer
type, public :: vardesc
  character(len=64)  :: name               !< Variable name in a NetCDF file
  character(len=48)  :: units              !< Physical dimensions of the variable
  character(len=240) :: longname           !< Long name of the variable
  character(len=8)   :: hor_grid           !< Horizontal grid:  u, v, h, q, Cu, Cv, T, Bu, or 1
  character(len=8)   :: z_grid             !< Vertical grid:  L, i, or 1
  character(len=8)   :: t_grid             !< Time description: s, p, or 1
  character(len=64)  :: cmor_field_name    !< CMOR name
  character(len=64)  :: cmor_units         !< CMOR physical dimensions of the variable
  character(len=240) :: cmor_longname      !< CMOR long name of the variable
  real               :: conversion         !< for unit conversions, such as needed to
                                           !! convert from intensive to extensive
end type vardesc

integer, public :: stdout = stdout_iso  !< standard output unit
integer, public :: stderr = stderr_iso  !< standard output unit

contains

!> Routine creates a new NetCDF file.  It also sets up
!! structures that describe this file and variables that will
!! later be written to this file. Type for describing a variable, typically a tracer
subroutine create_file(unit, filename, vars, novars, fields, threading, timeunit, G, dG, GV, checksums)
  integer,               intent(out)   :: unit       !< unit id of an open file or -1 on a
                                                     !! nonwriting PE with single file output
  character(len=*),      intent(in)    :: filename   !< full path to the file to create
  type(vardesc),         intent(in)    :: vars(:)    !< structures describing fields written to filename
  integer,               intent(in)    :: novars     !< number of fields written to filename
  type(fieldtype),       intent(inout) :: fields(:)  !< array of fieldtypes for each variable
  integer, optional,     intent(in)    :: threading  !< SINGLE_FILE or MULTIPLE
  real, optional,        intent(in)    :: timeunit   !< length of the units for time [s]. The
                                                     !! default value is 86400.0, for 1 day.
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if the new file uses any
                                                     !! vertical grid axes.
  integer(kind=8), optional,      intent(in)    :: checksums(:,:)  !< checksums of vars

  logical        :: use_lath, use_lonh, use_latq, use_lonq, use_time
  logical        :: use_layer, use_int, use_periodic
  logical        :: one_file, domain_set
  type(axistype) :: axis_lath, axis_latq, axis_lonh, axis_lonq
  type(axistype) :: axis_layer, axis_int, axis_time, axis_periodic
  type(axistype) :: axes(4)
  type(MOM_domain_type), pointer :: Domain => NULL()
  type(domain1d) :: x_domain, y_domain
  integer        :: numaxes, pack, thread, k
  integer        :: isg, ieg, jsg, jeg, IsgB, IegB, JsgB, JegB
  integer        :: var_periods, num_periods=0
  real, dimension(:), allocatable :: period_val
  real, pointer, dimension(:) :: &
    gridLatT => NULL(), & ! The latitude or longitude of T or B points for
    gridLatB => NULL(), & ! the purpose of labeling the output axes.
    gridLonT => NULL(), gridLonB => NULL()
  character(len=40) :: time_units, x_axis_units, y_axis_units
  character(len=8)  :: t_grid, t_grid_read

  use_lath  = .false. ; use_lonh     = .false.
  use_latq  = .false. ; use_lonq     = .false.
  use_time  = .false. ; use_periodic = .false.
  use_layer = .false. ; use_int      = .false.

  thread = SINGLE_FILE
  if (PRESENT(threading)) thread = threading

  domain_set = .false.
  if (present(G)) then
    domain_set = .true. ; Domain => G%Domain
    gridLatT => G%gridLatT ; gridLatB => G%gridLatB
    gridLonT => G%gridLonT ; gridLonB => G%gridLonB
    x_axis_units = G%x_axis_units ; y_axis_units = G%y_axis_units
    isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg
    IsgB = G%IsgB ; IegB = G%IegB ; JsgB = G%JsgB ; JegB = G%JegB
  elseif (present(dG)) then
    domain_set = .true. ; Domain => dG%Domain
    gridLatT => dG%gridLatT ; gridLatB => dG%gridLatB
    gridLonT => dG%gridLonT ; gridLonB => dG%gridLonB
    x_axis_units = dG%x_axis_units ; y_axis_units = dG%y_axis_units
    isg = dG%isg ; ieg = dG%ieg ; jsg = dG%jsg ; jeg = dG%jeg
    IsgB = dG%IsgB ; IegB = dG%IegB ; JsgB = dG%JsgB ; JegB = dG%JegB
  endif

  one_file = .true.
  if (domain_set) one_file = (thread == SINGLE_FILE)

  if (one_file) then
    call open_file(unit, filename, OVERWRITE_FILE, NETCDF_FILE, threading=thread)
  else
    call open_file(unit, filename, OVERWRITE_FILE, NETCDF_FILE, MOM_domain=Domain)
  endif

! Define the coordinates.
  do k=1,novars
    select case (vars(k)%hor_grid)
      case ('h') ; use_lath = .true. ; use_lonh = .true.
      case ('q') ; use_latq = .true. ; use_lonq = .true.
      case ('u') ; use_lath = .true. ; use_lonq = .true.
      case ('v') ; use_latq = .true. ; use_lonh = .true.
      case ('T')  ; use_lath = .true. ; use_lonh = .true.
      case ('Bu') ; use_latq = .true. ; use_lonq = .true.
      case ('Cu') ; use_lath = .true. ; use_lonq = .true.
      case ('Cv') ; use_latq = .true. ; use_lonh = .true.
      case ('1') ! Do nothing.
      case default
        call MOM_error(WARNING, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized hor_grid "//trim(vars(k)%hor_grid))
    end select
    select case (vars(k)%z_grid)
      case ('L') ; use_layer = .true.
      case ('i') ; use_int = .true.
      case ('1') ! Do nothing.
      case default
        call MOM_error(FATAL, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized z_grid "//trim(vars(k)%z_grid))
    end select
    t_grid = adjustl(vars(k)%t_grid)
    select case (t_grid(1:1))
      case ('s', 'a', 'm') ; use_time = .true.
      case ('p') ; use_periodic = .true.
        if (len_trim(t_grid(2:8)) <= 0) call MOM_error(FATAL, &
          "MOM_io create_file: No periodic axis length was specified in "//&
          trim(vars(k)%t_grid) // " in the periodic axes of variable "//&
          trim(vars(k)%name)//" in file "//trim(filename))
        var_periods = -9999999
        t_grid_read = adjustl(t_grid(2:8))
        read(t_grid_read,*) var_periods
        if (var_periods == -9999999) call MOM_error(FATAL, &
          "MOM_io create_file: Failed to read the number of periods from "//&
          trim(vars(k)%t_grid) // " in the periodic axes of variable "//&
          trim(vars(k)%name)//" in file "//trim(filename))
        if (var_periods < 1) call MOM_error(FATAL, "MOM_io create_file: "//&
           "variable "//trim(vars(k)%name)//" in file "//trim(filename)//&
           " uses a periodic time axis, and must have a positive "//&
           "value for the number of periods in "//vars(k)%t_grid )
        if ((num_periods > 0) .and. (var_periods /= num_periods)) &
          call MOM_error(FATAL, "MOM_io create_file: "//&
            "Only one value of the number of periods can be used in the "//&
            "create_file call for file "//trim(filename)//".  The second is "//&
            "variable "//trim(vars(k)%name)//" with t_grid "//vars(k)%t_grid )

        num_periods = var_periods
      case ('1') ! Do nothing.
      case default
        call MOM_error(WARNING, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized t_grid "//trim(vars(k)%t_grid))
    end select
  enddo

  if ((use_lath .or. use_lonh .or. use_latq .or. use_lonq)) then
    if (.not.domain_set) call MOM_error(FATAL, "create_file: "//&
      "An ocean_grid_type or dyn_horgrid_type is required to create a file with a horizontal coordinate.")

    call get_domain_components(Domain, x_domain, y_domain)
  endif
  if ((use_layer .or. use_int) .and. .not.present(GV)) call MOM_error(FATAL, &
    "create_file: A vertical grid type is required to create a file with a vertical coordinate.")

  if (use_lath) &
    call write_metadata(unit, axis_lath, name="lath", units=y_axis_units, longname="Latitude", &
                        cartesian='Y', domain=y_domain, data=gridLatT(jsg:jeg))

  if (use_lonh) &
    call write_metadata(unit, axis_lonh, name="lonh", units=x_axis_units, longname="Longitude", &
                        cartesian='X', domain=x_domain, data=gridLonT(isg:ieg))

  if (use_latq) &
    call write_metadata(unit, axis_latq, name="latq", units=y_axis_units, longname="Latitude", &
                        cartesian='Y', domain=y_domain, data=gridLatB(JsgB:JegB))

  if (use_lonq) &
    call write_metadata(unit, axis_lonq, name="lonq", units=x_axis_units, longname="Longitude", &
                        cartesian='X', domain=x_domain, data=gridLonB(IsgB:IegB))

  if (use_layer) &
    call write_metadata(unit, axis_layer, name="Layer", units=trim(GV%zAxisUnits), &
                        longname="Layer "//trim(GV%zAxisLongName), cartesian='Z', &
                        sense=1, data=GV%sLayer(1:GV%ke))

  if (use_int) &
    call write_metadata(unit, axis_int, name="Interface", units=trim(GV%zAxisUnits), &
                        longname="Interface "//trim(GV%zAxisLongName), cartesian='Z', &
                        sense=1, data=GV%sInterface(1:GV%ke+1))

  if (use_time) then ; if (present(timeunit)) then
    ! Set appropriate units, depending on the value.
    if (timeunit < 0.0) then
      time_units = "days" ! The default value.
    elseif ((timeunit >= 0.99) .and. (timeunit < 1.01)) then
      time_units = "seconds"
    elseif ((timeunit >= 3599.0) .and. (timeunit < 3601.0)) then
      time_units = "hours"
    elseif ((timeunit >= 86399.0) .and. (timeunit < 86401.0)) then
      time_units = "days"
    elseif ((timeunit >= 3.0e7) .and. (timeunit < 3.2e7)) then
      time_units = "years"
    else
      write(time_units,'(es8.2," s")') timeunit
    endif

    call write_metadata(unit, axis_time, name="Time", units=time_units, longname="Time", cartesian='T')
  else
    call write_metadata(unit, axis_time, name="Time", units="days", longname="Time", cartesian= 'T')
  endif ; endif

  if (use_periodic) then
    if (num_periods <= 1) call MOM_error(FATAL, "MOM_io create_file: "//&
      "num_periods for file "//trim(filename)//" must be at least 1.")
    ! Define a periodic axis with unit labels.
    allocate(period_val(num_periods))
    do k=1,num_periods ; period_val(k) = real(k) ; enddo
    call write_metadata(unit, axis_periodic, name="Period", units="nondimensional", &
                        longname="Periods for cyclical varaiables", cartesian='T', data=period_val)
    deallocate(period_val)
  endif

  do k=1,novars
    numaxes = 0
    select case (vars(k)%hor_grid)
      case ('h')  ; numaxes = 2 ; axes(1) = axis_lonh ; axes(2) = axis_lath
      case ('q')  ; numaxes = 2 ; axes(1) = axis_lonq ; axes(2) = axis_latq
      case ('u')  ; numaxes = 2 ; axes(1) = axis_lonq ; axes(2) = axis_lath
      case ('v')  ; numaxes = 2 ; axes(1) = axis_lonh ; axes(2) = axis_latq
      case ('T')  ; numaxes = 2 ; axes(1) = axis_lonh ; axes(2) = axis_lath
      case ('Bu') ; numaxes = 2 ; axes(1) = axis_lonq ; axes(2) = axis_latq
      case ('Cu') ; numaxes = 2 ; axes(1) = axis_lonq ; axes(2) = axis_lath
      case ('Cv') ; numaxes = 2 ; axes(1) = axis_lonh ; axes(2) = axis_latq
      case ('1') ! Do nothing.
      case default
        call MOM_error(WARNING, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized hor_grid "//trim(vars(k)%hor_grid))
    end select
    select case (vars(k)%z_grid)
      case ('L') ; numaxes = numaxes+1 ; axes(numaxes) = axis_layer
      case ('i') ; numaxes = numaxes+1 ; axes(numaxes) = axis_int
      case ('1') ! Do nothing.
      case default
        call MOM_error(FATAL, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized z_grid "//trim(vars(k)%z_grid))
    end select
    t_grid = adjustl(vars(k)%t_grid)
    select case (t_grid(1:1))
      case ('s', 'a', 'm') ; numaxes = numaxes+1 ; axes(numaxes) = axis_time
      case ('p')           ; numaxes = numaxes+1 ; axes(numaxes) = axis_periodic
      case ('1') ! Do nothing.
      case default
        call MOM_error(WARNING, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized t_grid "//trim(vars(k)%t_grid))
    end select

    pack = 1
    if (present(checksums)) then
      call write_metadata(unit, fields(k), axes(1:numaxes), vars(k)%name, vars(k)%units, &
                          vars(k)%longname, pack=pack, checksum=checksums(k,:))
    else
      call write_metadata(unit, fields(k), axes(1:numaxes), vars(k)%name, vars(k)%units, &
                          vars(k)%longname, pack=pack)
    endif
  enddo

  if (use_lath) call write_field(unit, axis_lath)
  if (use_latq) call write_field(unit, axis_latq)
  if (use_lonh) call write_field(unit, axis_lonh)
  if (use_lonq) call write_field(unit, axis_lonq)
  if (use_layer) call write_field(unit, axis_layer)
  if (use_int) call write_field(unit, axis_int)
  if (use_periodic) call write_field(unit, axis_periodic)

end subroutine create_file


!> This routine opens an existing NetCDF file for output.  If it
!! does not find the file, a new file is created.  It also sets up
!! structures that describe this file and the variables that will
!! later be written to this file.
subroutine reopen_file(unit, filename, vars, novars, fields, threading, timeunit, G, dG, GV)
  integer,               intent(out)   :: unit       !< unit id of an open file or -1 on a
                                                     !! nonwriting PE with single file output
  character(len=*),      intent(in)    :: filename   !< full path to the file to create
  type(vardesc),         intent(in)    :: vars(:)    !< structures describing fields written to filename
  integer,               intent(in)    :: novars     !< number of fields written to filename
  type(fieldtype),       intent(inout) :: fields(:)  !< array of fieldtypes for each variable
  integer, optional,     intent(in)    :: threading  !< SINGLE_FILE or MULTIPLE
  real, optional,        intent(in)    :: timeunit   !< length of the units for time [s]. The
                                                     !! default value is 86400.0, for 1 day.
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if a new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if a new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if a new file uses any
                                                     !! vertical grid axes.

  type(MOM_domain_type), pointer :: Domain => NULL()
  character(len=200) :: check_name, mesg
  integer :: length, ndim, nvar, natt, ntime, thread
  logical :: exists, one_file, domain_set

  thread = SINGLE_FILE
  if (PRESENT(threading)) thread = threading

  check_name = filename
  length = len(trim(check_name))
  if (check_name(length-2:length) /= ".nc") check_name = trim(check_name)//".nc"
  if (thread /= SINGLE_FILE) check_name = trim(check_name)//".0000"

  inquire(file=check_name,EXIST=exists)

  if (.not.exists) then
    call create_file(unit, filename, vars, novars, fields, threading, timeunit, &
                     G=G, dG=dG, GV=GV)
  else

    domain_set = .false.
    if (present(G)) then
      domain_set = .true. ; Domain => G%Domain
    elseif (present(dG)) then
      domain_set = .true. ; Domain => dG%Domain
    endif

    one_file = .true.
    if (domain_set) one_file = (thread == SINGLE_FILE)

    if (one_file) then
      call open_file(unit, filename, APPEND_FILE, NETCDF_FILE, threading=thread)
    else
      call open_file(unit, filename, APPEND_FILE, NETCDF_FILE, MOM_domain=Domain)
    endif
    if (unit < 0) return

    call get_file_info(unit, ndim, nvar, natt, ntime)

    if (nvar == -1) then
      write (mesg,*) "Reopening file ",trim(filename)," apparently had ",nvar,&
                     " variables. Clobbering and creating file with ",novars," instead."
      call MOM_error(WARNING,"MOM_io: "//mesg)
      call create_file(unit, filename, vars, novars, fields, threading, timeunit, G=G, GV=GV)
    elseif (nvar /= novars) then
      write (mesg,*) "Reopening file ",trim(filename)," with ",novars,&
                     " variables instead of ",nvar,"."
      call MOM_error(FATAL,"MOM_io: "//mesg)
    endif

    if (nvar > 0) call get_file_fields(unit, fields(1:nvar))

    ! Check for inconsistent field names...
!    do i=1,nvar
!      call get_field_atts(fields(i), name)
!      !if (trim(name) /= trim(vars%name)) then
!      !  write (mesg, '("Reopening file ",a," variable ",a," is called ",a,".")',&
!      !         trim(filename), trim(vars%name), trim(name))
!      !  call MOM_error(NOTE, "MOM_io: "//trim(mesg))
!      !endif
!    enddo
  endif

end subroutine reopen_file


!> This function determines how many time levels a variable has.
function num_timelevels(filename, varname, min_dims) result(n_time)
  character(len=*),  intent(in) :: filename   !< name of the file to read
  character(len=*),  intent(in) :: varname    !< variable whose number of time levels
                                              !! are to be returned
  integer, optional, intent(in) :: min_dims   !< The minimum number of dimensions a variable must have
                                              !! if it has a time dimension.  If the variable has 1 less
                                              !! dimension than this, then 0 is returned.
  integer :: n_time                           !< number of time levels varname has in filename

  logical :: found
  character(len=200) :: msg
  character(len=nf90_max_name) :: name
  integer :: ncid, nvars, status, varid, ndims, n
  integer, allocatable :: varids(:)
  integer, dimension(nf90_max_var_dims) :: dimids

  n_time = -1
  found = .false.

  status = NF90_OPEN(filename, NF90_NOWRITE, ncid)
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING,"num_timelevels: "//&
        " Difficulties opening "//trim(filename)//" - "//trim(NF90_STRERROR(status)))
    return
  endif

  status = NF90_INQUIRE(ncid, nVariables=nvars)
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING,"num_timelevels: "//&
        " Difficulties getting the number of variables in file "//&
        trim(filename)//" - "//trim(NF90_STRERROR(status)))
    return
  endif

  if (nvars < 1) then
    call MOM_error(WARNING,"num_timelevels: "//&
        " There appear not to be any variables in "//trim(filename))
    return
  endif

  allocate(varids(nvars))

  status = nf90_inq_varids(ncid, nvars, varids)
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING,"num_timelevels: "//&
        " Difficulties getting the variable IDs in file "//&
        trim(filename)//" - "//trim(NF90_STRERROR(status)))
    deallocate(varids) ; return
  endif

  do n = 1,nvars
    status = nf90_inquire_variable(ncid, varids(n), name=name)
    if (status /= NF90_NOERR) then
      call MOM_error(WARNING,"num_timelevels: "//&
          " Difficulties getting a variable name in file "//&
          trim(filename)//" - "//trim(NF90_STRERROR(status)))
    endif

    if (trim(lowercase(name)) == trim(lowercase(varname))) then
      if (found) then
        call MOM_error(WARNING,"num_timelevels: "//&
          " Two variables match the case-insensitive name "//trim(varname)//&
          " in file "//trim(filename)//" - "//trim(NF90_STRERROR(status)))
      else
        varid = varids(n) ; found = .true.
      endif
    endif
  enddo

  deallocate(varids)

  if (.not.found) then
    call MOM_error(WARNING,"num_timelevels: "//&
        " variable "//trim(varname)//" was not found in file "//trim(filename))
    return
  endif

  status = nf90_inquire_variable(ncid, varid, ndims = ndims)
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING,"num_timelevels: "//trim(NF90_STRERROR(status))//&
      " Getting number of dimensions of "//trim(varname)//" in "//trim(filename))
    return
  endif

  if (present(min_dims)) then
    if (ndims < min_dims-1) then
      write(msg, '(I3)') min_dims
      call MOM_error(WARNING, "num_timelevels: variable "//trim(varname)//&
        " in file "//trim(filename)//" has fewer than min_dims = "//trim(msg)//&
        " dimensions.")
    elseif (ndims == min_dims - 1) then
      n_time = 0 ; return
    endif
  endif

  status = nf90_inquire_variable(ncid, varid, dimids = dimids(1:ndims))
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING,"num_timelevels: "//trim(NF90_STRERROR(status))//&
      " Getting last dimension ID for "//trim(varname)//" in "//trim(filename))
    return
  endif

  status = nf90_Inquire_Dimension(ncid, dimids(ndims), len=n_time)
  if (status /= NF90_NOERR) call MOM_error(WARNING,"num_timelevels: "//&
      trim(NF90_STRERROR(status))//" Getting number of time levels of "//&
      trim(varname)//" in "//trim(filename))

end function num_timelevels


!> Returns a vardesc type whose elements have been filled with the provided
!! fields.  The argument name is required, while the others are optional and
!! have default values that are empty strings or are appropriate for a 3-d
!! tracer field at the tracer cell centers.
function var_desc(name, units, longname, hor_grid, z_grid, t_grid, &
                  cmor_field_name, cmor_units, cmor_longname, conversion, caller) result(vd)
  character(len=*),           intent(in) :: name               !< variable name
  character(len=*), optional, intent(in) :: units              !< variable units
  character(len=*), optional, intent(in) :: longname           !< variable long name
  character(len=*), optional, intent(in) :: hor_grid           !< variable horizontal staggering
  character(len=*), optional, intent(in) :: z_grid             !< variable vertical staggering
  character(len=*), optional, intent(in) :: t_grid             !< time description: s, p, or 1
  character(len=*), optional, intent(in) :: cmor_field_name    !< CMOR name
  character(len=*), optional, intent(in) :: cmor_units         !< CMOR physical dimensions of variable
  character(len=*), optional, intent(in) :: cmor_longname      !< CMOR long name
  real            , optional, intent(in) :: conversion         !< for unit conversions, such as needed to
                                                               !! convert from intensive to extensive
  character(len=*), optional, intent(in) :: caller             !< calling routine?
  type(vardesc)                          :: vd                 !< vardesc type that is created

  character(len=120) :: cllr
  cllr = "var_desc"
  if (present(caller)) cllr = trim(caller)

  call safe_string_copy(name, vd%name, "vd%name", cllr)

  vd%longname = "" ; vd%units = ""
  vd%hor_grid = 'h' ; vd%z_grid = 'L' ; vd%t_grid = 's'

  vd%cmor_field_name  =  ""
  vd%cmor_units       =  ""
  vd%cmor_longname    =  ""
  vd%conversion       =  1.0

  call modify_vardesc(vd, units=units, longname=longname, hor_grid=hor_grid, &
                      z_grid=z_grid, t_grid=t_grid,                          &
                      cmor_field_name=cmor_field_name,cmor_units=cmor_units, &
                      cmor_longname=cmor_longname, conversion=conversion, caller=cllr)

end function var_desc


!> This routine modifies the named elements of a vardesc type.
!! All arguments are optional, except the vardesc type to be modified.
subroutine modify_vardesc(vd, name, units, longname, hor_grid, z_grid, t_grid, &
                 cmor_field_name, cmor_units, cmor_longname, conversion, caller)
  type(vardesc),              intent(inout) :: vd              !< vardesc type that is modified
  character(len=*), optional, intent(in)    :: name            !< name of variable
  character(len=*), optional, intent(in)    :: units           !< units of variable
  character(len=*), optional, intent(in)    :: longname        !< long name of variable
  character(len=*), optional, intent(in)    :: hor_grid        !< horizontal staggering of variable
  character(len=*), optional, intent(in)    :: z_grid          !< vertical staggering of variable
  character(len=*), optional, intent(in)    :: t_grid          !< time description: s, p, or 1
  character(len=*), optional, intent(in)    :: cmor_field_name !< CMOR name
  character(len=*), optional, intent(in)    :: cmor_units      !< CMOR physical dimensions of variable
  character(len=*), optional, intent(in)    :: cmor_longname   !< CMOR long name
  real            , optional, intent(in)    :: conversion      !< for unit conversions, such as needed
                                                               !! to convert from intensive to extensive
  character(len=*), optional, intent(in)    :: caller          !< calling routine?

  character(len=120) :: cllr
  cllr = "mod_vardesc"
  if (present(caller)) cllr = trim(caller)

  if (present(name))      call safe_string_copy(name, vd%name, "vd%name", cllr)

  if (present(longname))  call safe_string_copy(longname, vd%longname, &
                               "vd%longname of "//trim(vd%name), cllr)
  if (present(units))     call safe_string_copy(units, vd%units,       &
                               "vd%units of "//trim(vd%name), cllr)
  if (present(hor_grid))  call safe_string_copy(hor_grid, vd%hor_grid, &
                               "vd%hor_grid of "//trim(vd%name), cllr)
  if (present(z_grid))    call safe_string_copy(z_grid, vd%z_grid,     &
                               "vd%z_grid of "//trim(vd%name), cllr)
  if (present(t_grid))    call safe_string_copy(t_grid, vd%t_grid,     &
                               "vd%t_grid of "//trim(vd%name), cllr)

  if (present(cmor_field_name)) call safe_string_copy(cmor_field_name, vd%cmor_field_name, &
                                     "vd%cmor_field_name of "//trim(vd%name), cllr)
  if (present(cmor_units))      call safe_string_copy(cmor_units, vd%cmor_units, &
                                     "vd%cmor_units of "//trim(vd%name), cllr)
  if (present(cmor_longname))   call safe_string_copy(cmor_longname, vd%cmor_longname, &
                                     "vd%cmor_longname of "//trim(vd%name), cllr)

end subroutine modify_vardesc

!> This function returns the CMOR standard name given a CMOR longname, based on
!! the standard pattern of character conversions.
function cmor_long_std(longname) result(std_name)
  character(len=*), intent(in) :: longname  !< The CMOR longname being converted
  character(len=len(longname)) :: std_name  !< The CMOR standard name generated from longname

  integer :: k

  std_name = lowercase(longname)

  do k=1, len_trim(std_name)
    if (std_name(k:k) == ' ') std_name(k:k) = '_'
  enddo

end function cmor_long_std

!> This routine queries vardesc
subroutine query_vardesc(vd, name, units, longname, hor_grid, z_grid, t_grid, &
                         cmor_field_name, cmor_units, cmor_longname, conversion, caller)
  type(vardesc),              intent(in)  :: vd                 !< vardesc type that is queried
  character(len=*), optional, intent(out) :: name               !< name of variable
  character(len=*), optional, intent(out) :: units              !< units of variable
  character(len=*), optional, intent(out) :: longname           !< long name of variable
  character(len=*), optional, intent(out) :: hor_grid           !< horizontal staggering of variable
  character(len=*), optional, intent(out) :: z_grid             !< verticle staggering of variable
  character(len=*), optional, intent(out) :: t_grid             !< time description: s, p, or 1
  character(len=*), optional, intent(out) :: cmor_field_name    !< CMOR name
  character(len=*), optional, intent(out) :: cmor_units         !< CMOR physical dimensions of variable
  character(len=*), optional, intent(out) :: cmor_longname      !< CMOR long name
  real            , optional, intent(out) :: conversion         !< for unit conversions, such as needed to
                                                                !! convert from intensive to extensive
  character(len=*), optional, intent(in)  :: caller             !< calling routine?


  character(len=120) :: cllr
  cllr = "mod_vardesc"
  if (present(caller)) cllr = trim(caller)

  if (present(name))      call safe_string_copy(vd%name, name,         &
                               "vd%name of "//trim(vd%name), cllr)
  if (present(longname))  call safe_string_copy(vd%longname, longname, &
                               "vd%longname of "//trim(vd%name), cllr)
  if (present(units))     call safe_string_copy(vd%units, units,       &
                               "vd%units of "//trim(vd%name), cllr)
  if (present(hor_grid))  call safe_string_copy(vd%hor_grid, hor_grid, &
                               "vd%hor_grid of "//trim(vd%name), cllr)
  if (present(z_grid))    call safe_string_copy(vd%z_grid, z_grid,     &
                               "vd%z_grid of "//trim(vd%name), cllr)
  if (present(t_grid))    call safe_string_copy(vd%t_grid, t_grid,     &
                               "vd%t_grid of "//trim(vd%name), cllr)

  if (present(cmor_field_name)) call safe_string_copy(vd%cmor_field_name, cmor_field_name, &
                                     "vd%cmor_field_name of "//trim(vd%name), cllr)
  if (present(cmor_units))      call safe_string_copy(vd%cmor_units, cmor_units,          &
                                     "vd%cmor_units of "//trim(vd%name), cllr)
  if (present(cmor_longname))   call safe_string_copy(vd%cmor_longname, cmor_longname, &
                                     "vd%cmor_longname of "//trim(vd%name), cllr)

end subroutine query_vardesc

!> Copies a string
subroutine safe_string_copy(str1, str2, fieldnm, caller)
  character(len=*),           intent(in)  :: str1    !< The string being copied
  character(len=*),           intent(out) :: str2    !< The string being copied into
  character(len=*), optional, intent(in)  :: fieldnm !< The name of the field for error messages
  character(len=*), optional, intent(in)  :: caller  !< The calling routine for error messages

  if (len(trim(str1)) > len(str2)) then
    if (present(fieldnm) .and. present(caller)) then
      call MOM_error(FATAL, trim(caller)//" attempted to copy the overly long"//&
        " string "//trim(str1)//" into "//trim(fieldnm))
    else
      call MOM_error(FATAL, "safe_string_copy: The string "//trim(str1)//&
                     " is longer than its intended target.")
    endif
  endif
  str2 = trim(str1)
end subroutine safe_string_copy

!> Returns a name with "%#E" or "%E" replaced with the ensemble member number.
function ensembler(name, ens_no_in) result(en_nm)
  character(len=*),  intent(in) :: name       !< The name to be modified
  integer, optional, intent(in) :: ens_no_in  !< The number of the current ensemble member
  character(len=len(name)) :: en_nm  !< The name encoded with the ensemble number

  ! This function replaces "%#E" or "%E" with the ensemble number anywhere it
  ! occurs in name, with %E using 4 or 6 digits (depending on the ensemble size)
  ! and %#E using # digits, where # is a number from 1 to 9.

  character(len=len(name)) :: tmp
  character(10) :: ens_num_char
  character(3)  :: code_str
  integer :: ens_no
  integer :: n, is, ie

  en_nm = trim(name)
  if (index(name,"%") == 0) return

  if (present(ens_no_in)) then
    ens_no = ens_no_in
  else
    ens_no = get_ensemble_id()
  endif

  write(ens_num_char, '(I10)') ens_no ; ens_num_char = adjustl(ens_num_char)
  do
    is = index(en_nm,"%E")
    if (is == 0) exit
    if (len(en_nm) < len(trim(en_nm)) + len(trim(ens_num_char)) - 2) &
      call MOM_error(FATAL, "MOM_io ensembler: name "//trim(name)// &
      " is not long enough for %E expansion for ens_no "//trim(ens_num_char))
    tmp = en_nm(1:is-1)//trim(ens_num_char)//trim(en_nm(is+2:))
    en_nm = tmp
  enddo

  if (index(name,"%") == 0) return

  write(ens_num_char, '(I10.10)') ens_no
  do n=1,9 ; do
    write(code_str, '("%",I1,"E")') n

    is = index(en_nm,code_str)
    if (is == 0) exit
    if (ens_no < 10**n) then
      if (len(en_nm) < len(trim(en_nm)) + n-3) call MOM_error(FATAL, &
        "MOM_io ensembler: name "//trim(name)//" is not long enough for %E expansion.")
      tmp = en_nm(1:is-1)//trim(ens_num_char(11-n:10))//trim(en_nm(is+3:))
    else
      call MOM_error(FATAL, "MOM_io ensembler: Ensemble number is too large "//&
          "to be encoded with "//code_str//" in "//trim(name))
    endif
    en_nm = tmp
  enddo ; enddo

end function ensembler

!> Initialize the MOM_io module
subroutine MOM_io_init(param_file)
  type(param_file_type), intent(in) :: param_file  !< structure indicating the open file to
                                                   !! parse for model parameter values.

  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_io" ! This module's name.

  call log_version(param_file, mdl, version)

end subroutine MOM_io_init

!> \namespace mom_io
!!
!!   This file contains a number of subroutines that manipulate
!!  NetCDF files and handle input and output of fields.  These
!!  subroutines, along with their purpose, are:
!!
!!   * create_file: create a new file and set up structures that are
!!       needed for subsequent output and write out the coordinates.
!!   * reopen_file: reopen an existing file for writing and set up
!!       structures that are needed for subsequent output.
!!   * open_input_file: open the indicated file for reading only.
!!   * close_file: close an open file.
!!   * synch_file: flush the buffers, completing all pending output.
!!
!!   * write_field: write a field to an open file.
!!   * write_time: write a value of the time axis to an open file.
!!   * read_data: read a variable from an open file.
!!   * read_time: read a time from an open file.
!!
!!   * name_output_file: provide a name for an output file based on a
!!       name root and the time of the output.
!!   * find_input_file: find a file that has been previously written by
!!       MOM and named by name_output_file and open it for reading.
!!
!!   * handle_error: write an error code and quit.



end module MOM_io
