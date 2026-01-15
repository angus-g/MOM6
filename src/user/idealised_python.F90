module idealised_python

use, intrinsic :: iso_c_binding

use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : FATAL, MOM_error
use MOM_file_parser, only : get_param, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public idealised_python_topography, idealised_python_TS, idealised_python_velocity

character(len=len("idealised_python")) :: mdl = "idealised_python"
logical :: python_initialised = .false.

type(c_ptr), bind(C, name="PyObject_None") :: PyObject_None
integer(kind=c_int), bind(C, name="_NPY_FLOAT64") :: NPY_FLOAT64
integer(kind=c_int), bind(C, name="_NPY_ARRAY_F_CONTIGUOUS") :: NPY_ARRAY_F_CONTIGUOUS

interface
   function array_shim() bind(C)
     import :: c_int

     integer(kind=c_int) :: array_shim
   end function
end interface

interface
   function PyArray_DATA(arr) bind(C, name="_PyArray_DATA")
     import :: c_ptr

     type(c_ptr), value, intent(in) :: arr
     type(c_ptr) :: PyArray_DATA
   end function
end interface

interface
   function PyArray_DescrFromType(typenum) bind(C, name="_PyArray_DescrFromType")
     import :: c_int, c_ptr

     integer(kind=c_int), value, intent(in) :: typenum
     type(c_ptr) :: PyArray_DescrFromType
   end function PyArray_DescrFromType
end interface

interface
   function PyArray_FromAny(op, dtype, min_depth, max_depth, requirements, context) bind(C, name="_PyArray_FromAny")
     import :: c_int, c_ptr

     type(c_ptr), value, intent(in) :: op, dtype, context
     integer(kind=c_int), value, intent(in) :: min_depth, max_depth, requirements
     type(c_ptr) :: PyArray_FromAny
   end function PyArray_FromAny
end interface

interface
   subroutine Py_DECREF(obj) bind(C, name="Py_DecRef")
     import :: c_ptr

     type(c_ptr), value, intent(in) :: obj
   end subroutine Py_DECREF
end interface

interface
   function PyErr_Occurred() bind(C, name="PyErr_Occurred")
     import :: c_ptr

     type(c_ptr) :: PyErr_Occurred
   end function PyErr_Occurred
end interface

interface
   subroutine PyErr_Print() bind(C, name="PyErr_Print")
   end subroutine PyErr_Print
end interface

interface
   function PyFloat_FromDouble(val) bind(C, name="PyFloat_FromDouble")
     import :: c_double, c_ptr

     real(kind=c_double), value, intent(in) :: val
     type(c_ptr) :: PyFloat_FromDouble
  end function PyFloat_FromDouble
end interface

interface
   function PyImport_Import(name) bind(C, name="PyImport_Import")
     import :: c_ptr

     type(c_ptr), value, intent(in) :: name
     type(c_ptr) :: PyImport_Import
   end function PyImport_Import
end interface

interface
   subroutine Py_INCREF(obj) bind(C, name="Py_IncRef")
     import :: c_ptr

     type(c_ptr), value, intent(in) :: obj
   end subroutine Py_INCREF
end interface

interface
   subroutine Py_Initialize() bind(C, name="Py_Initialize")
   end subroutine Py_Initialize
end interface

interface
   function PyList_Insert(list, index, obj) bind(C, name="PyList_Insert")
     import :: c_int, c_ptr, c_size_t

     type(c_ptr), value, intent(in) :: list
     integer(kind=c_size_t), value, intent(in) :: index
     type(c_ptr), value, intent(in) :: obj
     integer(kind=c_int) :: PyList_Insert
   end function PyList_Insert
end interface

interface
   function PyLong_FromLong(val) bind(C, name="PyLong_FromLong")
     import :: c_long, c_ptr

     integer(kind=c_long), value, intent(in) :: val
     type(c_ptr) :: PyLong_FromLong
   end function PyLong_FromLong
end interface

interface
   function PyObject_VectorcallMethod(name, args, nargs, kwnames) bind(C, name="PyObject_VectorcallMethod")
     import :: c_ptr, c_size_t

     type(c_ptr), value, intent(in) :: name
     integer(kind=c_size_t), value, intent(in) :: nargs
     type(c_ptr), dimension(nargs), intent(in) :: args
     type(c_ptr), value, intent(in) :: kwnames
     type(c_ptr) :: PyObject_VectorcallMethod
   end function PyObject_VectorcallMethod
end interface

interface
   function PySys_GetObject(name) bind(C, name="PySys_GetObject")
     import :: c_char, c_ptr

     character(kind=c_char), intent(in) :: name(*)
     type(c_ptr) :: PySys_GetObject
   end function PySys_GetObject
end interface

interface
   function PyTuple_GetItem(p, pos) bind(C, name="PyTuple_GetItem")
     import :: c_ptr, c_size_t

     type(c_ptr), value, intent(in) :: p
     integer(kind=c_size_t), value, intent(in) :: pos
     type(c_ptr) :: PyTuple_GetItem
   end function PyTuple_GetItem
end interface

interface
   function PyTuple_New(len) bind(C, name="PyTuple_New")
     import :: c_ptr, c_size_t

     integer(kind=c_size_t), value, intent(in) :: len
     type(c_ptr) :: PyTuple_New
   end function PyTuple_New
end interface

interface
   function PyTuple_SetItem(p, pos, o) bind(C, name="PyTuple_SetItem")
     import :: c_int, c_ptr, c_size_t

     type(c_ptr), value, intent(in) :: p, o
     integer(kind=c_size_t), value, intent(in) :: pos
     integer(kind=c_int) :: PyTuple_SetItem
   end function PyTuple_SetItem
end interface

interface
   function PyUnicode_FromString(str) bind(C, name="PyUnicode_FromString")
     import :: c_char, c_ptr

     character(kind=c_char), intent(in) :: str(*)
     type(c_ptr) :: PyUnicode_FromString
   end function PyUnicode_FromString
end interface

contains

subroutine idealised_python_topography(D, G, param_file, max_depth)
  type(dyn_horgrid_type), intent(in) :: G
  real, dimension(G%isd:G%ied,G%jsd:G%jed), intent(out) :: D
  type(param_file_type), intent(in) :: param_file
  real, intent(in) :: max_depth

  character(len=20) :: module_name, topo_func
  type(c_ptr) :: topo_mod

  call python_init

  call get_param(param_file, mdl, "PYTHON_TOPO_MODULE", module_name, &
       "Python module containing topography definition function", fail_if_missing=.true.)
  call get_param(param_file, mdl, "PYTHON_TOPO_FUNC", topo_func, &
       "Python function f(isd, ied, jsd, jed, max_depth) returning\n"//&
       "an array of the topography depths.", fail_if_missing=.true.)

  call load_module(module_name, topo_mod)
  call run_topo_func(topo_mod, topo_func, G, max_depth, D)
  call unload_module(topo_mod)
end subroutine idealised_python_topography

subroutine idealised_python_velocity(u, v, G, GV, param_file, just_read)
  type(ocean_grid_type), intent(in) :: G
  type(verticalGrid_type), intent(in) :: GV
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(out) :: u
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(out) :: v
  type(param_file_type), intent(in) :: param_file
  logical, intent(in) :: just_read

  character(len=20) :: module_name, velocity_func
  type(c_ptr) :: velocity_mod

  if (just_read) return

  call python_init

  call get_param(param_file, mdl, "PYTHON_VELOCITY_MODULE", module_name, &
       "Python module containing velocity definition function", fail_if_missing=.true.)
  call get_param(param_file, mdl, "PYTHON_VELOCITY_FUNC", velocity_func, &
       "Python function f((IsdB,IedB,jsd,jed),(isd,ied,JsdB,JedB),nz) returning\n"//&
       "(u,v) arrays of the velocity.", fail_if_missing=.true.)

  call load_module(module_name, velocity_mod)
  call run_velocity_func(velocity_mod, velocity_func, G, GV, u, v)
  call unload_module(velocity_mod)
end subroutine idealised_python_velocity

subroutine idealised_python_TS(T, S, G, GV, param_file, just_read)
  type(ocean_grid_type), intent(in) :: G
  type(verticalGrid_type), intent(in) :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T, S
  type(param_file_type), intent(in) :: param_file
  logical, intent(in) :: just_read

  character(len=20) :: module_name, ts_func
  type(c_ptr) :: ts_mod

  if (just_read) return

  call python_init

  call get_param(param_file, mdl, "PYTHON_TS_MODULE", module_name, &
       "Python module containing temperature/salinity definition function", fail_if_missing=.true.)
  call get_param(param_file, mdl, "PYTHON_TS_FUNC", ts_func, &
       "Python function f(isd, ied, jsd, jed) returning\n"//&
       "(T,S) arrays of the temperature/salinity.", fail_if_missing=.true.)

  call load_module(module_name, ts_mod)
  call run_ts_func(ts_mod, ts_func, G, GV, T, S)
  call unload_module(ts_mod)
end subroutine idealised_python_TS

subroutine python_init()
  integer(kind=c_int) :: ret
  integer(kind=c_size_t) :: index = 0
  type(c_ptr) :: empty_str, sys_path

  if (.not. python_initialised) then
    call Py_Initialize()
    sys_path = PySys_GetObject("path" // c_null_char)
    if (.not. c_associated(sys_path)) call MOM_error(FATAL, "sys.path not associated")
    call Py_INCREF(sys_path)
    empty_str = PyUnicode_FromString(c_null_char)
    if (.not. c_associated(empty_str)) call MOM_error(FATAL, "error creating empty string")
    ret = PyList_Insert(sys_path, index, empty_str)

    if (ret /= 0) then
      call PyErr_Print
      call MOM_error(FATAL, "unable to modify sys.path")
    end if

    call Py_DECREF(empty_str)
    call Py_DECREF(sys_path)

    ret = array_shim()
    if (ret < 0) then
      call MOM_error(FATAL, "unable to import numpy api")
    end if

    python_initialised = .true.
  end if
end subroutine python_init

function py_string(name)
  character(len=*), intent(in) :: name
  type(c_ptr) :: py_string

  character(kind=c_char) :: cname(len_trim(name) + 1)
  integer :: i, lv
  type(c_ptr) :: err

  lv = len_trim(name)
  do concurrent (i=1:lv)
    cname(i) = name(i:i)
  end do
  cname(lv+1) = c_null_char

  py_string = PyUnicode_FromString(cname)
  err = PyErr_Occurred()
  if (c_associated(err)) then
    call PyErr_Print
    call MOM_error(FATAL, "exception creating python string")
  end if
end function py_string

subroutine load_module(name, topo_mod)
  character(len=*), intent(in) :: name
  type(c_ptr), intent(out) :: topo_mod
  type(c_ptr) :: py_name, err

  py_name = py_string(name)
  topo_mod = PyImport_Import(py_name)
  err = PyErr_Occurred()
  if (c_associated(err)) then
    call PyErr_Print
    call MOM_error(FATAL, "exception loading module " // name)
  end if

  call Py_DECREF(py_name)
end subroutine load_module

subroutine unload_module(py_mod)
  type(c_ptr), intent(in) :: py_mod

  call Py_DECREF(py_mod)
end subroutine unload_module

subroutine run_topo_func(topo_mod, name, G, max_depth, D)
  type(c_ptr), intent(in) :: topo_mod
  character(len=*), intent(in) :: name
  type(dyn_horgrid_type), intent(in) :: G
  real, intent(in) :: max_depth
  real, dimension(G%isd:G%ied,G%jsd:G%jed), intent(out) :: D

  integer :: i
  type(c_ptr) :: ret, err, arr, method_name, arrptr
  type(c_ptr), dimension(:), allocatable :: args
  real, dimension(:,:), pointer :: Dptr

  allocate(args(6))
  args(1) = topo_mod
  args(2) = PyLong_FromLong(int(G%isd, kind=c_long))
  args(3) = PyLong_FromLong(int(G%ied, kind=c_long))
  args(4) = PyLong_FromLong(int(G%jsd, kind=c_long))
  args(5) = PyLong_FromLong(int(G%jed, kind=c_long))
  args(6) = PyFloat_FromDouble(max_depth)

  method_name = py_string(name)
  ret = PyObject_VectorcallMethod(method_name, args(:), int(size(args), kind=c_size_t), PyObject_None)
  call Py_INCREF(ret)

  err = PyErr_Occurred()
  if (c_associated(err)) then
    call PyErr_Print
    call MOM_error(FATAL, "python interface raised exception")
  end if

  call Py_DECREF(method_name)

  do i = 2, size(args)
    call Py_DECREF(args(i))
  end do
  deallocate(args)

  arr = PyArray_FromAny(ret, PyArray_DescrFromType(NPY_FLOAT64), &
    int(2, kind=c_int), int(2, kind=c_int), NPY_ARRAY_F_CONTIGUOUS, c_null_ptr)
  call Py_INCREF(arr)

  err = PyErr_Occurred()
  if (c_associated(err)) then
    call PyErr_Print
    call MOM_error(FATAL, "exception while converting topo array")
  end if

  arrptr = PyArray_DATA(arr)
  call c_f_pointer(arrptr, Dptr, [G%ied-G%isd+1, G%jed-G%jsd+1])
  D(G%isc:G%iec,G%jsc:G%jec) = Dptr(1:G%iec-G%isc+1,1:G%jec-G%jsc+1)

  call Py_DECREF(arr)
  call Py_DECREF(ret)
end subroutine run_topo_func

subroutine run_velocity_func(velocity_mod, name, G, GV, u, v)
  type(c_ptr), intent(in) :: velocity_mod
  character(len=*), intent(in) :: name
  type(ocean_grid_type), intent(in) :: G
  type(verticalGrid_type), intent(in) :: GV
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(out) :: u
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(out) :: v

  integer :: i, iret
  type(c_ptr) :: ret, err, arr, method_name, arrptr
  type(c_ptr), dimension(:), allocatable :: args
  real, dimension(:,:,:), pointer :: uptr, vptr
  integer, dimension(4) :: bounds

  allocate(args(4))
  args(1) = velocity_mod
  args(2) = PyTuple_New(int(4, kind=c_size_t))
  args(3) = PyTuple_New(int(4, kind=c_size_t))

  bounds = [G%IsdB, G%IedB, G%jsd, G%jed]
  do i = 1, 4
    iret = PyTuple_SetItem(args(2), int(i-1, kind=c_size_t), &
         PyLong_FromLong(int(bounds(i), kind=c_long)))
    if (iret /= 0) call MOM_error(FATAL, "error setting u dim tuple")
  end do

  bounds = [G%isd, G%ied, G%JsdB, G%JedB]
  do i = 1, 4
    iret = PyTuple_SetItem(args(3), int(i-1, kind=c_size_t), &
         PyLong_FromLong(int(bounds(i), kind=c_long)))
    if (iret /= 0) call MOM_error(FATAL, "error setting v dim tuple")
  end do

  args(4) = PyLong_FromLong(int(GV%ke, kind=c_long))

  method_name = py_string(name)
  ret = PyObject_VectorcallMethod(method_name, args(:), int(size(args), kind=c_size_t), PyObject_None)
  call Py_INCREF(ret)

  err = PyErr_Occurred()
  if (c_associated(err)) then
    call PyErr_Print
    call MOM_error(FATAL, "python interface raised exception")
  end if

  call Py_DECREF(method_name)

  do i = 2, size(args)
    call Py_DECREF(args(i))
  end do
  deallocate(args)

  ! extract u array
  arr = PyTuple_GetItem(ret, int(0, kind=c_size_t))
  arr = PyArray_FromAny(arr, PyArray_DescrFromType(NPY_FLOAT64), &
       int(3, kind=c_int), int(3, kind=c_int), NPY_ARRAY_F_CONTIGUOUS, c_null_ptr)
  call Py_INCREF(arr)

  err = PyErr_Occurred()
  if (c_associated(err)) then
    call PyErr_Print
    call MOM_error(FATAL, "exception while converting u array")
  end if

  arrptr = PyArray_DATA(arr)
  call c_f_pointer(arrptr, uptr, [G%IedB-G%IsdB+1, G%jed-G%jsd+1, GV%ke])
  u(:,:,:) = uptr(:,:,:)
  call Py_DECREF(arr)
  uptr => null()

  ! extract v array
  arr = PyTuple_GetItem(ret, int(1, kind=c_size_t))
  arr = PyArray_FromAny(arr, PyArray_DescrFromType(NPY_FLOAT64), &
       int(3, kind=c_int), int(3, kind=c_int), NPY_ARRAY_F_CONTIGUOUS, c_null_ptr)
  call Py_INCREF(arr)

  err = PyErr_Occurred()
  if (c_associated(err)) then
    call PyErr_Print
    call MOM_error(FATAL, "exception while converting v array")
  end if

  arrptr = PyArray_DATA(arr)
  call c_f_pointer(arrptr, vptr, [G%ied-G%isd+1, G%JedB-G%JsdB+1, GV%ke])
  v(:,:,:) = vptr(:,:,:)
  call Py_DECREF(arr)

  vptr => null()

  call Py_DECREF(ret)
end subroutine run_velocity_func

subroutine run_ts_func(ts_mod, name, G, GV, T, S)
  type(c_ptr), intent(in) :: ts_mod
  character(len=*), intent(in) :: name
  type(ocean_grid_type), intent(in) :: G
  type(verticalGrid_type), intent(in) :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T, S

  integer :: i
  type(c_ptr) :: ret, err, arr, method_name, arrptr
  type(c_ptr), dimension(:), allocatable :: args
  real, dimension(:,:,:), pointer :: Tptr, Sptr

  allocate(args(6))
  args(1) = ts_mod
  args(2) = PyLong_FromLong(int(G%isd, kind=c_long))
  args(3) = PyLong_FromLong(int(G%ied, kind=c_long))
  args(4) = PyLong_FromLong(int(G%jsd, kind=c_long))
  args(5) = PyLong_FromLong(int(G%jed, kind=c_long))
  args(6) = PyLong_FromLong(int(GV%ke, kind=c_long))

  method_name = py_string(name)
  ret = PyObject_VectorcallMethod(method_name, args(:), int(size(args), kind=c_size_t), PyObject_None)
  call Py_INCREF(ret)

    err = PyErr_Occurred()
  if (c_associated(err)) then
    call PyErr_Print
    call MOM_error(FATAL, "python interface raised exception")
  end if

  call Py_DECREF(method_name)

  do i = 2, size(args)
    call Py_DECREF(args(i))
  end do
  deallocate(args)

  ! extract T array
  arr = PyTuple_GetItem(ret, int(0, kind=c_size_t))
  arr = PyArray_FromAny(arr, PyArray_DescrFromType(NPY_FLOAT64), &
    int(3, kind=c_int), int(3, kind=c_int), NPY_ARRAY_F_CONTIGUOUS, c_null_ptr)
  call Py_INCREF(arr)

  err = PyErr_Occurred()
  if (c_associated(err)) then
    call PyErr_Print
    call MOM_error(FATAL, "exception while converting T array")
  end if

  arrptr = PyArray_DATA(arr)
  call c_f_pointer(arrptr, Tptr, [G%ied-G%isd+1, G%jed-G%jsd+1, G%ke])
  T(G%isc:G%iec,G%jsc:G%jec,:) = Tptr(1:G%iec-G%isc+1,1:G%jec-G%jsc+1,:)
  call Py_DECREF(arr)
  Tptr => null()

  ! extract S array
  arr = PyTuple_GetItem(ret, int(1, kind=c_size_t))
  arr = PyArray_FromAny(arr, PyArray_DescrFromType(NPY_FLOAT64), &
    int(3, kind=c_int), int(3, kind=c_int), NPY_ARRAY_F_CONTIGUOUS, c_null_ptr)
  call Py_INCREF(arr)

  err = PyErr_Occurred()
  if (c_associated(err)) then
    call PyErr_Print
    call MOM_error(FATAL, "exception while converting S array")
  end if

  arrptr = PyArray_DATA(arr)
  call c_f_pointer(arrptr, Sptr, [G%ied-G%isd+1, G%jed-G%jsd+1, G%ke])
  S(G%isc:G%iec,G%jsc:G%jec,:) = Sptr(1:G%iec-G%isc+1,1:G%jec-G%jsc+1,:)
  call Py_DECREF(arr)
  Sptr => null()

  call Py_DECREF(ret)
end subroutine run_ts_func

end module idealised_python
