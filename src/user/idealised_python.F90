module idealised_python

use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : FATAL, MOM_error
use MOM_file_parser, only : get_param, param_file_type

implicit none ; private

#include <MOM_memory.h>

public idealised_python_topography

character(len=len("idealised_python")) :: mdl = "idealised_python"
logical :: python_initialised = .false.

interface
   subroutine Py_DECREF(obj) bind(C, name="Py_DecRef")
     use iso_c_binding, only : c_ptr

     type(c_ptr), value, intent(in) :: obj
   end subroutine Py_DECREF
end interface

interface
   function PyErr_Occurred() bind(C, name="PyErr_Occurred")
     use iso_c_binding, only : c_ptr

     type(c_ptr) :: PyErr_Occurred
   end function PyErr_Occurred
end interface

interface
   subroutine PyErr_Print() bind(C, name="PyErr_Print")
   end subroutine PyErr_Print
end interface

interface
   function PyImport_ImportModule(name) bind(C, name="PyImport_ImportModule")
     use iso_c_binding, only : c_char, c_ptr

     character(kind=c_char), intent(in) :: name(*)
     type(c_ptr) :: PyImport_ImportModule
   end function PyImport_ImportModule
end interface

interface
   subroutine Py_INCREF(obj) bind(C, name="Py_IncRef")
     use iso_c_binding, only : c_ptr

     type(c_ptr), value, intent(in) :: obj
   end subroutine Py_INCREF
end interface

interface
   subroutine Py_Initialize() bind(C, name="Py_Initialize")
   end subroutine Py_Initialize
end interface

interface
   function PyList_Insert(list, index, obj) bind(C, name="PyList_Insert")
     use iso_c_binding, only : c_int, c_ptr, c_size_t

     type(c_ptr), value, intent(in) :: list
     integer(kind=c_size_t), value, intent(in) :: index
     type(c_ptr), value, intent(in) :: obj
     integer(kind=c_int) :: PyList_Insert
   end function PyList_Insert
end interface

interface
   function PyObject_VectorcallMethod(obj, name, nargs, kwnames) bind(C, name="PyObject_VectorcallMethod")
     use iso_c_binding, only : c_char, c_ptr, c_size_t

     type(c_ptr), value, intent(in) :: obj
     character(kind=c_char), intent(in) :: name(*)
     integer(kind=c_size_t), value, intent(in) :: nargs
     type(c_ptr), value, intent(in) :: kwnames
     type(c_ptr) :: PyObject_VectorcallMethod
   end function PyObject_VectorcallMethod
end interface

interface
   function PySys_GetObject(name) bind(C, name="PySys_GetObject")
     use iso_c_binding, only : c_char, c_ptr

     character(kind=c_char), intent(in) :: name(*)
     type(c_ptr) :: PySys_GetObject
   end function PySys_GetObject
end interface

interface
   function PyUnicode_FromString(str) bind(C, name="PyUnicode_FromString")
     use iso_c_binding, only : c_char, c_ptr

     character(kind=c_char), intent(in) :: str(*)
     type(c_ptr) :: PyUnicode_FromString
   end function PyUnicode_FromString
end interface

contains

subroutine idealised_python_topography(D, G, param_file, max_depth)
  use iso_c_binding, only : c_ptr

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
  call run_topo_func(topo_mod, topo_func, G%isd, G%ied, G%jsd, G%jed, max_depth)
end subroutine idealised_python_topography

subroutine python_init()
  use iso_c_binding, only : c_associated, c_char, c_int, c_null_char, c_ptr, c_size_t, c_intptr_t

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
    python_initialised = .true.
  end if
end subroutine python_init

subroutine load_module(name, mod)
  use iso_c_binding, only : c_associated, c_char, c_null_char, c_ptr

  character(len=*), intent(in) :: name
  type(c_ptr), intent(out) :: mod
  character(kind=c_char) :: cname(len_trim(name) + 1)
  integer :: i, lv
  type(c_ptr) :: err

  lv = len_trim(name)
  do concurrent (i=1:lv)
    cname(i) = name(i:i)
  end do
  cname(lv+1) = c_null_char

  mod = PyImport_ImportModule(cname)
  err = PyErr_Occurred()

  if (c_associated(err)) then
    call PyErr_Print
    call MOM_error(FATAL, "python interface raised exception")
  end if
end subroutine load_module

subroutine run_topo_func(mod, name, isd, ied, jsd, jed, max_depth)
  use iso_c_binding, only : c_char, c_null_char, c_ptr

  type(c_ptr), intent(in) :: mod
  character(len=*), intent(in) :: name
  integer, intent(in) :: isd, ied, jsd, jed
  real, intent(in) :: max_depth
  character(kind=c_char) :: cname(len_trim(name) + 1)
  integer :: i, lv

  type(c_ptr) :: ret, err

  lv = len_trim(name)
  do concurrent (i=1:lv)
    cname(i) = name(i:i)
  end do
  cname(lv+1) = c_null_char

  ! needs argument array and kwargs...
  ret = PyObject_VectorcallMethod(mod, cname, int(0, kind=c_size_t), 
end subroutine run_topo_func
end module idealised_python
