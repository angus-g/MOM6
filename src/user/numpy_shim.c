#include <numpy/ndarrayobject.h>

int _NPY_FLOAT64 = NPY_FLOAT64;
int _NPY_ARRAY_F_CONTIGUOUS = NPY_ARRAY_F_CONTIGUOUS;

int array_shim() {
    return PyArray_ImportNumPyAPI();
}

PyArray_Descr *_PyArray_DescrFromType(int typenum) {
    return PyArray_DescrFromType(typenum);
}

PyObject *_PyArray_FromAny(PyObject *op, PyArray_Descr *dtype, int min_depth, int max_depth, int requirements, PyObject *context) {
    return PyArray_FromAny(op, dtype, min_depth, max_depth, requirements, context);
}

void *_PyArray_DATA(PyArrayObject *arr) {
    return PyArray_DATA(arr);
}
