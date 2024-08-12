// Copyright 2024 RIVERLANE LTD
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom
// the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
// DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
// OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#define Py_SSIZE_T_CLEAN


#include "pauli_decomposition.h"

#include <Python.h>
#include <numpy/arrayobject.h>

//
// Pauli coefficients
//

static PyObject* wrapper_pauli_coefficients(PyObject* self, PyObject* args){

    PyObject *in_m;

    if (!PyArg_ParseTuple(args, "O", &in_m))
        return NULL;

    if (!PyArray_Check(in_m)) {
        PyErr_SetString(PyExc_TypeError, "Provide an nd array!");
        return NULL;
    }

    PyArrayObject* in = (PyArrayObject*) in_m;

    double complex* pin = (double complex*)PyArray_DATA(in);
    uint32_t dim = PyArray_DIMS(in)[0];
    
    pauli_coefficients(dim, pin);
    
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject* wrapper_inverse_pauli_decomposition(PyObject* self, PyObject* args){

    PyObject *in_m;

    if (!PyArg_ParseTuple(args, "O", &in_m))
        return NULL;

    if (!PyArray_Check(in_m)) {
        PyErr_SetString(PyExc_TypeError, "Provide an nd array!");
        return NULL;
    }

    PyArrayObject* in = (PyArrayObject*) in_m;

    double complex* pin = (double complex*)PyArray_DATA(in);
    uint32_t dim = PyArray_DIMS(in)[0];

    inverse_pauli_decomposition(dim, pin);

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject* wrapper_pauli_strings(PyObject* self, PyObject* args){

    PyObject *in_m;
    uint32_t num_qubits;

    if (!PyArg_ParseTuple(args, "OI", &in_m, &num_qubits))
        return NULL;

   if (!PyArray_Check(in_m)) {
        PyErr_SetString(PyExc_TypeError, "First arg. provide an nd array!");
        return NULL;
    }

    PyArrayObject* in = (PyArrayObject*) in_m;
    char* out_str = (char*) PyArray_BYTES(in);

    uint32_t dim = 1 << num_qubits;
    uint32_t i, j;
    int32_t n;

    uint64_t id = 0;
    char *tmp_ptr = out_str;
    for(i = 0; i < dim; i++){
        for(j = 0; j < dim; j++){
            id = interleave_uint32_with_zeros(i) | (interleave_uint32_with_zeros(j) << 1);
            for(n=num_qubits-1; n >= 0; n--) {
                switch(id >> (2*n) & 0b11) {
                    case 0: *tmp_ptr++ = 'I'; break;
                    case 1: *tmp_ptr++ = 'X'; break;
                    case 3: *tmp_ptr++ = 'Y'; break;
                    case 2: *tmp_ptr++ = 'Z'; break;
                }
                // Numpy string use 4 bytes per character
                *tmp_ptr++ = 0;
                *tmp_ptr++ = 0;
                *tmp_ptr++ = 0;
            }

        }
    }

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject* wrapper_pauli_string_ij(PyObject* self, PyObject* args){

    uint32_t i;
    uint32_t j;
    uint32_t num_qubits;

    if (!PyArg_ParseTuple(args, "III", &i, &j, &num_qubits))
        return NULL;

    char out_str[num_qubits+1];
    pauli_string_ij(i, j, num_qubits, out_str);

    PyObject* py_str = PyUnicode_FromString(out_str);
    if (!py_str) {
        return NULL;
    }
    return py_str;
}


// 
// Pauli coefficients in lexicographic ordering
//

static PyObject* wrapper_pauli_coefficients_lexicographic_order(PyObject* self, PyObject* args){

    PyObject *in_m;

    if (!PyArg_ParseTuple(args, "O", &in_m))
        return NULL;

    if (!PyArray_Check(in_m)) {
        PyErr_SetString(PyExc_TypeError, "Provide an nd array!");
        return NULL;
    }

    PyArrayObject* in = (PyArrayObject*) in_m;


    double complex* pin = (double complex*)PyArray_DATA(in);
    uint32_t dim = PyArray_DIMS(in)[0];
    
    uint32_t num_qubits = 0;
    while(dim >> num_qubits > 0) {
        num_qubits++;
    }
    num_qubits--;

    pauli_coefficients_lexicographic_order(num_qubits, pin);
    
    Py_INCREF(Py_None);
    return Py_None;
}

//
// Pauli coefficients in same format as qiskit
//

static PyObject* pauli_coefficients_xz_phase(PyObject* self, PyObject* args){

    PyObject *in_m;
    PyObject *in_x;
    PyObject *in_z;
    PyObject *in_phase;
    uint32_t num_qubits;

    if (!PyArg_ParseTuple(args, "OOOOI", &in_m, &in_x, &in_z, &in_phase, &num_qubits))
        return NULL;

    if (!PyArray_Check(in_m)) {
        PyErr_SetString(PyExc_TypeError, "First arg. provide an nd array!");
        return NULL;
    }
    if (!PyArray_Check(in_x)) {
        PyErr_SetString(PyExc_TypeError, "Second arg. provide an nd array!");
        return NULL;
    }
    if (!PyArray_Check(in_z)) {
        PyErr_SetString(PyExc_TypeError, "Third arg. provide an nd array!");
        return NULL;
    }
    if (!PyArray_Check(in_phase)) {
        PyErr_SetString(PyExc_TypeError, "Fourth arg. provide an nd array!");
        return NULL;
    }

    const uint32_t dim = 1 << num_qubits;

    PyArrayObject* in = (PyArrayObject*) in_m;
    PyArrayObject* array_x = (PyArrayObject*) in_x;
    PyArrayObject* array_z = (PyArrayObject*) in_z;
    PyArrayObject* array_phase = (PyArrayObject*) in_phase;

    double complex* pin = (double complex*)PyArray_DATA(in);
    uint8_t* xp = (uint8_t*) PyArray_DATA(array_x);
    uint8_t* zp = (uint8_t*) PyArray_DATA(array_z);
    uint8_t* phasep = (uint8_t*) PyArray_DATA(array_phase);

    pauli_coefficients(dim, pin);

    uint32_t i, j, n;
    for(j = 0; j < dim; j++){
        for(i = 0; i < dim; i++){
           *phasep++ = __builtin_popcount(i&j) & 0b11;
           for (n = 0; n < num_qubits; n++) {
               *xp++ = (uint8_t) (j >> n) & 0b1;
               *zp++ = (uint8_t) (i >> n) & 0b1;
           }
        }
    }
    Py_INCREF(Py_None);
    return Py_None;
}


static PyMethodDef PauliMethods[] =
{
     {"wrapper_pauli_coefficients", wrapper_pauli_coefficients, METH_VARARGS, "calculate Pauli coefficients in symplectic order"},
     {"wrapper_inverse_pauli_decomposition", wrapper_inverse_pauli_decomposition, METH_VARARGS, "given Pauli coefficients in symplectic order, calculate original matrix"},
     {"wrapper_pauli_coefficients_lexicographic_order", wrapper_pauli_coefficients_lexicographic_order, METH_VARARGS, "calculate Pauli coefficients in lexicographic order"},
     {"pauli_coefficients_xz_phase", pauli_coefficients_xz_phase, METH_VARARGS, "calculate coefficients and phases"},
     {"wrapper_pauli_strings", wrapper_pauli_strings, METH_VARARGS, "calculate numpy[str] of Pauli strings in symplectic order"},
     {"wrapper_pauli_string_ij", wrapper_pauli_string_ij, METH_VARARGS, "calculate a Pauli stringr"},
     {NULL, NULL, 0, NULL}
 };


static struct PyModuleDef cModPy =
{
    PyModuleDef_HEAD_INIT,
    "pauli_lcu_module", "Pauli decomposition",
    -1,
    PauliMethods
};

PyMODINIT_FUNC
PyInit_pauli_lcu_module(void)
{
    import_array();
    return PyModule_Create(&cModPy);
}
