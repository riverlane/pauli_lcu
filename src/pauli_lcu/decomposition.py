# Copyright 2024 RIVERLANE LTD
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom
# the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


import numpy as np
from typing import Tuple
import pauli_lcu_module


def pauli_coefficients(matrix: np.ndarray):
    """
    Calculate Pauli coefficients of a matrix. This function overwrites
    the matrix with Pauli coefficients

    Parameters
    ----------
    matrix : ndarray
        2D array, `complex` type, C-contiguous

    Returns
    -------
    None

    """
    assert matrix.data.c_contiguous
    assert np.iscomplexobj(matrix), 'Provide a complex dtype'
    pauli_lcu_module.wrapper_pauli_coefficients(matrix)


def pauli_strings(num_qubits):
    """
    Calculate Pauli strings. The ordering corresponds to Pauli coefficients
    as calculated in pauli_coefficients function

    Parameters
    ----------
    num_qubits: int
        the number of qubits

    Returns
    -------
    ndarray
        dimension of (2^n, 2^n), and string dtype.
    """
    dim = 1 << num_qubits
    strings = np.empty(shape=(dim, dim), dtype=f"<U{num_qubits}")
    pauli_lcu_module.wrapper_pauli_strings(strings, num_qubits)
    return strings


def pauli_coefficients_xz_phase(matrix: np.ndarray):
    assert matrix.data.c_contiguous
    assert np.iscomplexobj(matrix)
    dim = matrix.shape[0]
    num_qubits = int(np.log2(dim))
    x = np.empty(shape=(1 << num_qubits << num_qubits, num_qubits), dtype=np.int8)
    z = np.empty(shape=(1 << num_qubits << num_qubits, num_qubits), dtype=np.int8)
    phase = np.empty(shape=(1 << num_qubits << num_qubits), dtype=np.int8)
    pauli_lcu_module.pauli_coefficients_xz_phase(matrix, x, z, phase, num_qubits)
    return x, z, phase


def pauli_coefficients_lexicographic(matrix: np.ndarray):
    """
    Calculate Pauli coefficients of a matrix. This function overwrites
    the matrix with Pauli coefficients.
    Pauli coefficients will be ordered according to lexicographic order of
    Pauli strings (for example, ['II', 'IX', 'IY', 'IZ', 'XI', ...]).
    This function is more expensive both in time and space
    since it requires reordering which is not done in place.

    Parameters
    ----------
    matrix : ndarray
        2D array, `complex` type, C-contiguous

    Returns
    -------
    None

    """
    assert matrix.data.c_contiguous
    assert np.iscomplexobj(matrix), 'Provide a complex dtype'
    pauli_lcu_module.wrapper_pauli_coefficients_lexicographic_order(matrix)


def inverse_pauli_decomposition(matrix: np.ndarray):
    """
    Given Pauli coefficients stored in 2D-array (matrix),
    restore the original matrix. matrix will be overwritten with new data

    Parameters
    ----------
    matrix : ndarray
        2D array, `complex` type, C-contiguous

    Returns
    -------
    None

    """
    assert matrix.data.c_contiguous
    assert np.iscomplexobj(matrix), 'Provide a complex dtype'
    pauli_lcu_module.wrapper_inverse_pauli_decomposition(matrix)


def pauli_string_ij(ij: Tuple[int, int], num_qubits: int):
    """
    Return Pauli string for a given indices i,j.

    Parameters
    ----------
    ij : tuple[int, int]
    num_qubits : int

    Returns
    -------
    str

    """
    assert len(ij) == 2
    i, j = ij
    assert 0 <= i < 1 << num_qubits
    assert 0 <= j < 1 << num_qubits

    return pauli_lcu_module.wrapper_pauli_string_ij(i, j, num_qubits)
