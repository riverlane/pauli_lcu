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

# Run the tests with:
# pytest test/test_pauli_lcu.py

import itertools

import numpy as np
import pauli_lcu
import pytest


def pauli_string_as_matrix(string):
    paulis = dict()
    paulis["I"] = np.array([[1, 0], [0, 1]])
    paulis["X"] = np.array([[0, 1], [1, 0]])
    paulis["Y"] = np.array([[0, -1j], [1j, 0]])
    paulis["Z"] = np.array([[1, 0], [0, -1]])

    matrix = np.array([[1]], dtype=complex)
    for p in string:
        matrix = np.kron(matrix, paulis[p])

    return matrix


@pytest.mark.parametrize("num_qubits", range(1, 6))
def test_decomposition(num_qubits: int):
    dim = 1 << num_qubits
    matrix_orig = np.random.random((dim, dim)) + 1j * np.random.random(
        (dim, dim)
    )
    matrix = np.copy(matrix_orig)

    pauli_lcu.pauli_coefficients(matrix)
    strings = pauli_lcu.pauli_strings(num_qubits)

    for s, coeff in np.nditer([strings, matrix]):
        matrix_orig -= coeff * pauli_string_as_matrix(str(s))

    assert np.allclose(matrix_orig, 0)


@pytest.mark.parametrize("num_qubits", range(1, 6))
def test_decomposition_lexicographic(num_qubits: int):
    dim = 1 << num_qubits
    matrix_orig = np.random.random((dim, dim)) + 1j * np.random.random(
        (dim, dim)
    )
    matrix = np.copy(matrix_orig)

    pauli_lcu.pauli_coefficients_lexicographic(matrix)
    matrix = matrix.flatten()
    strings = itertools.product("IXYZ", repeat=num_qubits)

    for s, coeff in zip(strings, matrix):
        matrix_orig -= coeff * pauli_string_as_matrix(s)

    assert np.allclose(matrix_orig, 0)


@pytest.mark.parametrize("n", range(1, 6))
def test_inverse_pauli_decomposition(n: int):
    dim = 1 << n
    matrix_orig = np.random.rand(dim, dim).astype(complex)
    matrix = matrix_orig.copy()
    # apply pauli decomp and then inverse. Should get the original matrix
    pauli_lcu.pauli_coefficients(matrix)
    pauli_lcu.inverse_pauli_decomposition(matrix)
    assert np.allclose(matrix, matrix_orig)


@pytest.mark.parametrize("n", range(1, 6))
def test_pauli_string_ij(n: int):
    pstrings = pauli_lcu.pauli_strings(n)
    it = np.nditer(pstrings, flags=["multi_index"])
    for x in it:
        assert x == pauli_lcu.pauli_string_ij(it.multi_index, n)


@pytest.mark.parametrize("num_qubits", range(1, 6))
def test_lex_indices(num_qubits: int):
    pstrings = pauli_lcu.pauli_strings(num_qubits)
    lex_strings = itertools.product("IXYZ", repeat=num_qubits)
    for id, lex_string in enumerate(lex_strings):
        assert pstrings[pauli_lcu.lex_indices(id)] == "".join(lex_string)
