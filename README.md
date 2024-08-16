# Pauli LCU

Pauli LCU is a package for calculating the Pauli Decomposition of a complex matrix $`A=\{a_{r,s}\}_{r,s=0}^{2^{n}-1}`$:

```math
    A=\sum_{r,s=0}^{2^n-1} \alpha_{r,s} \bigotimes_{j=0}^{n-1} P(r_j, s_j)
```

where

```math
    \alpha_{r,s}\coloneqq \frac{i^{-|r\wedge s|}}{2^n} \sum_{q=0}^{2^n-1}a_{q\oplus r,q}(H^{\otimes n})_{q,s}, \qquad P(r_j, s_j) = i^{r_j\wedge s_j} X^{r_j} Z^{s_j} \in \{I, X, Y, Z\},
```

and $r_j, s_j$ are the binary expansion coefficients of $r$ and $s$,  according
to the algorithm presented in "Timothy N. Georges, Bjorn K. Berntson, Christoph Sünderhauf, and Aleksei V. Ivanov, _Pauli Decomposition via the Fast Walsh-Hadamard Transform_ (2024)." 


## Installation
To install from PyPi:

    $ pip install pauli_lcu

To install Python bindings from source:

    $ git clone git@github.com:riverlane/pauli_lcu.git
    $ cd pauli_lcu
    $ pip install .

For use in C, `pauli_decomposition.h` can be used as a stand-alone library. 
Windows is not currently supported. In this case, we recommend to use WSL or Linux VM.

##  Examples for Python

Compute Pauli coefficients of a matrix:
```pycon
>>> import numpy as np
>>> from pauli_lcu import pauli_coefficients
>>> num_qubits = 1
>>> dim = 2 ** num_qubits
>>> matrix = np.arange(dim*dim).reshape(dim, dim).astype(complex)  # generate your matrix
>>> pauli_coefficients(matrix)   # calculate coefficients, stored in matrix, original matrix is overwritten
>>> matrix
array([[ 1.5+0.j , -1.5+0.j ],
       [ 1.5+0.j ,  0. -0.5j]])
```

One can also obtain the corresponding Pauli strings:
```pycon
>>> from pauli_lcu import pauli_strings
>>> strings = pauli_strings(num_qubits)   # calculate Pauli strings
>>> for p_string, p_coeff in np.nditer([strings, matrix]):
...     print(p_string, p_coeff)
...
I (1.5+0j)
Z (-1.5+0j)
X (1.5+0j)
Y -0.5j
```

For large matrices generating all Pauli strings is memory expensive. In this case, one can access a Pauli string as follows:

```pycon
>>> from pauli_lcu import pauli_string_ij
>>> it = np.nditer(matrix, flags=['multi_index'])
>>> for x in it:
...    print(pauli_string_ij(it.multi_index, num_qubits), x)
I (1.5+0j)
Z (-1.5+0j)
X (1.5+0j)
Y -0.5j
```
If you would like to access the array and Pauli strings in lexicographical order (eg $II$, $IX$, ...) then you can do:

```pycon
>>> from pauli_lcu import lex_indices
>>> for id in range(dim**2):
...     mult_ind = lex_indices(id)
...     print(pauli_string_ij(mult_ind, num_qubits), matrix[mult_ind])
I (1.5+0j)
X (1.5+0j)
Y -0.5j
Z (-1.5+0j)
```
Alternatively, you can use `pauli_lcu.pauli_coefficients_lexicographic()` which uses extra memory to return lexicographically ordered coefficients.


If you would like to restore original matrix elements apply inverse Pauli decomposition:
```pycon
>>> from pauli_lcu.decomposition import inverse_pauli_decomposition
>>> inverse_pauli_decomposition(matrix)
>>> matrix
array([[0.+0.j, 1.+0.j],
       [2.+0.j, 3.+0.j]])
```

We also added ZX-decomposition which can be used in Qiskit (you need it to install qiskit for this example):

```pycon
>>> from pauli_lcu.decomposition import pauli_coefficients_xz_phase
>>> pauli_coefficients_xz_phase(matrix)
>>> x, z, phase = pauli_coefficients_xz_phase(matrix)
>>> x 
array([[0],
       [1],
       [0],
       [1]], dtype=int8)
>>> from qiskit.quantum_info import PauliList
>>> PauliList.from_symplectic(z, x)
PauliList(['I', 'Z', 'X', 'Y'])
```

## References

If you find this package useful please cite our paper:

Timothy N. Georges, Bjorn K. Berntson, Christoph Sünderhauf, and Aleksei V. Ivanov, _Pauli Decomposition via the Fast Walsh-Hadamard Transform_ (2024).
