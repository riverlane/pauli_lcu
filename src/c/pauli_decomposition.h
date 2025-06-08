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


#include <complex.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

// Decompositions of matrices into Pauli strings.
//
// This C header file may be used as a stand-alone C library, or can be used through the Python&numpy bindings.
//

// interleave bits with zeros, see https://lemire.me/blog/2018/01/08/how-fast-can-you-bit-interleave-32-bit-integers/
static inline uint64_t interleave_uint32_with_zeros(uint32_t input)  {
    uint64_t word = input;
    word = (word ^ (word << 16)) & 0x0000ffff0000ffff;
    word = (word ^ (word << 8 )) & 0x00ff00ff00ff00ff;
    word = (word ^ (word << 4 )) & 0x0f0f0f0f0f0f0f0f;
    word = (word ^ (word << 2 )) & 0x3333333333333333;
    word = (word ^ (word << 1 )) & 0x5555555555555555;
    return word;
}

static inline void pauli_coefficients(uint32_t dim, double complex *data) {
    // Computes coefficients of Pauli decomposition.
    // O(dim^2 log dim) time complexity, O(1) space overhead
    // Parameters:
    //    dim - dimension of the matrix, 2**num_qubits
    //    data - dim x dim matrix (pointer to double complex)
    // Output:
    //    data - Pauli coefficients are written into data
    
    uint32_t i, j, hf;
    double complex v;
    double complex *a, *b;

    // XOR transform
    // #pragma omp parallel for private(a, b, v, i)
    for(j = 0; j < dim; j++) {
        b = data + j * dim;
        #pragma omp parallel for private(a, v)
        for(i = 0; i < j; i++) {
            // swap
            a = data + i * dim;
            v = b[i^j];
            b[i^j] = a[i^j];
            a[i^j] = v;
        }
    }

    #pragma omp parallel for private(a, b, v, i, hf)
    for(j = 0; j < dim; j++) {
        // Hadamard transform
        for(hf = 1; hf < dim; hf *= 2) {
            b = data + j * dim;
            while(b < data + j * dim + dim) {
                a = b;
                b += hf;
                for(i = 0; i < hf; i++){
                    v = *b;
                    *b++ = *a - v;
                    *a = *a + v;
                    a++;
                }
            }
        }

        double complex *row = data + j * dim;
        for (i = 0; i < dim; i++) {
            // Power of i factors
            switch((uint8_t) (__builtin_popcount(i&j) & 0b11)){
                case 1:
                    row[i] *= -I;
                    break;
                case 2:
                    row[i] *= -1.0;
                    break;
                case 3:
                    row[i] *= I;
                    break;
           }
           // Normalisation
            row[i] /= dim;
        }
    }
}

static inline void pauli_coefficients_lexicographic_order(uint32_t num_qubits, double complex *data) {
    // Computes coefficients of Pauli decomposition, in lexicographic order.
    // Parameters:
    //    num_qubits -- number of qubits.
    //    data - dim x dim matrix (pointer to double complex)
    // Output:
    //    data - Pauli coefficients are written into data in lexicographic order

    uint64_t dim = 1;
    dim <<= num_qubits;

    double complex *tmp = malloc((dim << num_qubits) * sizeof(double complex));
    memcpy(tmp, data, (dim << num_qubits) * sizeof(double complex));

    pauli_coefficients(dim, tmp);

    uint32_t i, j;
    for(i = 0; i < dim; i++) {
        for(j = 0; j < dim; j++) {
            uint64_t id = interleave_uint32_with_zeros(i^j) | (interleave_uint32_with_zeros(i) << 1);
            data[id] = tmp[(j << num_qubits) + i];
        }
    }
    free(tmp);
}

static inline void pauli_string_lexicographic_order(uint32_t id, int num_qubits, char *out) {
    // Computes the Pauli string corresponding to the id'th coefficient in lexicographic order
    // Parameters:
    //     id - integer specifying the required string
    //     num_qubits - integer for number of qubits
    // Output:
    //     out - char string with at least num_qubits+1 characters allocated

    for(num_qubits--; num_qubits >= 0; num_qubits--) {
        switch(id >> (2*num_qubits) & 0b11) {
            case 0: *out++ = 'I'; break;
            case 1: *out++ = 'X'; break;
            case 2: *out++ = 'Y'; break;
            case 3: *out++ = 'Z'; break;
        }
    }
    *out = '\0';
}

static inline void pauli_string_ij(uint32_t i, uint32_t j, uint32_t num_qubits, char *out) {
    // Computes the Pauli string corresponding to the (i,j)'th coefficient
    // Parameters:
    //     i, j - integers specifying the required string
    //     num_qubits - integer for number of qubits
    // Output:
    //     out - char string with at least num_qubits+1 characters allocated

    uint64_t id = interleave_uint32_with_zeros(i^j) | (interleave_uint32_with_zeros(j) << 1);
    pauli_string_lexicographic_order(id, num_qubits, out);
}

static inline void inverse_pauli_decomposition(uint32_t dim, double complex *data) {
    // Computes inverse Pauli decomposition.
    // O(dim^2 log dim) time complexity, O(1) space overhead
    // Parameters:
    //    dim - dimension of the matrix, 2**num_qubits
    //    data - dim x dim matrix (pointer to double complex)
    // Output:
    //    data - Pauli coefficients are written into data

    uint32_t i, j, hf;
    double complex v;
    double complex *a, *b, *c;


    a = data;
    // Power of i factors
    for(j = 0; j < dim; j++) {
        for(i = 0; i < dim; i++) {
            switch((uint8_t) (__builtin_popcount(i&j) & 0b11)){
                case 1:
                    *a *= I;
                    break;
                case 2:
                    *a *= -1.0;
                    break;
                case 3:
                    *a *= -I;
                    break;
           }
           a++;
        }
    }
    // Hadamard transform
    c = data;
    for(j = 0; j < dim; j++) {
        for(hf = 1; hf < dim; hf *= 2) {
            b = c;
            while(b < c + dim) {
                a = b;
                b += hf;
                for(i = 0; i < hf; i++){
                    v = *b;
                    *b++ = *a - v;
                    *a = *a + v;
                    a++;
                }
            }
        }
        c+=dim;
    }
    // XOR transform
    b = data;
    for(j = 0; j < dim; j++) {
        a = data;
        for(i = 0; i < j; i++) {
            // swap
            v = b[i^j];
            b[i^j] = a[i^j];
            a[i^j] = v;
            a += dim;
        }
        b += dim;
    }
}
