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


def _uint32_t_deinterleave_lowuint32(word):
    """
    deinterleave https://github.com/lemire/Code-used-on-Daniel-Lemire-s-blog/blob/master/2018/01/08/interleave.c
    """
    word = word
    word &= 0x5555555555555555
    word = (word ^ (word >> 1 )) & 0x3333333333333333
    word = (word ^ (word >> 2 )) & 0x0f0f0f0f0f0f0f0f
    word = (word ^ (word >> 4 )) & 0x00ff00ff00ff00ff
    word = (word ^ (word >> 8 )) & 0x0000ffff0000ffff
    word = (word ^ (word >> 16)) & 0x00000000ffffffff
    return word


def _interleave_uint32_with_zeros(input):
    """
    // interleave bits with zeros, see https://lemire.me/blog/2018/01/08/how-fast-can-you-bit-interleave-32-bit-integers/
    """

    word = input
    word = (word ^ (word << 16)) & 0x0000ffff0000ffff
    word = (word ^ (word << 8 )) & 0x00ff00ff00ff00ff
    word = (word ^ (word << 4 )) & 0x0f0f0f0f0f0f0f0f
    word = (word ^ (word << 2 )) & 0x3333333333333333
    word = (word ^ (word << 1 )) & 0x5555555555555555
    return word


def interleave(i, j):
    return _interleave_uint32_with_zeros(i) | (_interleave_uint32_with_zeros(j) << 1)


def deinterleave(word: int):
    i = _uint32_t_deinterleave_lowuint32(word)
    j = _uint32_t_deinterleave_lowuint32(word>>1)
    return (i, j)


def lex_indices(id: int):
    """

    Parameters
    ----------
    id : int
        id corresponding to a pauli string in lexicographical order

    Returns
    -------
    tuple[int, int]
        indices to access the same string in pauli_strings
    """
    i, j = deinterleave(id)
    return i ^ j, j
