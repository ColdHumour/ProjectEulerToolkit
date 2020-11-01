# -*- coding: utf-8 -*-
# distutils: language = c++

"""
simple_bigint.pyx

Cython extension of functions implementing big integers by c++ containers.
Function list:
    copy_bigint(sbi &a)
    eq_bigint(sbi &a, sbi &b)
    lt_bigint(sbi &a, sbi &b)
    add_bigint(sbi &a, sbi &b)
    rm1_bigint(sbi &a)
    mul_bigint(sbi a, sbi b)

@author: Jasper Wu
"""

from cpp_types cimport bool, int64, sbi

cdef int64 MAXINT = 2**63 - 1


cdef sbi copy_bigint(sbi &a):
    return sbi(a.first, a.second)

cdef bool eq_bigint(sbi &a, sbi &b):
    return a.first == b.first and a.second == b.second

cdef bool lt_bigint(sbi &a, sbi &b):
    if a.first != b.first:
        return a.first < b.first
    else:
        return a.second < b.second

cdef sbi add_bigint(sbi &a, sbi &b):
    cdef:
        sbi c = sbi(0, 0)
        int64 carry

    if a.second <= MAXINT - b.second:
        c.second = a.second + b.second
        carry = 0
    else:
        c.second = a.second - (MAXINT - b.second) - 1
        carry = 1

    c.first = a.first + b.first + carry
    return c

cdef sbi rm1_bigint(sbi &a):
    cdef sbi b = sbi(a.first // 2, a.second // 2)
    if a.first & 1:
        b = add_bigint(b, sbi(0, 2**62))
    return b

cdef sbi mul_bigint(sbi a, sbi b):
    cdef:
        sbi c, zero = sbi(0, 0)

    if not lt_bigint(a, b):
        c = a
        a = b
        b = c

    c = sbi(0, 0)
    while not eq_bigint(a, zero):
        if a.second & 1:
            c = add_bigint(c, b)
        b = add_bigint(b, b)
        a = rm1_bigint(a)
    return c
