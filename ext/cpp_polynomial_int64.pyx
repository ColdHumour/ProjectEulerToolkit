# -*- coding: utf-8 -*-
# distutils: language = c++

"""
cpp_polynomial_int64.pyx

Cython extension of functions implementing formulas via fast algorithms and c++ containers.
Function list:
    poly_truncate
    poly_add
    poly_mul
    poly_scalar_mul

@author: Jasper Wu
"""

cdef int64 poly_truncate(lvec &poly):
    """truncate zero terms until the coefficient of the highest power is nonzero"""

    while poly.size() > 1:
        if poly.back():
            break
        else:
            poly.pop_back()
    return 0


cdef lvec poly_add(lvec &poly1, lvec &poly2, int64 MOD=0):
    cdef:
        int64 l1 = poly1.size(), l2 = poly2.size(), i
        lvec poly3

    if l1 < l2:
        poly3 = lvec(poly2)
        for i in range(l1):
            poly3[i] += poly1[i]
            if MOD and poly3[i] >= MOD:
                poly3[i] -= MOD
    else:
        poly3 = lvec(poly1)
        for i in range(l2):
            poly3[i] += poly2[i]
            if MOD and poly3[i] >= MOD:
                poly3[i] -= MOD
    
    poly_truncate(poly3)
    return poly3


cdef lvec poly_mul(lvec &poly1, lvec &poly2, int64 max_length=0, int64 MOD=0):
    cdef:
        int64 l1 = poly1.size(), l2 = poly2.size(), i1, i2, c1, c2, c3
        lvec poly3

    if max_length:
        max_length = min(l1 + l2, max_length)
    else:
        max_length = l1 + l2

    poly3 = lvec(max_length, 0)
    for i1 in range(l1):
        if i1 == max_length:
            break

        c1 = poly1[i1]
        if c1:
            for i2 in range(l2):
                if i1 + i2 == max_length:
                    break

                c2 = poly2[i2]
                if c2:
                    if MOD:
                        poly3[i1+i2] += c1 * c2 % MOD
                        if poly3[i1+i2] >= MOD:
                            poly3[i1+i2] -= MOD
                    else:
                        poly3[i1+i2] += c1 * c2

    poly_truncate(poly3)
    return poly3


cdef lvec poly_scalar_mul(lvec &poly1, int64 t, int64 MOD=0):
    cdef:
        int64 i
        lvec poly2 = lvec(poly1)

    for i in range(poly2.size()):
        if MOD:
            poly2[i] = poly2[i] * t % MOD
        else:
            poly2[i] *= t

    poly_truncate(poly2)
    return poly2
