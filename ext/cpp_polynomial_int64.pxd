# -*- coding: utf-8 -*-
# distutils: language = c++

"""
cpp_polynomial_int64.pxd

Cython extension of functions implementing formulas via fast algorithms and c++ containers.
Function list:
    poly_truncate
    poly_add
    poly_mul
    poly_scalar_mul

@author: Jasper Wu
"""

from ProjectEulerToolkit.ext.cpp_types cimport int64, lvec

cdef int64 poly_truncate(lvec &poly)

cdef lvec poly_add(lvec &poly1, lvec &poly2, int64 MOD=*)

cdef lvec poly_mul(lvec &poly1, lvec &poly2, int64 max_length=*, int64 MOD=*)

cdef lvec poly_scalar_mul(lvec &poly1, int64 t, int64 MOD=*)
