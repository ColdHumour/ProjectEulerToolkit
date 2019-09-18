# -*- coding: utf-8 -*-
# distutils: language = c++

"""
c_formula_int64.pxd

Declaration file for _formula.pyx, functions can be only accessed by pyx

@author: Jasper Wu
"""

ctypedef long long int64

cdef int64 c_gcd_int64(int64 a, int64 b)

cdef (int64, int64, int64) c_extended_gcd_int64(int64 a, int64 b)

cdef int64 c_pow_int64(int64 a, int64 b, int64 m)

cdef int64 c_isqrt_int64(int64 n)

cdef short c_is_square_int64(int64 n)

cdef int64 c_inv_mod_int64(int64 n, int64 m)

cdef int64 c_sum_mod_int64(int64 n)

cdef int64 c_sum_power_series_mod_int64(int64 i, int64 n, int64 m)

cdef int64 c_sum_floor_int64(int64 n, int64 xmin, int64 xmax)
