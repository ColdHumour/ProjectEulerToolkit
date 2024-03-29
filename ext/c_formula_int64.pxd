# -*- coding: utf-8 -*-

"""
c_formula_int64.pxd

Declaration file for c_formula_int64.pyx, functions can be only accessed by pyx

@author: Jasper Wu
"""

from ProjectEulerToolkit.ext.cpp_types cimport int64

cdef int64 c_gcd_int64(int64 a, int64 b)

cdef int64 c_add_mod_int64(int64 MOD, int64 a0, int64 a1, int64 a2=*, int64 a3=*,
                           int64 a4=*, int64 a5=*, int64 a6=*, int64 a7=*,
                           int64 a8=*, int64 a9=*)

cdef int64 c_mul_mod_int64(int64 MOD, int64 a0, int64 a1, int64 a2=*, int64 a3=*,
                           int64 a4=*, int64 a5=*, int64 a6=*, int64 a7=*,
                           int64 a8=*, int64 a9=*)

cdef int64 c_pow_int64(int64 a, int64 b, int64 m)

cdef int64 c_isqrt_int64(int64 n)

cdef short c_is_square_int64(int64 n)

cdef int64 c_inv_mod_int64(int64 n, int64 m)

cdef int64 c_sum_over_mod_int64(int64 n)

cdef int64 c_sum_power_series_mod_int64(int64 i, int64 n, int64 m)

cdef int64 c_sum_floor_int64(int64 n, int64 xmin, int64 xmax)
