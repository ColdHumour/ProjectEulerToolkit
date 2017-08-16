# -*- coding: utf-8 -*-

"""
c_formula_int64.pxd

Declaration file for _formula.pyx, functions can be only accessed by pyx

@author: Jasper Wu
"""

cdef long long c_gcd_int64(long long a, long long b)

cdef long long c_pow_int64(long long a, long long b, long long m)

cdef long long c_isqrt_int64(long long n)

cdef short c_is_square_int64(long long n)

cdef long long c_sum_mod_int64(long long n)

cdef long long c_sum_power_series_mod_int64(long long i, long long n, long long m)
