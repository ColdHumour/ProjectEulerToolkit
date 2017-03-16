# -*- coding: utf-8 -*-

"""
c_formula_int64.pxd

Declaration file for _formula.pyx, functions can be only accessed by pyx

@author: Jasper Wu
"""

cdef unsigned long long c_gcd_uint64(unsigned long long a, unsigned long long b)

cdef unsigned long long c_pow_uint64(unsigned long long a, unsigned long long b, unsigned long long m)

cdef unsigned long long c_isqrt_uint64(unsigned long long n)

cdef short c_is_square_uint64(unsigned long long n)

cdef unsigned long long c_sum_mod_uint64(unsigned long long n)
