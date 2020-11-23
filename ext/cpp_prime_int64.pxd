# -*- coding: utf-8 -*-
# distutils: language = c++

"""
cpp_prime_int64.pxd

Declaration file for cpp_prime_int64.pyx, functions can be only accessed by pyx

@author: Jasper Wu
"""

from ProjectEulerToolkit.ext.cpp_types cimport bool, int64, lvec

cdef struct desc:
    int64 N
    int64 Nrt  # isqrt(N)
    bool flag  # N // Nrt == Nrt

cdef lvec get_primes(int64 N)

cdef desc create_desc_info(int64 N)

cdef lvec get_int_quotients(desc &info)

cdef int64 get_iqs_index(int64 n, desc &info)
