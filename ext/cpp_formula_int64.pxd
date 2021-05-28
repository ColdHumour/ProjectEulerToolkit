# -*- coding: utf-8 -*-
# distutils: language = c++

"""
cpp_formula_int64.pxd

Declaration file for cpp_formula_int64.pyx, functions can be only accessed by pyx

@author: Jasper Wu
"""

from ProjectEulerToolkit.ext.cpp_types cimport int64, lvec

cdef lvec extended_gcd(int64 a, int64 b)

cdef lvec tabulate_fac_mod(int64 n, int64 MOD)

cdef lvec tabulate_fac_inv(lvec &fac, int64 MOD)
