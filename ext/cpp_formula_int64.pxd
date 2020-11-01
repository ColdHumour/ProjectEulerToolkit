# -*- coding: utf-8 -*-
# distutils: language = c++

"""
cpp_formula_int64.pxd

Declaration file for cpp_formula_int64.pyx, functions can be only accessed by pyx

@author: Jasper Wu
"""

from ProjectEulerToolkit.ext.cpp_types cimport int64, lvec

cdef lvec cpp_extended_gcd_int64(int64 a, int64 b)
