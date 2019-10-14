# -*- coding: utf-8 -*-
# distutils: language = c++

"""
cpp_formula_int64.pxd

Declaration file for cpp_formula_int64.pyx, functions can be only accessed by pyx

@author: Jasper Wu
"""

from libcpp.vector cimport vector

cdef vector[long long] cpp_extended_gcd_int64(long long a, long long b)
