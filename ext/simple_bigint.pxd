# -*- coding: utf-8 -*-
# distutils: language = c++

"""
simple_bigint.pxd

Declaration file for simple_bigint.pyx, functions can be only accessed by pyx

@author: Jasper Wu
"""

from ProjectEulerToolkit.ext.cpp_types cimport bool, int64, sbi

cdef sbi copy_bigint(sbi &a)

cdef bool eq_bigint(sbi &a, sbi &b)

cdef bool lt_bigint(sbi &a, sbi &b)

cdef sbi add_bigint(sbi &a, sbi &b)

cdef sbi rm1_bigint(sbi &a)

cdef sbi mul_bigint(sbi a, sbi b)
