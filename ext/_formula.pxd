# -*- coding: utf-8 -*-

"""
_formula.pxd

Declaration file for _formula.pyx, functions can be only accessed by pyx

@author: Jasper Wu
"""

cdef unsigned long long _c_pow(unsigned long long a, 
                               unsigned long long b, 
                               unsigned long long n)