# -*- coding: utf-8 -*-
# distutils: language = c++

"""
cpp_formula_int64.pyx

Cython extension of functions implementing formulas via fast algorithms and c++ containers.
Function list:
    cpp_extended_gcd_int64

@author: Jasper Wu
"""

from libcpp.vector cimport vector

ctypedef vector[long long] lvec

cdef vector[long long] cpp_extended_gcd_int64(long long a, long long b):
    """
    return gcd(a, b), x, y that a*x + b*y = gcd(a, b)
    using Extended Euclid Algorithm, where a, b > 0
    """

    cdef:
        long long x, y, u, v, q, r, m, n
        vector[long long] output

    x = v = 0
    y = u = 1
    while a:
        q = b // a
        r = b - q * a
        m = x - u * q
        n = y - v * q

        b = a
        a = r
        x = u
        y = v
        u = m
        v = n

    output = lvec(3)  # don't know why vector[long long]() will raise syntax error
    output[0] = b
    output[1] = x
    output[2] = y
    return output
