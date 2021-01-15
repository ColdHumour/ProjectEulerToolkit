# -*- coding: utf-8 -*-
# distutils: language = c++

"""
cpp_formula_int64.pyx

Cython extension of functions implementing formulas via fast algorithms and c++ containers.
Function list:
    cpp_extended_gcd_int64

@author: Jasper Wu
"""

cdef lvec cpp_extended_gcd_int64(int64 a, int64 b):
    """
    return gcd(a, b), x, y that a*x + b*y = gcd(a, b)
    using Extended Euclid Algorithm, where a, b > 0
    """

    cdef:
        int64 x, y, u, v, q, r, m, n
        lvec output

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

    output = lvec()  # don't know why lvec() will raise syntax error
    output.push_back(b)
    output.push_back(x)
    output.push_back(y)
    # output[1] = x
    # output[2] = y
    return output
