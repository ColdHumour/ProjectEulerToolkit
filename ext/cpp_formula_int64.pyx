# -*- coding: utf-8 -*-
# distutils: language = c++

"""
cpp_formula_int64.pyx

Cython extension of functions implementing formulas via fast algorithms and c++ containers.
Function list:
    extended_gcd
    tabulate_fac_mod
    tabulate_fac_inv

@author: Jasper Wu
"""

from ProjectEulerToolkit.ext.c_formula_int64 cimport c_inv_mod_int64 as inv


cdef lvec extended_gcd(int64 a, int64 b):
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

    output = lvec()
    output.push_back(b)
    output.push_back(x)
    output.push_back(y)
    return output


cdef lvec tabulate_fac_mod(int64 n, int64 MOD):
    """return vector of i! % MOD for 0 <= i <= n"""

    cdef:
        int64 m
        lvec fac = lvec(n+1, 1)

    for m in range(2, n+1):
        fac[m] = fac[m-1] * m % MOD
    return fac


cdef lvec tabulate_fac_inv(lvec &fac, int64 MOD):
    """return vector of (i!)^(-1) % MOD for 0 <= i <= n"""

    cdef:
        int64 m
        lvec facinv = lvec(fac.size(), 1)

    for m in range(2, fac.size()):
        facinv[m] = inv(fac[m], MOD)
    return facinv
