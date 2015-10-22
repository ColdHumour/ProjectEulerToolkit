# -*- coding: utf-8 -*-

"""
_formula.pyx

Cython extension of functions implementing formulas via fast algorithms.
Function list: 
    _c_sum_mod, _cpower_mod

@author: Jasper Wu
"""

from libc.math cimport sqrt, pow


cdef unsigned long long _c_pow(unsigned long long a, 
                               unsigned long long b, 
                               unsigned long long n):
    """
    return (a ** b) % n, for other pyx only
    """
    
    cdef unsigned long long r

    if n == 1:
        return <unsigned long long>pow(a, b)

    r = a % n
    if r == 0 or r == 1:
        return r
    
    while b:
        if b & 1:
            r = (r * a) % n
        b >>= 1
        a = (a * a) % n
    r %= n
    return r

def _c_sum_mod(unsigned long long n):
    """
    return n%2 + n%3 + ... + n%(n-1)
    """
    
    cdef:
        unsigned long long s = 0, i = 0, imax = <unsigned long long>sqrt(n+1)
        unsigned long long a, b, c
    
    imax = (imax + 1) / 2
    for i in range(1, imax):
        a = n % (n/(i+1) + 1)
        b = n % (n/i) if i > 1 else 1
        c = (a-b) / i + 1
        s += b*c + i*(c-1)*c / 2

    for j in range(2, n/(i+1) + 1):
        s += n % j

    return s