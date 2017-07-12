# -*- coding: utf-8 -*-

"""
c_formula_uint64.pyx

Cython extension of functions implementing formulas via fast algorithms.
Function list:
    c_gcd_uint64
    c_pow_uint64
    c_isqrt_uint64
    c_is_square_uint64
    c_sum_mod_uint64

@author: Jasper Wu
"""

cdef unsigned long long c_gcd_uint64(unsigned long long a, unsigned long long b):
    """return gcd(a, b) , for cimport only"""

    while a and b:
        if a > b:
            a %= b
        else:
            b %= a
    if a:
        return a
    else:
        return b

cdef unsigned long long c_pow_uint64(unsigned long long a, unsigned long long b, unsigned long long m):
    """return (a ** b) % m, for cimport only"""

    cdef unsigned long long r

    if m == 1:
        return a**b

    a %= m
    if a == 0 or a == 1:
        return a

    if b == 0:
        return 1
    elif b == 1:
        return a
    elif b == 2:
        return (a * a) % m
    elif b == 3:
        return (((a * a) % m) * a) % m
    else:
        r = 1
        while b:
            if b & 1:
                r = (r * a) % m
            b >>= 1
            a = (a * a) % m
        r %= m
        return r

cdef unsigned long long c_isqrt_uint64(unsigned long long n):
    """
    return the nearest integer of sqrt(n), for cimport only
    for algorithm, see: http://www.codecodex.com/wiki/Calculate_an_integer_square_root
    """

    cdef:
        unsigned long long res = 0
        unsigned long long one

    if n >= (1 << 32):
        one = 1 << 62
    elif n >= (1 << 16):
        one = 1 << 30
    else:
        one = 1 << 14

    while one > n:
        one >>= 2

    while one:
        if n >= res + one:
            n -= res + one
            res += 2 * one
        res >>= 1
        one >>= 2
    return res

cdef short c_is_square_uint64(unsigned long long n):
    """return whether n is a perfect square, for cimport only"""

    cdef:
        unsigned long long s

    if n < 0:
        return 0

    s = c_isqrt_uint64(n)
    if s * s == n:
        return 1
    else:
        return 0

cdef unsigned long long c_sum_mod_uint64(unsigned long long n):
    """return n%2 + n%3 + ... + n%(n-1), for cimport only"""

    cdef:
        unsigned long long s = 0, i = 0, imax = c_isqrt_uint64(n+1)
        unsigned long long a, b, c

    imax = (imax + 1) >> 1
    for i in range(1, imax):
        a = n % (n//(i+1) + 1)
        b = n % (n//i) if i > 1 else 1
        c = (a-b) // i + 1
        s += b*c + i*(c-1)*c // 2

    for j in range(2, n//(i+1) + 1):
        s += n % j

    return s
