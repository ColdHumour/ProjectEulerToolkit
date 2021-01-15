# -*- coding: utf-8 -*-

"""
c_formula_int64.pyx

Cython extension of functions implementing formulas via fast algorithms.
Function list:
    c_gcd_int64
    c_pow_int64
    c_isqrt_int64
    c_is_square_int64
    c_inv_mod_int64
    c_sum_mod_int64
    c_sum_power_series_mod_int64
    c_sum_floor_int64

@author: Jasper Wu
"""

cdef int64 c_gcd_int64(int64 a, int64 b):
    """return gcd(a, b) , for cimport only"""

    if a < 0:
        a = -a
    if b < 0:
        b = -b

    if a == 0:
        return b
    elif b == 0:
        return a

    while a and b:
        if a > b:
            a %= b
        else:
            b %= a
    if a:
        return a
    else:
        return b

cdef int64 c_pow_int64(int64 a, int64 b, int64 m):
    """return (a ** b) % m, for cimport only"""

    cdef int64 r

    assert b >= 0

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

cdef int64 c_isqrt_int64(int64 n):
    """
    return the nearest integer of sqrt(n), for cimport only
    for algorithm, see: http://www.codecodex.com/wiki/Calculate_an_integer_square_root
    """
    

    cdef:
        int64 res = 0
        int64 one

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

cdef short c_is_square_int64(int64 n):
    """return whether n is a perfect square, for cimport only"""

    cdef:
        int64 s

    if n < 0:
        return 0

    s = c_isqrt_int64(n)
    if s * s == n:
        return 1
    else:
        return 0

cdef int64 c_inv_mod_int64(int64 n, int64 m):
    """return n^(-1) mod m using Extended Euclid Algorithm"""

    cdef:
        int64 m0 = m, x0 = 0, x1 = 1

    n %= m
    if n == 0 or m <= 0:
        return 0
    
    while n != 0:
        x0, x1 = x1, x0 - m // n * x1
        m, n = n, m % n

    if m == 1:
        return x0 % m0
    else:
        return 0

cdef int64 c_sum_mod_int64(int64 n):
    """return n%2 + n%3 + ... + n%(n-1), for cimport only"""

    cdef:
        int64 s = 0, i = 0, imax = c_isqrt_int64(n+1)
        int64 a, b, c

    imax = (imax + 1) >> 1
    for i in range(1, imax):
        a = n % (n//(i+1) + 1)
        b = n % (n//i) if i > 1 else 1
        c = (a-b) // i + 1
        s += b*c + i*(c-1)*c // 2

    for j in range(2, n//(i+1) + 1):
        s += n % j

    return s

cdef int64 c_sum_power_series_mod_int64(int64 i, int64 n, int64 m):
    """sum of x^i mod m from i=1 to n, for i = 0, 1, 2, 3"""

    cdef int64 res, r

    if i == 0:
        return n % m
    elif i == 1:
        if n & 1:
            return ((((n + 1) // 2) % m) * (n % m)) % m
        else:
            return (((n // 2) % m) * ((n + 1) % m)) % m
    elif i == 2:
        r = n % 6
        if r == 0:
            res = (n // 6) % m
            res = (res * ((n + 1) % m)) % m
            res = (res * ((2*n + 1) % m)) % m
        elif r == 1:
            res = n % m
            res = (res * (((n + 1) // 2) % m)) % m
            res = (res * (((2*n + 1) // 3) % m)) % m
        elif r == 2:
            res = (n // 2) % m
            res = (res * (((n + 1) // 3) % m)) % m
            res = (res * ((2*n + 1) % m)) % m
        elif r == 3:
            res = (n // 3) % m
            res = (res * (((n + 1) // 2) % m)) % m
            res = (res * ((2*n + 1) % m)) % m
        elif r == 4:
            res = (n // 2) % m
            res = (res * ((n + 1) % m)) % m
            res = (res * (((2*n + 1) // 3) % m)) % m
        else:
            res = n % m
            res = (res * (((n + 1) // 6) % m)) % m
            res = (res * ((2*n + 1) % m)) % m
        return res
    elif i == 3:
        if n & 1:
            res = ((((n + 1) // 2) % m) * (n % m)) % m
        else:
            res = (((n // 2) % m) * ((n + 1) % m)) % m
        return (res * res) % m
    else:
        return 0

cdef int64 c_sum_floor_int64(int64 n, int64 xmin, int64 xmax):
    """sum up n//x from x = xmin to xmax"""

    cdef:
        int64 nrt = c_isqrt_int64(n), res = 0
        int64 real_xmin, real_xmax, a0, a1, ub, lb

    if xmin > n:
        return 0
    if xmax > n:
        xmax = n

    if xmax <= nrt:
        for x in range(xmin, xmax+1):
            res += n // x
    elif xmin >= nrt:
        real_xmin = n // xmax
        real_xmax = n // xmin
        a0 = 0
        a1 = n // real_xmin
        for x in range(real_xmin, real_xmax+1):
            a0 = a1
            a1 = n // (x+1)
            ub = a0
            if a0 > xmax:
                ub = xmax
            lb = a1
            if a1 < xmin-1:
                lb = xmin-1
            res += (ub - lb) * x
    else:
        real_xmin = n // xmax
        if real_xmin > xmin:
            real_xmin = xmin

        a0 = 0
        a1 = n // real_xmin
        for x in range(real_xmin, nrt+1):
            a0 = a1
            a1 = n // (x+1)

            if x >= xmin:
                res += a0

            if a1 < xmax:
                ub = a0
                if a0 > xmax:
                    ub = xmax
                res += (ub - a1) * x

        if x == n // x:
            res -= x
    return res
