# -*- coding: utf-8 -*-

"""
c_prime_int64.pyx

Cython extension of functions using to dealing with prime-related problems.
Function list: 
    c_primes_list_int64
    c_mobius_list_int64
    c_factor_sieve_int64

@author: Jasper Wu
"""

import numpy as np
cimport numpy as np

from cpp_types cimport int64
from c_formula_int64 cimport c_isqrt_int64

ctypedef np.ndarray arr


cpdef arr[short, ndim=1] c_primes_list_int64(int64 n):
    """return primes list for primes < n"""
    
    cdef:
        arr[short, ndim=1] sieve = np.ones(n//3 + (n%6==2), dtype=np.int16)
        int64 i, k, imax = c_isqrt_int64(n)
        
    for i in range(1, imax//3+1):
        if sieve[i]:
            k = (3 * i + 1) | 1
            sieve[       k*k//3     ::2*k] = 0
            sieve[k*(k-2*(i&1)+4)//3::2*k] = 0
    return np.r_[2, 3, (3 * np.nonzero(sieve)[0][1:] + 1) | 1]


cpdef arr[int64, ndim=1] c_mobius_list_int64(int64 n):
    """return mobius function mu(k) for 0 <= k <= n"""

    cdef:
        arr[int64, ndim=1] sieve = np.ones(n+1, dtype=np.int64)
        arr[int64, ndim=1] plist
        int64 p, m
        int64 pmax = c_isqrt_int64(n)

    plist = c_primes_list_int64(pmax+1)
    for p in plist:
        sieve[::p] *= -p
        sieve[::p*p] = 0

    for m in range(1, n+1):
        if sieve[m]:
            if abs(sieve[m]) < m:
                sieve[m] *= -1

            if sieve[m] > 0:
                sieve[m] = 1
            else:
                sieve[m] = -1
    return sieve


cpdef arr[int64, ndim=1] c_factor_sieve_int64(int64 n):
    """return factor sieve for 0 <= k <= n"""

    cdef:
        arr[int64, ndim=1] sieve = np.ones(n+1, dtype=np.int64)
        int64 p

    for p in range(2, n):
        if p * p > n:
            break
        if sieve[p] == 1:
            sieve[p*p::p] = p
    return sieve
