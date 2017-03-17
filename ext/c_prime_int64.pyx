# -*- coding: utf-8 -*-

"""
_prime.pyx

Cython extension of functions using to dealing with prime-related problems.
Function list: 
    c_primes_list_int64
    c_mobius_list_int64
    c_factor_sieve_int64

@author: Jasper Wu
"""

import numpy as np
cimport numpy as cnp

from ProjectEulerToolkit.ext.c_formula_int64 cimport c_isqrt_int64


def c_primes_list_int64(long long n):
    """return primes list for primes < n"""
    
    cdef:
        cnp.ndarray[short, ndim=1] sieve = np.ones(n//3 + (n%6==2), dtype=np.int16)
        long long i, k, imax = c_isqrt_int64(n)
        
    for i in range(1, imax//3+1):
        if sieve[i]:
            k = (3 * i + 1) | 1
            sieve[       k*k//3     ::2*k] = 0
            sieve[k*(k-2*(i&1)+4)//3::2*k] = 0
    return np.r_[2, 3, (3 * np.nonzero(sieve)[0][1:] + 1) | 1]


def c_mobius_list_int64(long long n):
    """return mobius function mu(k) for 0 <= k <= n"""

    cdef:
        cnp.ndarray[long long, ndim=1] sieve = np.ones(n+1, dtype=np.int64)
        cnp.ndarray[long long, ndim=1] plist
        long long p, m
        long long pmax = c_isqrt_int64(n)

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


def c_factor_sieve_int64(long long n):
    """return factor sieve for 0 <= k <= n"""

    cdef:
        cnp.ndarray[long long, ndim=1] sieve = np.ones(n+1, dtype=np.int64)
        long long p

    for p in range(2, n):
        if p * p > n:
            break
        if sieve[p] == 1:
            sieve[p*p::p] = p
    return sieve
