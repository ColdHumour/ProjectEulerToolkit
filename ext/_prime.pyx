# -*- coding: utf-8 -*-

"""
_prime.pyx

Cython extension of functions using to dealing with prime-related problems.
Function list: 
    _c_primes_list

@author: Jasper Wu
"""

import numpy as np
cimport numpy as cnp

from libc.math cimport sqrt
from ProjectEulerToolkit.ext._formula cimport _c_pow


def _c_primes_list(unsigned long long n):
    """return primes list for primes < n"""
    
    cdef:
        cnp.ndarray[short, ndim=1] sieve = np.ones(n/3 + (n%6==2), 
                                                   dtype=np.int16)
        unsigned long long i, k
        unsigned long long imax = <unsigned long long>sqrt(<double>n)
        
    for i in range(1, imax/3+1):
        if sieve[i]:
            k = (3 * i + 1) | 1
            sieve[       k*k/3     ::2*k] = 0
            sieve[k*(k-2*(i&1)+4)/3::2*k] = 0
    return np.r_[2, 3, (3 * np.nonzero(sieve)[0][1:] + 1) | 1]


def _c_mobius_list(long long n):
    """return mobius function mu(k) for 0 <= k < n"""

    cdef:
        cnp.ndarray[long long, ndim=1] sieve = np.ones(n+1, dtype=np.int64)
        cnp.ndarray[long long, ndim=1] plist
        long long p, m
        long long pmax = <long long>sqrt(<double>n+1)

    plist = _c_primes_list(pmax)
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
