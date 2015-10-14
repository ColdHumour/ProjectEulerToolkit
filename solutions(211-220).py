# -*- coding: utf-8 -*-

"""
solutions(211-220).py

Some interesting solutions for problems 211-220 in Project Euler

@author: Jasper Wu
"""

import math

import ProjectEuler.combinatorics as pec
import ProjectEuler.formulas as pef
import ProjectEuler.prime as pep


# %load_ext Cython

# %%cython
# import numpy as np
# cimport numpy as np
# from libc.math cimport sqrt

# def is_square(unsigned long long n, unsigned long long v):
#     cdef unsigned long long s, t
#     t = n * n + v
#     s = int(sqrt(t))
#     return s * s == t

# def c_pe211(unsigned long long N):
#     cdef np.ndarray[unsigned long long, ndim=1] sieve
#     cdef unsigned long long m, n = 1, c = 1
#     cdef int p, q

#     sieve = np.ones(N+1, dtype=np.uint64)
#     p = int(sqrt(N))
    
#     for n in xrange(2, p+1):
#         sieve[n * n] += n * n
#         m, q = n + 1, N / n
#         while m <= q:
#             sieve[n * m] += n * n + m * m
#             m += 1
#         n += 1
        
#     for n in xrange(2, N):
#         if is_square(n, sieve[n]):
#             c += n
#     return c

def pe211(N=64000000):
    """
    Brute-force. Using Cython to accelerate and reduce memory.
    Wrap to record time. Final performance: ~16s, ~600M
    """

    # answer: 1922364685
    return c_pe211(N)