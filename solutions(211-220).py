# -*- coding: utf-8 -*-

"""
solutions(211-220).py

Some interesting solutions for problems 211-220 in Project Euler

@author: Jasper Wu
"""

import math

import ProjectEuler.combinatorics as pec
import ProjectEuler.formulas as pef
import ProjectEuler.generators as peg
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

# %load_ext Cython

# %%cython
# import numpy as np
# cimport numpy as np

# def sum_segments(list segs):
#     cdef unsigned short a, b, s, t
#     cdef long tlens
    
#     if not segs:
#         return 0
    
#     segs.sort()
#     s, t = segs[0]
#     tlens = t - s
#     for a,b in segs[1:]:
#         if a > t:
#             tlens += b - a
#             s, t = a, b
#         elif b > t:
#             tlens += b - t
#             t = b
#     return tlens

# def parse_recs(list reclist):
#     cdef dict idict = {}, pdict = {}, segs = {}
#     cdef set endpoints = set()
#     cdef unsigned short b, c, e, f, i, p
#     cdef long v
        
#     for b, _, e, _, _ in reclist:
#         endpoints.add(b)
#         endpoints.add(b+e)
    
#     for i,p in enumerate(sorted(endpoints)):
#         pdict[p] = i
#         idict[i] = p
    
#     for b, c, e, f, _ in reclist:
#         for i in range(pdict[b], pdict[b+e]):
#             p = idict[i]
#             if p in segs:
#                 segs[p].append((c, c+f))
#             else:
#                 segs[p] = [(c, c+f)]
    
#     v = 0
#     for i in range(len(idict)-1):
#         p = idict[i]
#         v += sum_segments(segs.get(p, [])) * (idict[i+1] - p)
#     return v

def pe212(N=50000):
    """
    Cut all cuboids into sub-cuboids with 1-width according to x-axis.
    Then solve the overlapping rectangles problem on the y-z plane.
    
    Record all y-axis endpoints of rectangles, and then merge the sub- 
    rectangles in each interval with same width.
    
    Cython provides huge savings on time and space. But still ~1.7G ~350s
    Need to be carefully reviewed.
    """
    
    def cube_gen():
        sgen = peg.pe_lagged_fibo_generator(6)
        while 1:
            a, b, c, d, e, f = next(sgen)
            a %= 10000
            b %= 10000
            c %= 10000
            d = 1 + (d % 399)
            e = 1 + (e % 399)
            f = 1 + (f % 399)
            yield a, b, c, d, e, f
    
    cg = cube_gen()
    cubdict = {}
    for _ in range(N):
        a, b, c, d, e, f = next(cg)
        cubdict.setdefault(a, []).append((b, c, e, f, d))
    
    a = v = 0
    while a < 10399:
        if a not in cubdict:
            a += 1
            continue

        for b, c, e, f, d in cubdict[a]:
            if d > 1:
                cubdict.setdefault(a+1, []).append((b, c, e, f, d-1))

        v += parse_recs(cubdict[a])
        a += 1
    
    # answer: 328968937309
    return v