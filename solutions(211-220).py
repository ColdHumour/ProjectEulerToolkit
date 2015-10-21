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

def pe213(N=30, B=50):
    """
    Compute out the B-step distribution of each flea.
    Then compute the vacant probability like birthday problem.
    """
    
    NN = N * N

    NEXTMOVE = {}
    for i in range(NN):
        if i < N:
            move = [i+N]
        elif i > NN - N - 1:
            move = [i-N]
        else:
            move = [i-N, i+N]

        if i % N == 0:
            move += [i+1]
        elif i % N == N - 1:    
            move += [i-1]
        else:
            move += [i-1, i+1]

        NEXTMOVE[i] = dict.fromkeys(move, 1./len(move))

    def jumpdist(n, b):
        dist = {n: 1}
        for _ in range(b):
            temp = {}
            for i,p in dist.items():
                for j,q in NEXTMOVE[i].items():
                    if j in temp:
                        temp[j] += p * q
                    else:
                        temp[j] = p * q
            dist = temp
        return dist

    grids = dict.fromkeys(range(NN), 1)
    for i in range(NN):
        jd = jumpdist(i, B)
        for i,p in jd.items():
            grids[i] *= 1 - p

    # answer: 330.721154
    return round(sum(grids.values()), 6)

# %load_ext Cython

# %%cython
# import numpy as np
# cimport numpy as np
# from libc.math cimport sqrt

# def primes_st_n(unsigned long long n):
#     """Return primes list for primes < n."""
    
#     cdef:
#         np.ndarray[short, ndim=1] sieve = np.ones(n/3 + (n%6==2), 
#                                               dtype=np.int16)
#         unsigned long long i, k
#         unsigned long long i2max = <unsigned long long>sqrt(n)
        
#     for i in range(1, i2max/3+1):
#         if sieve[i]:
#             k = (3 * i + 1) | 1
#             sieve[       k*k/3     ::2*k] = 0
#             sieve[k*(k-2*(i&1)+4)/3::2*k] = 0
#     return np.r_[2, 3, (3 * np.nonzero(sieve)[0][1:] + 1) | 1]

# def odd_phis_se_n(unsigned long long n):
#     """
#     Return Euler's totient funtion phi(x) for every odd number <= n.
    
#     Sieve method, which based on an array of length int((n-3)/2), with map 
#     i -> 2*i+3. Using facts:
    
#         (1) phi(p) = p - 1
#         (2) phi(p^k * q) = (p^k - p^(k-1)) * phi(q), where gcd(p, q) = 1
    
#     For even numbers, it's very fast to get using (2) when p = 2.
#     """
    
#     cdef:
#         unsigned long long ni = n / 2 - 1 + (n & 1)
#         np.ndarray[unsigned long long, ndim=1] sieve = np.zeros(ni, 
#                                                            dtype=np.uint64)
#         unsigned long long i, j, m, f, q
    
#     i = 0
#     while i < ni:
#         m = 2 * i + 3
#         if sieve[i] == 0:
#             sieve[i] = m - 1
#             j = 0
#             while (2*j+3) * m <= n:
#                 if sieve[j]:
#                     q = 2 * j + 3
#                     f = m - 1
#                     while q % m == 0:
#                         f *= m
#                         q //= m
                        
#                     if q == 1:
#                         sieve[2*i*j+3*(i+j)+3] = f
#                     else:
#                         sieve[2*i*j+3*(i+j)+3] = f * sieve[q/2-1]
#                 j += 1
#         i += 1
#     return sieve
        
# def c_pe214(long long n, short chain):
#     cdef:
#         np.ndarray[unsigned long long, ndim=1] phis = odd_phis_se_n(n)
#         np.ndarray[long long, ndim=1] primes = primes_st_n(n)
#         unsigned long long s, p
#         unsigned long i, f, p_amount = len(primes)
#         short j
    
#     for i in range(2, p_amount):
#         p = primes[i] - 1
#         for j in range(chain):
#             if p == 2:
#                 break
            
#             f = 1
#             while p & 1 == 0:
#                 f <<= 1
#                 p >>= 1

#             if f > 1:
#                 f >>= 1

#             if p == 1:
#                 p = f
#             else:
#                 p = f * phis[p/2-1]

#         if j == chain - 3:
#             s += primes[i]
#     return s

def pe214(n=40000000, c=25):
    """
    Construct primes list and odd phis list.
    Then check each prime with its phi chain.
    Sum those who satisfying chain length condition.
    """

    return c_pe214(n, c)

def pe215(width=32, height=10):
    """
    avail := {one-layer: [possible crack-free ones-layers]}
    state := {current layer: possible ways to achieve current layer}
    Put layer over layer and refresh state.
    """
    
    def one_layer_combinations(i, n):
        if n - i == 2 or n - i == 3:
            yield [i]
        elif n - i == 4:
            yield [i, i+2]
        else:
            for m in (i+2, i+3):
                for res in one_layer_combinations(m, n):
                    yield [i] + res if i else res

    avail = {}
    for b in one_layer_combinations(0, width):
        avail[tuple(b)] = []

    lines = avail.keys()
    for i,b1 in enumerate(lines):
        for b2 in lines[i:]:
            if not set(b1).intersection(b2):
                avail[b1].append(b2)
                avail[b2].append(b1)

    state = dict.fromkeys(lines, 1)
    for _ in range(height-1):
        new_stat = {}
        for b1,n in state.iteritems():
            for b2 in avail[b1]:
                if b2 in new_stat:
                    new_stat[b2] += n
                else:
                    new_stat[b2] = n
        state = new_stat
    
    # answer: 806844323190414
    return sum(state.values())

# %load_ext Cython

# %%cython
# import numpy as np
# cimport cython
# cimport numpy as np
# from libc.math cimport sqrt

# def isqrt(long n):
#     cdef long sn = <long>sqrt(n)
#     if sn * sn == n:
#         return sn
#     else:
#         return 0
    
# def sifted_primes_st_n(long n):
#     """Return primes list for primes < n."""
    
#     cdef:
#         np.ndarray[short, ndim=1] sieve
#         np.ndarray[long long, ndim=1] output
#         long i, k, i2max = <long>sqrt(n)
    
#     sieve = np.ones(n/3 + (n%6==2), dtype=np.int16)
#     for i in range(1, i2max/3+1):
#         if sieve[i]:
#             k = (3 * i + 1) | 1
#             sieve[       k*k/3     ::2*k] = 0
#             sieve[k*(k-2*(i&1)+4)/3::2*k] = 0
#     output = np.r_[(3 * np.nonzero(sieve)[0][1:] + 1) | 1]
#     return output[(output % 8 == 1) | (output % 8 == 7)]

# def tonelli_shanks(long n, long p):
#     """
#     Solve the equation
#         x^2 = n (mod p)
#     Details see: http://en.wikipedia.org/wiki/Tonelli-Shanks_algorithm
#     """
    
#     cdef:
#         long long p2 = (p-1)/2, r = isqrt(n)
#         long long z, c, i, t, t2
#         long q = p-1, s = 0
    
#     if r:
#         pass
#     elif p % 4 == 3:
#         r = <long long>pow(n, (p2+1)/2, p)
#     else:
#         for z in xrange(3, p2+1):
#             if <long long>pow(z, p2, p) == p - 1: # Euler's criterion
#                 break

#         while q & 1 == 0:
#             q >>= 1
#             s += 1

#         c = <long long>pow(z, q, p)
#         r = <long long>pow(n, (q+1)/2, p)
#         t = <long long>pow(n, q, p)
#         while t > 1:
#             t2 = t
#             for i in xrange(1, s):
#                 t2 = (t2 * t2) % p
#                 if t2 == 1:
#                     break

#             b = pow(c, 1 << (s-i-1), p)
#             r = (r * b) % p
#             t = (t * b * b) % p
#             c = (b * b) % p
#             s = i
        
#     if r > p2:
#         r = p - r
#     return r

# def c_pe216(long long N):
#     cdef:
#         np.ndarray[short, ndim=1] sieve
#         np.ndarray[long long, ndim=1] sps
#         unsigned long i, d, p
#         long long b
    
#     sieve = np.ones(N+1, dtype=np.int16)
#     sps = sifted_primes_st_n(<long>sqrt(2*N*N+1))
#     d = len(sps)
#     for i in range(d):
#         p = sps[i]
#         b = tonelli_shanks((p + 1) / 2, p)
#         sieve[p-b::p] = 0
#         sieve[p+b::p] = 0
#     return sieve.sum() - 2

def pe216(N=50000000):
    """
    (2n^2 - 1) has factor p when n = k*p + b, where: 
        (1) p is odd prime and
        (2) 2b^2 = 1 (mod p)
    Using Second Supplementary Law of quadratic residue, we can find that
    2b^2 = 1 (mod p) has solution only when p = 1, 7 (mod 8)
    Then applying Tonelli-Shanks algorithm to work out b.
    The remaining work is just simple sieve-method.
    """
    
    # answer: 5437849
    return c_pe216(N)