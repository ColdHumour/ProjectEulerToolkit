# -*- coding: utf-8 -*-
# distutils: language = c++

"""
cpp_prime_int64.pyx

Cython extension of functions implementing primes related algorithms and c++ containers.
Function list:
    get_primes(int64 N)
    create_desc_info(int64 N)
    get_int_quotients(desc &info)
    get_iqs_index(int64 n, desc &info)

@author: Jasper Wu
"""

from libcpp.vector cimport vector as vec

from cpp_types cimport bool, int64, lvec
from c_formula_int64 cimport c_isqrt_int64 as isqrt

cdef struct desc:
    int64 N
    int64 Nrt  # isqrt(N)
    bool flag  # N // Nrt == Nrt


cdef lvec get_primes(int64 N):
    """linear sieve for primes <= N"""

    cdef:
        int64 n, p
        vec[bool] sieve = vec[bool](N+1, 1)
        lvec primes = lvec()
    
    sieve[0] = sieve[1] = 0
    for n in range(2, N+1):
        if sieve[n]:
            primes.push_back(n)

        for p in primes:
            if n * p > N:
                break

            sieve[n * p] = 0
            if n % p == 0:
                break
    return primes


cdef desc create_desc_info(int64 N):
    cdef desc info = desc(N, isqrt(N), 0)
    info.flag = N // info.Nrt == info.Nrt
    return info


cdef lvec get_int_quotients(desc &info):
    """generate all integer quotient of N, or all possible N // m"""
    
    cdef lvec iqs = lvec()

    for n in range(1, info.Nrt+1):
        iqs.push_back(info.N // n)

    if not info.flag:
        iqs.push_back(info.Nrt)

    for n in range(info.Nrt-1, 0, -1):
        iqs.push_back(n)

    return iqs


cdef inline int64 get_iqs_index(int64 n, desc &info):
    """
    return index of n in array N, N//2, ..., N//m, m, m-1 ..., 1
    where m = int(sqrt(N))
    isqrt_flag is flag of whether m == N // m
    """

    if n > info.Nrt:
        return info.N // n - 1
    elif info.flag:
        return 2 * info.Nrt - n - 1
    else:
        return 2 * info.Nrt - n
