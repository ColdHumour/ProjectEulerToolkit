# -*- coding: utf-8 -*-
# distutils: language = c++

"""
cpp_prime_int64.pyx

Cython extension of functions implementing primes related algorithms and c++ containers.
Function list:
    get_primes(int64 N)
    get_factor_sieve(int64 N)

    get_mobius_vec(int64 N)
    get_mertens_vec(int64 N)
    get_mertens(int64 n, int64 MOD, int64 L, lvec &Mvec, llmap &Mcache)

    create_desc_info(int64 N)
    get_int_quotients(desc &info)
    get_iqs_index(int64 n, desc &info)

@author: Jasper Wu
"""

from libcpp.vector cimport vector as vec

from ProjectEulerToolkit.ext.c_formula_int64 cimport c_isqrt_int64 as isqrt


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


cdef lvec get_factor_sieve(int64 N):
    """factor sieve for 0 <= k <= N"""

    cdef:
        lvec sieve = lvec(N+1, 1)
        int64 p, m

    for p in range(2, N):
        if p * p > N:
            break
        if sieve[p] == 1:
            m = p * p
            while m <= N:
                sieve[m] = p
                m += p
    return sieve


cdef lvec get_mobius_vec(int64 N):
    """return mobius function mu(k) for 0 <= k <= N"""

    cdef:
        lvec sieve = lvec(N+1, 1)
        lvec plist = get_primes(isqrt(N))
        int64 p, p2, m

    for p in plist:
        m = p
        p2 = p * p
        while m <= N:
            if sieve[m]:
                if m % p2:
                    sieve[m] *= -p
                else:
                    sieve[m] = 0
            m += p

    sieve[0] = 0
    for m in range(1, N+1):
        if sieve[m]:
            if abs(sieve[m]) < m:
                sieve[m] = -1 if sieve[m] > 0 else 1
            else:
                sieve[m] = 1 if sieve[m] > 0 else -1
    return sieve


cdef lvec get_mertens_vec(int64 N):
    """
    return mertens function M(k) for 0 <= k <= N
    M(k) = sum mu(i) for 1 <= i <= k
    """

    cdef:
        lvec sieve = lvec(N+1, 1)
        lvec plist = get_primes(isqrt(N))
        int64 p, p2, m

    for p in plist:
        m = p
        p2 = p * p
        while m <= N:
            if sieve[m]:
                if m % p2:
                    sieve[m] *= -p
                else:
                    sieve[m] = 0
            m += p

    sieve[0] = 0
    for m in range(1, N+1):
        if sieve[m]:
            if abs(sieve[m]) < m:
                p = -1 if sieve[m] > 0 else 1
            else:
                p = 1 if sieve[m] > 0 else -1
        else:
            p = 0
        sieve[m] = sieve[m-1] + p
    return sieve


cdef int64 get_mertens(int64 n, int64 MOD, int64 L, lvec &Mvec, llmap &Mcache):
    """
    return mertens function M(n)
    some parameters must be defined outside
        int64 L = <int64>((<double>N/log(log(N)))**(2/3.))
        lvec Mvec = get_mertens_vec(L)
        llmap Mcache = llmap()
    """

    cdef:
        int64 nrt, b0, b1, d, r, res

    if n <= L:
        return Mvec[n]

    if Mcache.find(n) == Mcache.end():
        nrt = isqrt(n)
        res = (1 + n//2 - n) % MOD
        b1 = n // 2
        for d in range(2, nrt+1):
            b0 = b1
            b1 = n // (d+1)

            r = get_mertens(n//d, MOD, L, Mvec, Mcache)
            res -= r
            if res < 0:
                res += MOD

            r = get_mertens(d, MOD, L, Mvec, Mcache)
            res -= r * (b0 - b1) % MOD
            if res < 0:
                res += MOD

        if nrt == n // nrt:
            res += r
            if res > MOD:
                res -= MOD

        Mcache[n] = res
    return Mcache[n]



# ---------------------- tool functions for recursions like pi(x) ------------------------

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
