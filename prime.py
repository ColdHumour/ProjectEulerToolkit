# -*- coding: utf-8 -*-

"""
prime.py

Functions using to dealing with prime-related problems.
Function list: 
    primes_list
    is_prime
    is_coprime
    prime_divisor_decomp
    all_divisors

@author: Jasper Wu
"""


import random
import numpy as np

from ProjectEulerToolkit.formula import cprod, gcd, sqrt

try:
    from gmpy2 import is_prime
except:
    is_prime = _is_prime

try:
    from ProjectEulerToolkit.ext._prime import _c_primes_list
    primes_list = _c_primes_list
except:
    primes_list = _primes_list


# Supplementry Implementation
def _primes_list(n):
    """Input n>=6, Returns a array of primes, 2 <= p < n"""

    sieve = np.ones(n/3 + (n%6 == 2), dtype=np.bool)
    for i in xrange(1, int(sqrt(n))/3+1):
        if sieve[i]:
            k = (3 * i + 1) | 1
            sieve[          k * k / 3         ::2*k] = False
            sieve[k * (k - 2 * (i&1) + 4) / 3 ::2*k] = False
    return np.r_[2, 3, ((3 * np.nonzero(sieve)[0][1:] + 1) | 1)]

def _mr_decompose(n):
    exponentOfTwo = 0
    while n % 2 == 0:
        n /= 2
        exponentOfTwo += 1
    return exponentOfTwo, n

def _mr_isWitness(possibleWitness, p, exponent, remainder):
    possibleWitness = pow(possibleWitness, remainder, p)
    if possibleWitness == 1 or possibleWitness == p - 1:
        return False
    for _ in xrange(exponent):
        possibleWitness = pow(possibleWitness, 2, p)
        if possibleWitness == p - 1:
            return False
    return True

def _aks_expand_x_1(n): 
    c = 1
    for i in xrange(n/2 + 1):
        c *= (n-i) / (i+1)
        yield c

def _is_prime(p, accuracy=100, how='mr'):
    """
    Miller-Rabin primality test
    https://en.wikipedia.org/wiki/Miller-Rabin_primality_test

    AKS primality test
    https://en.wikipedia.org/wiki/AKS_primality_test
    """
    
    if p < 2: return False
    if p == 2 or p == 3: return True
    
    if how == 'mr':
        numTries = 0
        exponent, remainder = _mr_decompose(p - 1)

        for _ in xrange(accuracy):
            possibleWitness = random.randint(2, p - 2)
            if _mr_isWitness(possibleWitness, p, exponent, remainder):
                return False
        return True
    elif how == 'aks':
        for i in _aks_expand_x_1(p):
            if i % p:
                return False
        return True

def is_coprime(a, b):
    return gcd(a, b) == 1

def _pollard_rho(n, rand=False):
    """
    Pollard rho prime factorization algorithm
    https://en.wikipedia.org/wiki/Pollard's_rho_algorithm
    """

    f = lambda x, c: x*x + c
    if not rand:
        x, c = 1, 1
    else:
        x, c = random.randrange(2, 1e6), random.randrange(2, 1e6)
    
    y, d = x, 1
    while d == 1 and d != n:
        x = f(x, c) % n
        y = f(y, c) % n
        y = f(y, c) % n
        d = gcd(y-x, n)
    return d

P10K = primes_list(10000)
P10Kset = set(P10K)
def prime_divisor_decomp(n, rand=False):
    dlist, clist = [], []
    
    # 奇偶性判断
    c = 0
    while n % 2 == 0:
        n /= 2
        c += 1
    if c:
        dlist.append(2)
        clist.append(c)
    
    # 首先用10000以内的小素数试除
    for p in iter(P10K):
        c = 0
        while n % p == 0:
            n /= p
            c += 1
        if c:
            dlist.append(p)
            clist.append(c)
            
        if n == 1:
            return zip(dlist, clist)
        
        if n in P10Kset: # set的in操作复杂度<=O(log(n))
            dlist.append(n)
            clist.append(1)
            return zip(dlist, clist)

    # 然后用Pollard rho方法生成素因子
    while 1:
        if n == 1:
            return zip(dlist, clist)
        
        if is_prime(n):
            dlist.append(n)
            clist.append(1)
            return zip(dlist, clist)
    
        p = _pollard_rho(n, rand)
        c = 0
        while n % p == 0:
            n /= p
            c += 1
        dlist.append(p)
        clist.append(c)

def all_divisors(n, rand=False):
    if n == 1:
        return [1]

    primefactors = prime_divisor_decomp(n, rand)
    d = len(primefactors)
    clist = [0] * d
    output = []
    while 1:
        output.append(cprod([primefactors[i][0]**clist[i] for i in range(d)]))
        k = 0
        while 1:
            clist[k] += 1
            if clist[k] <= primefactors[k][1]:
                break
            clist[k] = 0
            k += 1
            if k >= d:
                return sorted(output)