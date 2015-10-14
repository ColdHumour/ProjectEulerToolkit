# -*- coding: utf-8 -*-

"""
prime.py

Functions using to dealing with primes related problems.
Function list: 
    P10K, p100k, p1m, p10m
    atkin_sieve
    isprime
    iscoprime
    primeDivisorDecomp
    divisorDecomp

@author: Jasper Wu
"""

import os
import random
from math import sqrt
from fractions import gcd

from . formulas import cprod
from . constant import P10K, P10Kset


BASE_DIR = os.path.dirname(__file__)


def p100k():
    file_path = os.path.join(BASE_DIR, 'p100k')
    plist = [int(p.strip()) for p in file(file_path).readlines()]
    return plist

def p1m():
    file_path = os.path.join(BASE_DIR, 'p1m')
    plist = [int(p.strip()) for p in file(file_path).readlines()]
    return plist

def p10m():
    file_path = os.path.join(BASE_DIR, 'p10m')
    plist = [int(p.strip()) for p in file(file_path).readlines()]
    return plist

def atkin_sieve(limit=1000000):
    """
    Sieve of Atkin, which is a much stronger version than sieve of Eratosthenes
    https://en.wikipedia.org/wiki/Sieve_of_Atkin
    """

    plist = [0] * limit

    # n = 3x^2 + y^2 section
    x = 3
    for i in xrange(0, 12*int(sqrt((limit-1)/3)), 24):
        x += i
        y_limit = int(12*sqrt(limit-x)-36)
        n = x + 16
        for j in xrange(-12, y_limit+1, 72):
            n += j
            plist[n] = not plist[n]

        n = x + 4
        for j in xrange( 12, y_limit+1, 72):
            n += j
            plist[n] = not plist[n]

    # n = 4x^2 + y^2 section
    x = 0
    for i in xrange(4, 8*int(sqrt((limit-1)/4))+4, 8):
        x += i
        n = x + 1
        if x % 3:
            for j in xrange(0, 4*int(sqrt(limit-x))-3, 8):
                n += j
                plist[n] = not plist[n]
        else:
            y_limit = 12 * int(sqrt(limit-x)) - 36
            
            n = x + 25
            for j in xrange(-24, y_limit+1, 72):
                n += j
                plist[n] = not plist[n]

            n = x + 1
            for j in xrange( 24, y_limit+1, 72):
                n += j
                plist[n] = not plist[n]

    # n = 3x^2 - y^2 section
    x = 1
    for i in xrange(3, int(sqrt(limit/2))+1, 2):
        x += 4 * i - 4
        n = 3 * x
        if n > limit:
            y = (int(sqrt(n-limit)) >> 2) << 2
            n -= y * y
            s = 4 * y + 4
        else:
            s = 4

        for j in xrange(s, 4*i, 8):
            n -= j
            if n <= limit and n % 12 == 11:
                plist[n] = not plist[n]

    x = 0
    for i in xrange(2, int(sqrt(limit/2))+1, 2):
        x += 4 * i - 4
        n = 3 * x

        if n > limit:
            y = ((int(sqrt(n-limit)) >> 2) << 2) - 1
            n -= y * y
            s = 4 * y + 4
        else:
            n -= 1
            s = 0

        for j in xrange(s, 4*i, 8):
            n -= j
            if n <= limit and n % 12 == 11:
                plist[n] = not plist[n]

    # eliminate squares        
    for n in xrange(5, int(sqrt(limit))+1, 2):
        if plist[n]:
            for k in range(n*n, limit, n*n):
                plist[k] = False

    return [2,3] + filter(plist.__getitem__, xrange(5,limit,2))

def mr_decompose(n):
    exponentOfTwo = 0
    while n % 2 == 0:
        n /= 2
        exponentOfTwo += 1
    return exponentOfTwo, n

def mr_isWitness(possibleWitness, p, exponent, remainder):
    possibleWitness = pow(possibleWitness, remainder, p)
    if possibleWitness == 1 or possibleWitness == p - 1:
        return False
    for _ in xrange(exponent):
        possibleWitness = pow(possibleWitness, 2, p)
        if possibleWitness == p - 1:
            return False
    return True

def aks_expand_x_1(n): 
    c = 1
    for i in xrange(n/2 + 1):
        c *= (n-i) / (i+1)
        yield c

def isprime(p, accuracy=100, how='mr'):
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
        exponent, remainder = mr_decompose(p - 1)

        for _ in xrange(accuracy):
            possibleWitness = random.randint(2, p - 2)
            if mr_isWitness(possibleWitness, p, exponent, remainder):
                return False
        return True
    elif how == 'aks':
        for i in aks_expand_x_1(p):
            if i % p:
                return False
        return True

def iscoprime(a, b):
    return gcd(a, b) == 1

def pollard_rho(n, rand=False):
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
        
def primeDivisorDecomp(n, rand=False):
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
        
        if isprime(n):
            dlist.append(n)
            clist.append(1)
            return zip(dlist, clist)
    
        p = pollard_rho(n, rand)
        c = 0
        while n % p == 0:
            n /= p
            c += 1
        dlist.append(p)
        clist.append(c)

def divisorDecomp(n, rand=False):
    if n == 1:
        return [1]

    primefactors = primeDivisorDecomp(n, rand)
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


import numpy
def primesfrom2to(n):
    """ Input n>=6, Returns a array of primes, 2 <= p < n """
    sieve = numpy.ones(n/3 + (n%6==2), dtype=numpy.bool)
    for i in xrange(1,int(n**0.5)/3+1):
        if sieve[i]:
            k=3*i+1|1
            sieve[       k*k/3     ::2*k] = False
            sieve[k*(k-2*(i&1)+4)/3::2*k] = False
    return numpy.r_[2,3,((3*numpy.nonzero(sieve)[0][1:]+1)|1)]