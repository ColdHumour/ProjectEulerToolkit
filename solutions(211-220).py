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


def pe211(N=64000000):
    """
    Brute-force via sieve-method. 
    Slow and using too much cache. Need to optimize.
    """
    
    from gmpy2 import is_square
    
    p = int(math.sqrt(N))
    sieve = [0] * (N + 1)

    n = 1
    while n <= p:
        sieve[n * n] += n * n
        
        m, q = n + 1, N / n
        while m <= q:
            sieve[n * m] += n * n + m * m
            m += 1

        n += 1

    c = 0
    for i,n in enumerate(sieve):
        if is_square(n):
            c += i

    # answer: 1922364685
    return c