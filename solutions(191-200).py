# -*- coding: utf-8 -*-

"""
solutions(191-200).py

Some interesting solutions for problems 191-200 in Project Euler

@author: Jasper Wu
"""

import math

import ProjectEuler.combinatorics as pec
import ProjectEuler.formulas as pef
import ProjectEuler.prime as pep


def pe191(N=30):
    """
    Keep numbers of state (k, t) when string length is d.
    (k, t) := string containing only 'O' and 'A', 
              with k single 'A's and end with t 'A's.
    
    Map relation:
        d: (k, 0) -> d+1: (k, 0), (k+1, 1) 
        d: (k, 1) -> d+1: (k, 0), (k-1, 2)
        d: (k, 2) -> d+1: (k, 0)
    
    When string length is N-4, replace single 'A' to 'AALAA'.
    When string length is N-3, replace single 'A' to 'ALAA' or 'AALA'.
    When string length is N-1, insert 'L' to N intervals.
    """
    
    states = {(0, 0): 1, (1, 1): 1}
    late = 0
    for d in range(2, N+1):
        new_stat = {}
        for (k, t), n in states.iteritems():
            if t == 0:
                for s in [(k+1, 1), (k, 0)]: 
                    if s in new_stat:
                        new_stat[s] += n
                    else:
                        new_stat[s] = n
            elif t == 1:
                for s in [(k-1, 2), (k, 0)]: 
                    if s in new_stat:
                        new_stat[s] += n
                    else:
                        new_stat[s] = n
            else:
                s = (k, 0) 
                if s in new_stat:
                    new_stat[s] += n
                else:
                    new_stat[s] = n
        states = new_stat

        if d == N - 4:
            late += sum(k * n for (k, _), n in states.iteritems())
        if d == N - 3:
            late += sum(2 * k * n for (k, _), n in states.iteritems())
        if d == N - 1:
            late += N * sum(states.values())
    
    # answer: 1918080160
    return late + sum(states.values())

def pe192(N=100000):
    from math import sqrt
    from fractions import Fraction

    def sign(a, b, n):
        """return the sign of a/b - sqrt(n)"""
        return a * a - n * b * b > 0

    def cfracDenom(n, sn, precision=10**12):
        """
        Basically continuous fraction convergence.
        It still need some tricks to find the best approximation
        between two known best approximations. For details, see
        https://en.wikipedia.org/wiki/Continued_fraction#Best_rational_approximations
        """

        p, q, a = 0, 1, int(sn), 
        x1, x2, y1, y2 = 1, a, 0, 1
        while y2 < precision:
            p = a * q - p
            q = Fraction(n - p * p, q)
            a = int((p + sn) / q)
            x0, x1, x2 = x1, x2, a*x2+x1
            y0, y1, y2 = y1, y2, a*y2+y1

        ybest = y1
        if a % 2 == 0:
            xtmp = (a/2) * x1 + x0
            ytmp = (a/2) * y1 + y0
            if ytmp < precision:
                sn1 = sign(x1, y1, n)
                sn2 = sign(xtmp, ytmp, n)
                if sn1 and sn2 and (x1 * ytmp > xtmp * y1):
                    ybest = ytmp
                elif sn1 and sign(x1*ytmp + xtmp*y1, y1*ytmp, 4*n):
                    ybest = ytmp
                elif sn2 and not sign(x1*ytmp + xtmp*y1, y1*ytmp, 4*n):
                    ybest = ytmp
                elif (not sn1) and (not sn2) and (x1 * ytmp < xtmp * y1):
                    ybest = ytmp
            else:
                return ybest

        for aa in range(a/2+1, a):
            ytmp = aa * y1 + y0
            if ytmp < precision:
                ybest = ytmp
            else:
                return ybest
        return ybest

    c = 0
    for n in xrange(2, 100001):
        sn = sqrt(n)
        if int(sn)**2 < n:
            c += cfracDenom(n, sn)
    
    # answer: 57060635927998347
    return c

def pe193():
    """
    Using a formula of square-free numbers and mobius function. 
    See reference articles PE193.
    """
    
    N = 2**50 - 1
    sN = 33554431
    c = N - N / 4
    
    p_sieve  = [True, True] + [False, True] * 16777215
    mu_sieve = [1, 1] + [-1, 1, 0, 1] * 8388607 + [-1, 1]
    for i in xrange(3, sN+1):
        if p_sieve[i]:
            c -= N / (i * i)

            k = i + i
            while k <= sN:
                p_sieve[k] = False
                if k % (i * i):
                    mu_sieve[k] = -mu_sieve[k]
                else:
                    mu_sieve[k] = 0
                k += i
        else:
            c += mu_sieve[i] * (N / (i * i))

    # answer: 684465067343069
    return c

def pe194(a=25, b=75, c=1984, m=10**8):
    """
    =-----------------------------------------------------------=
    |     Type     |             Given (a, b)                   |
    |--------------+--------------------------------------------|
    | a---b  a---b |                                            |
    | |cxc|  |cxc  | 1 * (n-1) * (n-2)                          |
    | b---a  b---a |                                            |
    |--------------+--------------------------------------------|
    | a---b  a---b |                                            |
    | |cxd|  |cxd  | 1 * (n-2) * (n-2)(n-3)                     |
    | b---a  b---a |                                            |
    |--------------+--------------------------------------------|
    | a---b  a---b |                                            |
    | |dxd|  |dxd  | 2(n-2) * (n-1) * (n-3)                     |
    | b---c  b---c |                                            |
    |--------------+--------------------------------------------|
    | a---b  a---b |                                            |
    | |dxe|  |dxe  | 2(n-2) * (n-2) * [(n-2) + (n-3)(n-3)]      |
    | b---c  b---c |                                            |
    |--------------+--------------------------------------------|
    | a---b  a---b |                                            |
    | |exe|  |exe  | (n-2)(n-3) * (n-1) * (n-4)                 |
    | c---d  c---d |                                            |
    |--------------+--------------------------------------------|
    | a---b  a---b |                                            |
    | |exf|  |exf  | (n-2)(n-3) * (n-2) * [2(n-2) + (n-4)(n-3)] |
    | c---d  c---d |                                            |
    |--------------+--------------------------------------------|
    | a---b        |                                            |
    | |dxd         | (n-2) * (n-1)(n-3)                         |
    | c---b        |                                            |
    |--------------+--------------------------------------------|
    | a---b        |                                            |
    | |dxe         | (n-2) * (n-2) * [(n-1) + (n-2)(n-3)]       |
    | c---b        |                                            |
    =-----------------------------------------------------------=

    Summary:
    A(n): (n - 2) * (n^4 - 7 n^3 + 20 n^2 - 29 n + 19)
    B(n): (n - 2)^3 * (n^2 - 2 n + 3)
    N(a, b, c) = c*(c-1) * C(a, a+b) * A(n)^a * B(n)^b
    """
    
    def A(n):
        return (n-2) * (n*n*n*n - 7*n*n*n + 20*n*n - 29*n + 19)

    def B(n):
        return (n-2)*(n-2)*(n-2) * (n*n - 2*n + 3)
    
    count = c * (c-1)
    count %= m
    count *= (pec.C(a+b, a) % m)
    count %= m
    count *= pow(A(c), a, m) * pow(B(c), b, m)
    count %= m
    
    # answer: 61190912
    return count