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