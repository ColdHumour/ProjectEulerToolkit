# -*- coding: utf-8 -*-

"""
solutions(201-210).py

Some interesting solutions for problems 201-210 in Project Euler

@author: Jasper Wu
"""

import math

import ProjectEuler.combinatorics as pec
import ProjectEuler.formulas as pef
import ProjectEuler.prime as pep


def pe201():
    """
    states[i] := all sums using i numbers among nlist[:i]
    For states[i], if m is unique, then states[i][m] = 1, else >1.
    """
    
    nlist = [n*n for n in range(1, 101)]
    d = len(nlist) / 2
    states = [{0: 1}] + [{} for _ in range(d)]

    for i,n in enumerate(nlist[:-1]):
        for k in range(min(d-1, i), max(0, i-d)-1, -1):
            for m0, v in states[k].iteritems():
                m1 = m0 + n
                if v >= 2:
                    states[k+1][m1] = 2
                else:    
                    states[k+1][m1] = states[k+1].get(m1, 0) + 1

    n = nlist[-1]
    for m0,v in states[d-1].iteritems():
        m1 = m0 + n
        if v >= 2:
            states[d][m1] = 2
        else:    
            states[d][m1] = states[d].get(m1, 0) + 1

    # answer: 115039000
    return sum(k for k,v in states[-1].iteritems() if v==1)

def pe202():
    """
     | \| \| \| \| \| \|   Expand triangle mirrors infinitely, then every 
     +--+--+--+--+--+--+-  possible path is a C-C line from bottomleft to
    C|\B|\A|\C|\B|\A|\C|\  somewhere upright and doesn't cross any vertex.
     | \| \| \| \| \| \|   
     +--+--+--+--+--+--+-  From the table can easity see that C-C lines are
    B|\A|\C|\B|\A|\C|\B|\  symmetric and lie on y = x + 3k. Without loss of
     | \| \| \| \| \| \|   generality, we just need to count all (a, b) with
     +--+--+--+--+--+--+-  a > b.
    A|\C|\B|\A|\C|\B|\A|\  
     | \| \| \| \| \| \|   Also, from (0, 0) to (a, b), the amount of cross
     +--+--+--+--+--+--+-  points with horizontal, vertical and diagonal
    C  B  A  C  B  A  C    lines is (a-1) + (b-1) + (a+b-1) = 2a + 2b - 3.
    
    Therefore to bounce N times, we need to find all (a, b) that:
        (1) a + b = (N + 3) / 2
        (2) gcd(a, b) = 1
        (3) a - b = 3k
    
    Recording to reference articles PE202, the number of solutions of (1), (2)
    and (3) is:
    
        sum m(d) * f(n/d) for all d satisfying d|n
    
    where m(x) if mobius function, f(x) = floor(x/2) - floor(x/3)
    """
    
    def f(x):
        return x / 2 - x / 3
    
    N = 12017639147
    N = (N + 3) / 2 # 6008819575
    
    # 6008819575 = 5*5 * 11 * 17 * 23 * 29 * 41 * 47
    pfactors = (5, 11, 17, 23, 29, 41, 47)
    count = f(N)
    for i in range(1, len(pfactors)+1):
        for pcomb in pec.combinations(pfactors, i):
            d = pef.cprod(pcomb)
            if i % 2:
                count -= f(N / d)
            else:
                count += f(N / d)
    
    # answer: 1209002624
    return count * 2

def pe203(N=51):
    """
    Use prime divisors to decide whether a combination number is square-free.
    """
    
    def oprfac(n1, n2, direction=1):
        for n,a in n2.iteritems():
            if n in n1:
                n1[n] += direction * a
            else:
                n1[n] = a
        return n1

    NFACS = {n: dict(pep.primeDivisorDecomp(n)) for n in range(2, N)}

    def squarefreeC(n):
        c, f, sfc = n, NFACS[n], set()
        if all(v==1 for v in f.values()):
            sfc.add(c)

        for i in range(2, n/2+1):
            c = c * (n+1-i) / i
            f = oprfac(f, NFACS[n+1-i], 1)
            f = oprfac(f, NFACS[i], -1)
            if all(v==1 for v in f.values()):
                sfc.add(c)
        return sfc

    distinct_sf = {1, 2}
    for i in xrange(3, N):
        distinct_sf |= squarefreeC(i)

    # answer: 34029210557338
    return sum(list(distinct_sf))