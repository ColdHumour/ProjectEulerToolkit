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
     | \| \| \| \| \| \|   Expand triangle mirrors infinitely,
     +--+--+--+--+--+--+-  then every possible path is a C-C
    C|\B|\A|\C|\B|\A|\C|\  line from bottomleft, and doesn't
     | \| \| \| \| \| \|   cross any vertex.
     +--+--+--+--+--+--+-
    B|\A|\C|\B|\A|\C|\B|\  From the table can easity see that
     | \| \| \| \| \| \|   C-C lines are symmetric and lie on
     +--+--+--+--+--+--+-  y = x + 3k.
    A|\C|\B|\A|\C|\B|\A|\  
     | \| \| \| \| \| \|   Also, from (0, 0) to (a, b), there 
     +--+--+--+--+--+--+-  are following amount of cross points,
    C  B  A  C  B  A  C    to horizontal, vertical and diagonal
    
    lines: (a-1) + (b-1) + (a+b-1) = 2a + 2b - 3 
    
    Therefore to bounce N times, we need to find all (a, b) that:
    (1) a + b = (N + 3) / 2
    (2) gcd(a, b) = 1
    (3) a - b = 3k
    
    Recording to reference articles PE202, the number of solutions
    of (1), (2) and (3) is:
    
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