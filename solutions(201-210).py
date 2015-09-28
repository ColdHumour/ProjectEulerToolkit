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