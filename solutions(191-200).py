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