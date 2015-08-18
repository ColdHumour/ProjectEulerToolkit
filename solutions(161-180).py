# -*- coding: utf-8 -*-

"""
solutions(161-180).py

Some interesting solutions for problems 161-180 in Project Euler

@author: Jasper Wu
"""

import ProjectEuler.formulas as pef


def pe161():
    """
    Dynamic programming. 
    Number 0-107 from upleft to bottomright row by row. Then a subset of (0~107) is a possible state. Define possible states at i as all subsets containing i. 
    For all possible states at (i-1), try to add 6 triominoes to them to cover the ith grid, then record the new state as possible states at i, and accumulate counters.
    """

    triominoes = [
        lambda x: (x/9 < 10) and {x, x+9, x+18},
        lambda x: (x%9 < 7)  and {x, x+1, x+2},
        lambda x: (x/9 < 11) and (x%9 < 8) and {x, x+1, x+9},
        lambda x: (x/9 < 11) and (x%9 < 8) and {x, x+1, x+10},
        lambda x: (x/9 < 11) and (x%9 < 8) and {x, x+9, x+10},
        lambda x: (x/9 < 11) and (x%9) and {x, x+8, x+9},
    ]
    
    def transform(state, n):
        tm_avail = filter(None, [t(n) for t in triominoes])
        for tm in tm_avail:
            if not (state & tm):
                new_state = state | tm
                i = n
                while i+1 in new_state:
                    i += 1
                yield new_state, i+1
    
    all_states = {tuple(): [0, 1]}
    for i in range(106):
        temp = {}
        for state, (j, c) in all_states.iteritems():
            if j != i:
                if state in temp:
                    temp[state][1] += c
                else:
                    temp[state] = [j, c]
            else:
                for state_new, k in transform(set(state), j):
                    state_new = tuple(sorted(state_new))
                    if state_new in temp:
                        temp[state_new][1] += c
                    else:
                        temp[state_new] = [k, c]
        all_states = temp
    
    # answer: 20574308184277971
    return all_states.values()[0][-1]

def pe162():
    """
    Easy combinatorics problem. Work out a close-form formula.
    Amount of strings with length n and 0, 1, A appears at least once: 
    F(n) = 15*16**(n-1) - 43*15**(n-1) + 41*14**(n-1) - 13**n, n >= 3
    """

    c = 0
    for n in range(3, 17):
        c += 15*16**(n-1) - 43*15**(n-1) + 41*14**(n-1) - 13**n

    slist = '0123456789ABCDEF'
    s = ''
    while c:
        c, q = divmod(c, 16)
        s = slist[q] + s

    # answer: 3D58725572C62302
    return s

def pe163(n=36):
    """
    Divide and conquer for all kinds of triangles.
    Can summarize out some close-form formula. See reference articles.
    """

    def backdiff3_120(n):
        return 3 * ( 1 + [0,1,1,1][n % 4] + 2*[0,1][n % 2] )

    def backdiff3_90(n):
        return 6 * ( 1 + [1,2][n % 2] + [0,1,1][n % 3] + [0,2,1,1,2][n % 5] )

    def backdiff3_60(n):
        return 1 + [0,1][n % 2] + 6 * [0,1,1][n % 3] + 2 * (n>1)

    d0 = d1 = d2 = 0
    for i in range(n):
        d3 = backdiff3_120(i) + backdiff3_60(i) + backdiff3_90(i)
        d2 = d2 + d3
        d1 = d1 + d2
        d0 = d0 + d1
    
    # answer: 343047
    return d0

def pe164(limit=20):
    """
    Dynamic programming. 
    Define possible states at i as all left two digits when the number is i-digit. 
    For all possible states at (i-1), try to expand to i-digit number while satisfying constraints, then record the new state as possible states at i, and accumulate counters.
    """

    situation = {}
    for i in range(10):
        for j in range(i, 10-i):
            situation[(i, j)] = situation[(j, i)] = 1

    for _ in range(3, limit):
        temp = {}
        for (i, j), n in situation.iteritems():
            for k in range(10-i-j):
                if (j, k) in temp:
                    temp[(j, k)] += n
                else:
                    temp[(j, k)] = n
        situation = temp

    c = 0
    for (i, j), n in situation.iteritems():  
        if i + j < 9: 
            for k in range(1, 10-i-j):
                c += n

    # answer: 378158756814587
    return c

def pe165():
    """
    Brute-force 
    """

    def intersection(line1, line2):
        x1, y1, x2, y2 = line1
        x3, y3, x4, y4 = line2

        a, b = x2 - x1, y2 - y1
        c, d = x3 - x1, y3 - y1
        e, f = x4 - x1, y4 - y1

        t1 = c * b - a * d
        t2 = e * b - a * f

        if t1 * t2 >= 0:
            return 0

        g, h = x4 - x2, y4 - y2
        i, j = x4 - x3, y4 - y3

        t1 = e * j - i * f
        t2 = g * j - i * h

        if t1 * t2 >= 0:
            return 0
        
        denom = a * j - b * i
        ynum = b * c * j + a * j * y1 - b * i * y3
        xnum = a * (c * j - d * i) + x1 * denom
        g = pef.ggcd((denom, ynum, xnum))
        return (ynum/g, xnum/g, denom/g)
    
    lines, t, c = [], [], set()
    s = 290797
    for _ in xrange(20000):
        s = (s*s) % 50515093
        t.append(s % 500)
        if len(t) == 4:
            for l in lines:
                c.add(intersection(t, l))
            lines.append(tuple(t))
            t = []
    
    # answer: 2868868
    return len(c) - 1