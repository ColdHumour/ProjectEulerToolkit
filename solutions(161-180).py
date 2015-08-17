# -*- coding: utf-8 -*-

"""
solutions(161-180).py

Some interesting solutions for problems 161-180 in Project Euler

@author: Jasper Wu
"""


def pe161():
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