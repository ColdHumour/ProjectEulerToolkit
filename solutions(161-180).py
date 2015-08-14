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