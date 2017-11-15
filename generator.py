# -*- coding: utf-8 -*-

"""
generator.py

Functions using to generate numbers defind by PE.
Function list:
    pe_lagged_fibo_generator

@author: Jasper Wu
"""

from collections import deque


def pe_lagged_fibo_generator(batch=1, mod=1000000):
    """Lagged Fibonacci Generator"""

    seq, out = deque([], maxlen=55), []
    while 1:
        if len(seq) < 55:
            k = len(seq) + 1
            x = (100003 - 200003 * k + 300007 * k*k*k) % mod
        else:
            x = (seq[-55] + seq[-24]) % mod

        seq.append(x)
        out.append(x)
        if len(out) == batch:
            yield out
            out = []
