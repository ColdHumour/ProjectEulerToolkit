# -*- coding: utf-8 -*-

"""
solutions(181-190).py

Some interesting solutions for problems 181-190 in Project Euler

@author: Jasper Wu
"""

import math

import ProjectEuler.combinatorics as pec
import ProjectEuler.formulas as pef


def pe181_recursive():
    # P(n, k): ways of partition n into k parts, order doesn't matter
    #
    # p(n, k) = p(n-1, k-1) + p(n-k, k) for n > 0, k > 0
    #
    # p(n, 1) = 1, p(n, k) = 0 for all k > n
    #
    # (can be calculated ahead)
    
    Pdict = [[0], [0, 1], [0, 1, 1]]
    for n in range(3, 101):
        m = n/2
        r = [0]
        for k in range(1, m+1):
            if k == 1:
                r.append(Pdict[n-k][k])
            else:
                r.append(Pdict[n-1][k-1] + Pdict[n-k][k])
        Pdict.append(r + Pdict[n-1][m:])

    def P(n, k):
        if k > n:
            return 0
        elif k == 1 or k == n:
            return 1
        else:
            return Pdict[n][k]

    # A(a, b, k): partition a into m1 parts, b into m2 parts, where m1 + m2 = k
    #
    # A(a, b, k) = P(a, m) * P(b, k-m) for max(1, k-b) <= m <= min(a, k-1)
    # s.t. a + b >= k
    #
    # A(0, 0, k) = 1, A(a, 0, k) = P(a, k)
    # A(a, b, k) = A(b, a, k)
    #
    # (can be calculated ahead)
    
    Adict = {}
    for a in range(41):
        for b in range(61):
            for k in range(101):
                if a + b == 0:
                    Adict[(a, b, k)] = 1
                elif a * b == 0:
                    Adict[(a, b, k)] = P(a+b, k)
                else:
                    Adict[(a, b, k)] = sum(P(a, m) * P(b, k-m) for m in range(max(1, k-b), min(a, k-1)+1))

    def A(a, b, k):
        if a > b:
            a, b = b, a
        return Adict[(a, b, k)]

    
    # B(a, b, k): group a and b into m groups, where m <= k
    #
    # B(a, b, k) = Q(a, b, m) for 1 <= m <= min(k, a+b)
    # s.t. a + b >= k
    #
    # B(a, b, 0) = 1 if (a == b == 0) else 0
    # B(a, b, k) = B(b, a, k)
    
    Bdict, Qdict = {}, {}
    
    def B(a, b, k):
        if a + b == 0:
            return 1
        elif (a, b, k) in Bdict:
            return Bdict[(a, b, k)]
        elif a * b == 0:
            Bdict[(a, b, k)] = sum(P(a+b, m) for m in range(1, min(k, a+b)+1))
            return Bdict[(a, b, k)]
        elif a > b:
            a, b = b, a

        Bdict[(a, b, k)] = sum(Q(a, b, m) for m in range(1, min(k, a+b)+1))
        return Bdict[(a, b, k)]

    # Q(a, b, k): group a and b into k groups, which can be decomposed into 2 steps
    #             (1) group (a-i) and (b-j) using A(i, j, k-m)
    #             (2) group i and j using B(i, j, m)
    #
    # Q(a, b, k) = A(a, b, k) + 
    #              A(a-i, b-j, k-m) * B(i-m, j-m, m) for i, j, m
    # s.t. 1 <= m <= min(a, k)
    #      m <= i <= min(a, a+b-k)
    #      m <= j <= min(b, a+b-k+m-i)
    #
    # Q(0, 0, k) = 1, Q(a, b, 0) = 0 (a+b>0), Q(a, b, 1) = 1
    # Q(a, b, k) = Q(b, a, k)
    
    def Q(a, b, k):
        if a == 0 and b == 0:
            return 1
        elif k == 1:
            return 1
        elif a + b == k:
            return 1
        elif a > b:
            a, b = b, a

        if (a, b, k) in Qdict:
            return Qdict[(a, b, k)]

        if a * b == 0:
            Qdict[(a, b, k)] = P(a+b, k)
            return Qdict[(a, b, k)]

        c = A(a, b, k)
        for m in range(1, min(a, k)+1):
            for i in range(m, min(a, a+b-k) + 1):
                for j in range(m, min(b, a+b-k+m-i) + 1):
                    c += A(a-i, b-j, k-m) * B(i-m, j-m, m)

        Qdict[(a, b, k)] = c
        return c

    c = 2
    for k in range(2, 100):
        for a in range(1, 41):
            for b in range(a, 61):
                B(a, b, k)
                Q(a, b, k)
                if a == 40 and b == 60:
                    c += Q(40, 60, k)
    
    # time: 16 minutes
    # answer: 83735848679360680
    return c

def pe181_dp(a=40, b=60):
    Q = [[0 for j in range(b+1)] for i in range(a+1)]
    
    Q[0][0] = 1
    for i in range(a+1):
        for j in range(b+1):
            for k in range(i, a+1):
                for l in range(j, b+1):
                    Q[k][l] += Q[k-i][l-j]
    return Q[a][b] / 2