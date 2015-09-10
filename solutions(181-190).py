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

def pe182():
    """
    Exact number of unconcealed messages: 
        (1 + gcd(e-1, p-1)) * (1 + gcd(e-1, p-1))
        
    (1009 - 1) * (3643 - 1) = 2^5 * 3^3 * 7 * 607
    Then e = 6 * k + 1 or 5, e % 7 != 0 and e % 607 != 0
    """
    
    from fractions import gcd

    p, q = 1009, 3643
    emin, edict = p * q, {}
    for k in xrange(611856):
        for e in (6 * k + 1, 6 * k + 5):
            if e % 7 and e % 607:
                u = (1 + gcd(e-1, p-1)) * (1 + gcd(e-1, q-1))
                emin = min(emin, u)
                if u in edict:
                    edict[u] += e
                else:
                    edict[u] = e
    
    # answer: 399788195976
    return edict[emin]

def pe183(N=10000):
    """
    When N < delta(k), P(N, k) - P(N, k+1) > 0 always holds.
    Hence for all N between delta(k-1) and delta(k), P(N, k) is max.
    """

    from fractions import gcd

    def delta(k):
        return (1 + 1./k)**k * (k + 1)

    kdict = {}
    n, k = 5, 2
    while n <= N:
        m = int(delta(k))
        kdict[k] = range(n, min(N, m)+1)
        n = m + 1
        k += 1

    c = 0
    for k,nlist in kdict.iteritems():
        for n in nlist:
            denom = k / gcd(n, k)
            while denom % 2 == 0:
                denom /= 2
            while denom % 5 == 0:
                denom /= 5

            c += -n if denom == 1 else n
    
    # answer: 48861552
    return c

def pe184(r=105):
    """
    Find all points in first quardrant, then group them by slope (a.k.a. lines).
    Count how many points between each line and y-axis.
    There 5 different types of triangles. Count them out.
    """
    
    from fractions import gcd

    def count(p, r):
        xp, yp = p

        xmax = int(math.sqrt(1. * r*r / (xp*xp + yp*yp)) * xp)
        yk = 1. * yp * xmax / xp
        if xmax*xmax + yk*yk < r*r:
            xmax += 1

        c = 0
        for x in range(1, xmax):
            ymax = int(math.sqrt(r*r - x*x))
            if x*x + ymax*ymax == r*r:
                ymax -= 1

            ymin = int(1. * x * yp / xp) + 1
            c += ymax - ymin + 1
        return c

    nq1, nline = 0, {}
    for x in range(1, r):
        for y in range(1, int(math.sqrt(r*r - x*x)) + 1):
            if x*x + y*y < r*r:
                nq1 += 1
                g = gcd(x, y)
                xr = x / g
                yr = y / g

                if (xr, yr) in nline:
                    nline[(xr, yr)] += 1
                else:
                    nline[(xr, yr)] = 1

    # a(*4): one apex on y+, one apex on x+
    # b(*8): one apex on y+, one apex in 1st quardrant
    # c(*4): one apex on y+, one apex in 2nd quardrant, one apex in 3rd quardrant
    # d(*4): two apexes in first quardrant
    # e(*2): three apexes in different quardrants

    a = (r-1) * (r-1) * nq1
    b = c = d = e = 0

    nty, ntx = {}, {}
    for (x, y), n in nline.iteritems():
        nty[(x, y)] = count((x, y), r)
        ntx[(x, y)] = nq1 - n - nty[(x, y)]

        b += (r-1) * n * nty[(x, y)]
        c += (r-1) * n * nq1
        e += 2 * n * n * ntx[(x, y)]

    for (x1, y1), (x2, y2) in pec.combinations(nline.keys(), 2):
        if 1.*y2/x2 > 1.*y1/x1:
            x1, y1, x2, y2 = x2, y2, x1, y1

        coeff = nline[(x1, y1)] * nline[(x2, y2)]
        n = nty[(x2, y2)] - nty[(x1, y1)] - nline[(x1, y1)]

        d += coeff * n
        e += 2 * coeff * (n + nline[(x2, y2)] + 2 * ntx[(x2, y2)])

    # answer: 1725323624056
    return a*4 + b*8 + c*4 + d*4 + e*2