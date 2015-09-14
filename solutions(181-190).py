# -*- coding: utf-8 -*-

"""
solutions(181-190).py

Some interesting solutions for problems 181-190 in Project Euler

@author: Jasper Wu
"""

import math

import ProjectEuler.combinatorics as pec
import ProjectEuler.formulas as pef
import ProjectEuler.prime as pep


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

def pe185():
    """
    Modelling as an integer programming problem.
    Then using PuLP to solve it. It's really fast, just 0.24 seconds. 
    For details, see https://pythonhosted.org/PuLP/index.html
    """
    
    from pulp import LpProblem, LpVariable, LpMinimize, LpInteger, lpSum, value

    constraints = [
        ('2321386104303845', 0),
        ('3847439647293047', 1),
        ('3174248439465858', 1),
        ('8157356344118483', 1),
        ('6375711915077050', 1),
        ('6913859173121360', 1),
        ('4895722652190306', 1),
        ('5616185650518293', 2),
        ('4513559094146117', 2),
        ('2615250744386899', 2),
        ('6442889055042768', 2),
        ('2326509471271448', 2),
        ('5251583379644322', 2),
        ('2659862637316867', 2),
        ('5855462940810587', 3),
        ('9742855507068353', 3),
        ('4296849643607543', 3),
        ('7890971548908067', 3),
        ('8690095851526254', 3),
        ('1748270476758276', 3),
        ('3041631117224635', 3),
        ('1841236454324589', 3)
    ]

    VALs = map(str, range(10))
    LOCs = map(str, range(16))
    choices = LpVariable.dicts("Choice", (LOCs, VALs), 0, 1, LpInteger)

    prob = LpProblem("pe185", LpMinimize)
    prob += 0, "Arbitrary Objective Function"

    for s in LOCs:
        prob += lpSum([choices[s][v] for v in VALs]) == 1, ""

    for c, n in constraints:
        prob += lpSum([choices[str(i)][v] for i,v in enumerate(c)]) == n, ""

    prob.writeLP("pe185.lp")
    prob.solve()
    res = int(''.join(v for s in LOCs for v in VALs if value(choices[s][v])))

    # answer: 4640261571849533
    return res

def pe186(PM=524287):
    """
    Union-Find
    """
    
    idxmap = {i:i for i in xrange(1000000)}
    size   = {i:1 for i in xrange(1000000)}
    
    def find(p):
        q = idxmap[p]
        if q != p:
            q = idxmap[p] = find(q)
        return q
    
    def pairgen():
        cache, flag = [], True
        for d in range(1, 56):
            cache.append((100003 - 200003*d + 300007*d*d*d) % 1000000)
            if d % 2 == 0:
                yield cache[-2:]
            
        while 1:
            cache.append((cache[0] + cache[31]) % 1000000)
            del cache[0]
            
            flag = not flag
            if flag:
                yield cache[-2:]
    
    PG = pairgen()
    
    count = 0
    while size[find(PM)] < 990000:
        a, b = PG.next()
        if a != b:
            count += 1
            a, b = find(a), find(b)
            if a == b:
                pass
            elif size[a] < size[b]:
                idxmap[a] = b
                size[b] += size[a]
            else:
                idxmap[b] = a
                size[a] += size[b]
                
    # answer: 2325629
    return count

def pe187(N=100000000):
    """
    Sieve method to enumerate primes up to N/2
    """

    from collections import OrderedDict

    UB = N/2
    cdict = OrderedDict()
    for p in pep.P10K:
        cdict[p] = 0

    sieve = [1] * UB
    for n in range(2, UB):
        if sieve[n]:
            x = n + n
            while x < UB:
                sieve[x] = 0
                x += n

            for p in cdict:
                if p <= min(n, int(N/n)):
                    cdict[p] += 1
                else:
                    break

    # answer: 17427258
    return sum(cdict.values())

def pe188():
    """
    By Euler Theorem, phi(100000000) = 40000000,
    hence 1777^40000000 = 1 (mod 100000000).
    Then just power modular.
    """
    
    n = 1777
    for _ in range(1854):
        n = pow(1777, n, 40000000)
        
    # answer: 95962097
    return pow(1777, n, 100000000)

def pe189(N=8):
    """
    For a row with maximum length n (a1 ~ an, b1 ~ b(n-1))
           .-----.-----.-----.  ...
          / \ b1/ \ b2/ \ b3/ \
         / a1\ / a2\ / a3\ / a4\
        .-----.-----.-----.-----.  ...
    1) Enumerate all possible [a1, a2, ..., an], 3**n kinds.
    2) Fill [b1, b2, ..., b(n-1)], at most 2**(n-1) kinds.
    3) Sum up all ways can be paired with [b1, b2, ..., b(n-1)], record it to state [a1, a2, ..., an]
    4) Use all states [a1, a2, ..., an], generate all possible pairing state [c1, c2, ..., cn], at most 2**n kinds. Then record map [c...] -> [[a...], [a...], ...], which is used in 3) of next loop.  
    """

    def rowgen(d):
        if d == 1:
            for i in 'RGB':
                yield i
        else:
            for i in 'RGB':
                for remains in rowgen(d-1):
                    yield i + remains

    def rowfill(row):
        if len(row) == 2:
            for i in [s for s in 'RGB' if s != row[0] and s != row[1]]: 
                yield i
        else:
            for i in [s for s in 'RGB' if s != row[0] and s != row[1]]:
                for remains in rowfill(row[1:]):
                    yield i + remains

    def rowavail(row):
        if len(row) == 1:
            for i in [s for s in 'RGB' if s != row[0]]: 
                yield i
        else:
            for i in [s for s in 'RGB' if s != row[0]]:
                for remains in rowavail(row[1:]):
                    yield i + remains                
                    
    states = {'R': 1, 'G': 1, 'B': 1}
    avails = {'R': 2, 'G': 2, 'B': 2}
    for d in range(2, N):
        new_stat = {}
        for bottom in rowgen(d):
            for top in rowfill(bottom):
                if bottom in new_stat:
                    new_stat[bottom] += avails[top]
                else:
                    new_stat[bottom] = avails[top]
        states = new_stat
        
        avails = {}
        for r1,n in states.iteritems():
            for r2 in rowavail(r1):
                if r2 in avails:
                    avails[r2] += n
                else:
                    avails[r2] = n
    
    c = 0
    for bottom in rowgen(N):
        for top in rowfill(bottom):
            c += avails[top]

    # answer: 10834893628237824
    return c

def pe190():
    """
    Pm reaches max when x1 = x2 / 2 = x3 / 3 = ... = xm / m.
    Hence x1 = 2 / (m + 1), then the close-form of Pm can be deducted.
    """
    
    c = 0
    for m in range(2, 16):
        x = 1.
        for i in range(1, m+1):
            x *= (2. * i / (m+1))**i
        c += int(x)
    
    # answer: 371048281
    return c