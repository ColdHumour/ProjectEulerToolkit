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

def pe192(N=100000):
    from math import sqrt
    from fractions import Fraction

    def sign(a, b, n):
        """return the sign of a/b - sqrt(n)"""
        return a * a - n * b * b > 0

    def cfracDenom(n, sn, precision=10**12):
        """
        Basically continuous fraction convergence.
        It still need some tricks to find the best approximation
        between two known best approximations. For details, see
        https://en.wikipedia.org/wiki/Continued_fraction#Best_rational_approximations
        """

        p, q, a = 0, 1, int(sn), 
        x1, x2, y1, y2 = 1, a, 0, 1
        while y2 < precision:
            p = a * q - p
            q = Fraction(n - p * p, q)
            a = int((p + sn) / q)
            x0, x1, x2 = x1, x2, a*x2+x1
            y0, y1, y2 = y1, y2, a*y2+y1

        ybest = y1
        if a % 2 == 0:
            xtmp = (a/2) * x1 + x0
            ytmp = (a/2) * y1 + y0
            if ytmp < precision:
                sn1 = sign(x1, y1, n)
                sn2 = sign(xtmp, ytmp, n)
                if sn1 and sn2 and (x1 * ytmp > xtmp * y1):
                    ybest = ytmp
                elif sn1 and sign(x1*ytmp + xtmp*y1, y1*ytmp, 4*n):
                    ybest = ytmp
                elif sn2 and not sign(x1*ytmp + xtmp*y1, y1*ytmp, 4*n):
                    ybest = ytmp
                elif (not sn1) and (not sn2) and (x1 * ytmp < xtmp * y1):
                    ybest = ytmp
            else:
                return ybest

        for aa in range(a/2+1, a):
            ytmp = aa * y1 + y0
            if ytmp < precision:
                ybest = ytmp
            else:
                return ybest
        return ybest

    c = 0
    for n in xrange(2, 100001):
        sn = sqrt(n)
        if int(sn)**2 < n:
            c += cfracDenom(n, sn)
    
    # answer: 57060635927998347
    return c

def pe193():
    """
    Using a formula of square-free numbers and mobius function. 
    See reference articles PE193.
    """
    
    N = 2**50 - 1
    sN = 33554431
    c = N - N / 4
    
    p_sieve  = [True, True] + [False, True] * 16777215
    mu_sieve = [1, 1] + [-1, 1, 0, 1] * 8388607 + [-1, 1]
    for i in xrange(3, sN+1):
        if p_sieve[i]:
            c -= N / (i * i)

            k = i + i
            while k <= sN:
                p_sieve[k] = False
                if k % (i * i):
                    mu_sieve[k] = -mu_sieve[k]
                else:
                    mu_sieve[k] = 0
                k += i
        else:
            c += mu_sieve[i] * (N / (i * i))

    # answer: 684465067343069
    return c

def pe194(a=25, b=75, c=1984, m=10**8):
    """
    =-----------------------------------------------------------=
    |     Type     |             Given (a, b)                   |
    |--------------+--------------------------------------------|
    | a---b  a---b |                                            |
    | |cxc|  |cxc  | 1 * (n-1) * (n-2)                          |
    | b---a  b---a |                                            |
    |--------------+--------------------------------------------|
    | a---b  a---b |                                            |
    | |cxd|  |cxd  | 1 * (n-2) * (n-2)(n-3)                     |
    | b---a  b---a |                                            |
    |--------------+--------------------------------------------|
    | a---b  a---b |                                            |
    | |dxd|  |dxd  | 2(n-2) * (n-1) * (n-3)                     |
    | b---c  b---c |                                            |
    |--------------+--------------------------------------------|
    | a---b  a---b |                                            |
    | |dxe|  |dxe  | 2(n-2) * (n-2) * [(n-2) + (n-3)(n-3)]      |
    | b---c  b---c |                                            |
    |--------------+--------------------------------------------|
    | a---b  a---b |                                            |
    | |exe|  |exe  | (n-2)(n-3) * (n-1) * (n-4)                 |
    | c---d  c---d |                                            |
    |--------------+--------------------------------------------|
    | a---b  a---b |                                            |
    | |exf|  |exf  | (n-2)(n-3) * (n-2) * [2(n-2) + (n-4)(n-3)] |
    | c---d  c---d |                                            |
    |--------------+--------------------------------------------|
    | a---b        |                                            |
    | |dxd         | (n-2) * (n-1)(n-3)                         |
    | c---b        |                                            |
    |--------------+--------------------------------------------|
    | a---b        |                                            |
    | |dxe         | (n-2) * (n-2) * [(n-1) + (n-2)(n-3)]       |
    | c---b        |                                            |
    =-----------------------------------------------------------=

    Summary:
    A(n): (n - 2) * (n^4 - 7 n^3 + 20 n^2 - 29 n + 19)
    B(n): (n - 2)^3 * (n^2 - 2 n + 3)
    N(a, b, c) = c*(c-1) * C(a, a+b) * A(n)^a * B(n)^b
    """
    
    def A(n):
        return (n-2) * (n*n*n*n - 7*n*n*n + 20*n*n - 29*n + 19)

    def B(n):
        return (n-2)*(n-2)*(n-2) * (n*n - 2*n + 3)
    
    count = c * (c-1)
    count %= m
    count *= (pec.C(a+b, a) % m)
    count %= m
    count *= pow(A(c), a, m) * pow(B(c), b, m)
    count %= m
    
    # answer: 61190912
    return count

def pe195(R=1053779):
    def CWgen(a, b):
        """Using Calkin-Wilf Tree to generate co-prime pairs."""

        return [(a, a+b), (a+b, b)]

    def rgen(m, n):
        """
        Using co-prime pairs to generate primitive Eisenstein triplets.
        See https://en.wikipedia.org/wiki/Integer_triangle
        """

        if m > 2 * n:
            a, b = 2*m*n - n*n, m*m - n*n
        else:
            a, b = m*m - n*n, 2*m*n - n*n

        g = pef.ggcd((a, b, m*m - m*n + n*n))
        return g, sqrt(3) / 2 * (m-n) * n / g

    count, to_search = 0, [(1, 3), (3, 2)]
    triangles = []

    while to_search:
        m, n = to_search.pop()
        if m < n:
            m, n = n, m

        g, r = rgen(m, n)
        if g * r > R * 3:
            continue  
        elif r <= R:
            count += int(R / r)
        to_search += CWgen(m, n)

    # answer: 75085391
    return count / 2

def pe196():
    """
    Basically enumerate all primes from row n-2 to row n+2.
    Using 6k+1, 6k+5 rule and gmpy2 to accelerate.
    Details see https://gmpy2.readthedocs.org/en/latest/index.html
    """

    from gmpy2 import is_prime
    
    MIDDLE = {
        (1, 0): lambda x,r: {x-r: (x-2, ), x+r: (x+2*r, )},
        (1, 1): lambda x,r: {x-r+1: (x-2*r+2, ), x+r-1: (x-2, )},
        (1, 2): lambda x,r: {x-r: (x-2*r+2, x-2), x-r+2: (x-2*r+4, )},
        (1, 3): lambda x,r: {x-r+1: (x-2*r+4, ), x+r+1: ()},
        (1, 4): lambda x,r: {x-r+2: (), x+r: (x+2*r+2, )},
        (1, 5): lambda x,r: {x+r-1: (x-2, x+2*r), x+r+1: (x+2*r+2, )},
        (5, 0): lambda x,r: {x-r: (x-2*r+2, ), x-r+2: (x+2, ), x+r: (x+2*r, x+2*r+2)},
        (5, 1): lambda x,r: {x-r+1: (x-2*r+2, x-2*r+4), x+r-1: (x+2*r, ), x+r+1: (x+2, )},
        (5, 2): lambda x,r: {x-r+2: (x-2*r+4, x+2), x+r: (x+2*r+2, )},
        (5, 3): lambda x,r: {x+r-1: (x+2*r, )},
        (5, 4): lambda x,r: {x-r: (x-2*r+2, )},
        (5, 5): lambda x,r: {x-r+1: (x-2*r+4, ), x+r+1: (x+2, x+2*r+2)}
    }

    def S(n):
        res, xmin = 0, (n-1) * n / 2
        for i in xrange(1+(xmin%2), n+1, 2):
            x = xmin + i
            if x % 3 and x % 7 and x % 11 and is_prime(x):
                univ = MIDDLE[(x%6, n%6)](x, n)        
                check = filter(is_prime, univ.keys())
                d = len(check)

                if d > 1:
                    res += x
                elif d:
                    try:
                        for y in check:
                            for z in univ[y]:
                                if is_prime(z):
                                    res += x
                                    raise ValueError
                    except:
                        pass
        return res
    
    # answer: 322303240771079935
    return S(5678027) + S(7208785)

def pe197():
    """
    Oscillating between 1.029461839 and 0.681175878. Since
    f(0.681175878) = 1.029461839
    f(1.029461839) = 0.681175878
    """
    
    def f(x):
        return int(2**(30.403243784-x*x)) * 1e-9

    x = 1
    for _ in xrange(1000):
        x = f(x)
        
    # answer: 1.710637717
    return x + f(x)

def pe198(N=10**8):
    """
    Every number equals to (a + b) / 2 is an ambiguous number, 
    where a and b are two adjacent leaves of Farey sequence.
    Farey sequence can be implemented as a binary tree, beginning
    with [(0, 1), (1, 1)] and insert (a+c, b+d) between (a, b) 
    and (c, d). When keeping b+d <= n, then it's all reduced 
    rational numbers 0 <= p/q <= 1 with 0 <= p <= q <= n.
    """

    to_search = [(0, 1, 1, 1)]
    count = 0
    while to_search:
        a, b, c, d = to_search.pop()
        e, f = a+c, b+d
        
        x, y = a*d+b*c, 2*b*d
        if x > 1:
            g = gcd(x, y)
            x /= g
            y /= g

        if a * 100 < b and y <= N:
            to_search.append((a, b, e, f))
            if e * 100 < f:
                count += 1
                to_search.append((e, f, c, d))
    
    # answer: 52374475
    return count

def pe199(N=10):
    from math import sqrt

    def DT(k1, k2, k3):
        """
        Descartes's Theorem, for details see 
        https://en.wikipedia.org/wiki/Descartes'_theorem
        
        Notice that the curvature of a circle is -1/r when 
        inscribed and 1/r when tangent.
        """

        return k1 + k2 + k3 + 2 * sqrt(k1*k2 + k1*k3 + k2*k3)

    # special cases
    k1 = 1 / (2*sqrt(3) - 3) # 3 largest circles
    k2 = DT(-1, k1, k1)      # 3 2nd largest circles
    k3 = DT(k1, k1, k1)      # 1 circle in the middle
    s = 3 / (k1 * k1) + 3 / (k2 * k2) + 1 / (k3 * k3)

    states1 = [(k1, k2)] # inscribed the outer circle
    states2 = [(k1, k1, k2), (k1, k1, k3)] # 3 kissing circles
    for _ in range(N-1):
        new_stat1, new_stat2 = [], []
        for k1, k2 in states1:
            k4 = DT(-1, k1, k2)
            s += 6 / (k4 * k4)
            new_stat1.extend([(k1, k4), (k2, k4)])
            new_stat2.append((k1, k2, k4))

        for k1, k2, k3 in states2:
            k4 = DT(k1, k2, k3)
            if k1 == k2:
                s += 3 / (k4 * k4)
                new_stat2.extend([(k1, k1, k4), (k1, k3, k4)])
            else:
                s += 6 / (k4 * k4)
                new_stat2.extend([(k1, k2, k4), 
                                  (k1, k3, k4),
                                  (k2, k3, k4)])
        states1, states2 = new_stat1, new_stat2
    
    # answer: 0.00396087
    return round(1 - s, 8)