# -*- coding: utf-8 -*-

"""
solutions(171-180).py

Some interesting solutions for problems 171-180 in Project Euler

@author: Jasper Wu
"""

import math

import ProjectEuler.combinatorics as pec
import ProjectEuler.formulas as pef


def pe171(N=20):
    """
    Make transformation from state to state. 
    Record states at each level.
    """
    
    digits = range(10)
    squares = [n**2 for n in range(41)]

    nsum = 45
    sdict = {n**2: (n, 1) for n in range(1, 10)}
    for _ in range(N-1):
        temp = {}
        for n, (s1, i1) in sdict.iteritems():
            for m in digits:
                t = n + m * m            
                s0, i0 = temp.get(t, (0, 0))
                temp[t] = (s0 + 10 * s1 + m * i1, i0 + i1)
        sdict = temp

        for n in squares:
            nsum += sdict.get(n, (0, 0))[0]
        nsum %= 10**9
    
    # answer: 142989277
    return nsum

def pe172(N=18):
    """
    Split 18 digits to 10 piles, each pile contains at most 3 digits.
    Get all permutations of all partitions, find out whether it contains 0.
    And then count the possible permutations of those digits.
    """
    
    from copy import deepcopy

    def countnum(partition):
        ap = pec.MP(partition) # all permutations

        c = 0
        for n in set(partition):
            # count how many permutations containing n 0s
            temp = deepcopy(partition)
            
            temp.remove(n)
            k = pec.MP([temp.count(x) for x in set(temp)])

            if n == 0: # if not containing any 0
                c += k * ap
            else:      # else get rid off numbers starting with 0
                temp.append(n-1)
                c += k * (ap - pec.MP(temp))
        return c

    c = 0
    for r in pec.allPartitions(10, N, init=0):
        if all([n < 4 for n in r]):
            c += countnum(r)
    
    # answer: 227485267000992000
    return c

def pe173():
    """
    Easy to calculate a upper bound and a counting formula.
    """
    
    c = 0
    for n in range(1, 250000):
        c += int((math.sqrt(n*n + 1000000) - n) / 2)
        
    # answer: 1572729
    return c

def pe174():
    """
    Use the result of pe173.
    """
    
    c, H = 0, {}
    for n in range(1, 250000):
        for m in range(1, int((math.sqrt(n*n + 1000000) - n) / 2) + 1):
            x = n + 2 * m
            y = x * x - n * n
            if y in H:
                H[y] += 1
            else:
                H[y] = 1
        
        # save space
        for i in range(4*n, 4*n+4):
            if H.get(i, 11) <= 10:
                c += 1
                del H[i]
                
    for v in H.values():
        c += v <= 10
    
    # answer: 209566
    return c

def pe175(a=987654321, b=123456789):
    """
    Using the amazing relationship between Calkin-Wilf Tree,
    Stern's Diatomic Series and hyperbinary representation.
    """
    
    from collections import deque
    from fractions import gcd
    
    g = gcd(a, b)
    a /= g
    b /= g
    
    s = deque([])
    while a and b:
        if a > b:
            s.appendleft(a / b)
            a %= b
        else:
            s.appendleft(b / a)
            b %= a
    
    if a:
        s[0] -= 1
        s.appendleft(1)
        
    # answer: 1,13717420,8
    return ','.join(map(str, s))

def pe176(P=47547):
    """
    N = 2**a0 * p1**a1 * p2**a2 * ... * pk**ak
    A = (2*a1 + 1) * (2*a2 + 1) * ... * (2*ak + 1) 
    P(N) = (A - 1) / 2,              if a0 = 0;
         = ((2*a0 - 1) * A - 1) / 2, if a0 > 0.
    """
    
    ODDPRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    
    def F(plist):
        return pef.cprod([2*n+1 for n in plist])

    cmin, info = 10**40, []
    for i in range(1, 100):
        for j in range(1, 10):
            for r in pec.allPartitions(j, i):
                r = filter(None, r)[::-1]
                a = F(r)
                if (a - 1) / 2 == P:
                    c = 1
                    for x, y in zip(ODDPRIMES[:len(r)], r):
                        c *= x**y
                        if c > cmin:
                            break
                    if c < cmin:
                        cmin, info = c, [0] + r
                elif (P*2 + 1) % a == 0:
                    b = (P*2 + 1) / a
                    if b % 2 == 1 and b / 2 + 1 <= 40:
                        c = 2**(b / 2 + 1)
                        for x, y in zip(ODDPRIMES[:len(r)], r):
                            c *= x**y
                            if c > cmin:
                                break
                        if c < cmin:
                            cmin, info = c, [b / 2 + 1] + r
                elif (a - 1) / 2 > 47547:
                    break

    # answer: 96818198400000 = 2**10 * 3**6 * 5**5 * 7**3 * 11**2
    return cmin, info

def pe177():
    """
       A *       For a quarilateral, if we get 4 angles BAC, DAC, ACB, ACD,
        /|\      we can set AC to 1, and calculate AB and AD by Sine Theorem.
       / | \     Then using Cosine Theorem, BD = AB**2 + AD**2 - 2*AB*AD*cosA.
    B *--+--* D  Using Cosine Theorem again, we can get angle ABD, then the
       \ | /     left 3 corners CBD, ADB, CDB can be easily calculated.
        \|/
         * C
    """
    
    from math import pi, sqrt, sin, cos, asin, acos

    D2R = pi / 180
    DIVS = {n:[(i, n-i) for i in xrange(1, n/2+1)] for n in xrange(2, 179)}

    def solve(a, b, c, d):
        # solve 4 remaining angles
        
        ra = a * D2R
        rb = b * D2R
        rc = c * D2R
        rd = d * D2R

        e, f = sin(rc)*100, sin(rd)*100 
        g, h, i = sin(ra + rc)*100, sin(rb + rd)*100, cos(ra + rb)

        AB2, AD2 = e*e*h*h, f*f*g*g
        BD2 = AB2 + AD2 - 2*e*f*g*h*i
        cosX = (AB2 + BD2 - AD2) / (2*e*h*sqrt(BD2))

        if cosX <= -1 or cosX >= 1:
            return 0

        x1 = acos(cosX) / D2R
        if abs(x1 - round(x1)) > 1e-9:
            return 0

        x1 = round(x1)
        x2 = 180 - a - c - x1
        x3 = x2 + c - b
        x4 = x1 + a - d

        # rearrange x1~x4 to get the different representation 
        # of the quadrilateral
        if x3 + x4 < x1 + x2:
            x1, x2, x3, x4 = x3, x4, x1, x2
        if x1 + x2 == x3 + x4 and min(x3, x4) < min(x1, x2):
            x1, x2, x3, x4 = x3, x4, x1, x2
        if x1 > x2:
            x1, x2, x3, x4 = x2, x1, x4, x3
        if x1 == x2 and x3 > x4:
            x3, x4 = x4, x3

        return x1, x2, x3, x4

    count, visited = 0, set()
    for n1 in xrange(2, 179):
        for n2 in xrange(n1, 179):
            for i,(a,b) in enumerate(DIVS[n1]):
                # carefully distinguish different situation
                if n1 == n2:
                    SN2 = DIVS[n1][i:]
                else:
                    SN2 = DIVS[n2]

                for c,d in SN2:
                    if a + c > 178:
                        break

                    if b + d <= 178:
                        if (a, b, c, d) in visited:
                            visited.remove((a, b, c, d))
                        else:
                            sol = solve(a, b, c, d)
                            if sol:
                                count += 1
                                if sol != (a, b, c, d):
                                    visited.add(sol)

                    # sometimes it needs to check (a, b, d, c)
                    if a < b and c < d and a + d <= 178 and b + c <= 178:
                        if (a, b, d, c) in visited:
                            visited.remove((a, b, d, c))
                        else:
                            sol = solve(a, b, d, c)
                            if sol:
                                count += 1
                                if sol != (a, b, d, c):
                                    visited.add(sol)
    assert len(visited) == 0

    # answer: 129325
    return count

def pe178(N=40):
    """
    state i: (start digit, end digit, min digit, max digit) of i-digit number
    state i -> state i+1: add one digit to the tail of i-digit numbers
    """
    
    c = 0
    state = {(i, i, i, i): 1 for i in range(1, 10)}
    for k in range(N-1):
        temp = {}
        for (i, j, a, b), n in state.iteritems():
            if j == 0:
                if (i, 1, a, b) in temp:
                    temp[(i, 1, a, b)] += n
                else:
                    temp[(i, 1, a, b)] = n
            elif j == 9:
                a = min(a, 8)
                if (i, 8, a, 9) in temp:
                    temp[(i, 8, a, 9)] += n
                else:
                    temp[(i, 8, a, 9)] = n
            else:
                aa = min(a, j-1)
                if (i, j-1, aa, b) in temp:
                    temp[(i, j-1, aa, b)] += n
                else:
                    temp[(i, j-1, aa, b)] = n

                bb = max(b, j+1)
                if (i, j+1, a, bb) in temp:
                    temp[(i, j+1, a, bb)] += n
                else:
                    temp[(i, j+1, a, bb)] = n
        state = temp

        # k = 8 means add 9 digits, therefore 10 digits already
        if k > 7:
            for (i, j, a, b), n in state.iteritems():
                if (a, b) == (0, 9):
                    c += n

    # answer: 126461847755
    return c

def pe179(N=10**7):
    """
    factors sieve method
    """
    
    c = 0
    ndiv = [2] * (N + 1)
    for i in xrange(2, N):
        if ndiv[i] == ndiv[i+1]:
            c += 1
        j = i + i
        while j < N + 1:
            ndiv[j] += 1
            j += i
    
    # answer: 986262
    return c

def pe180(order=35):
    """
    Simplify the expression to x**n + y**n = z**n.
    By Fermat Last Theorem, n = -2, -1, 1, 2. Then brute-force.
    """
    
    from fractions import Fraction
    from math import sqrt

    RNS = sorted([Fraction(a, b) for b in range(2, order+1) for a in range(1, b)])

    fset = set()
    for i,x in enumerate(RNS):
        for y in RNS[i:]:
            x1, y1 = x**(-1), y**(-1)
            x2, y2 = x*x, y*y
            x3, y3 = x2**(-1), y2**(-1)

            s = x + y

            # n = 1
            if s.numerator < s.denominator < 36:
                t = s + s
                fset.add((t.numerator, t.denominator))

            # n = -1
            z = x1 + y1
            if z.denominator < z.numerator < 36:
                t = s + z**(-1)
                fset.add((t.numerator, t.denominator))

            # n = 2
            z = x2 + y2
            a2, b2 = z.numerator, z.denominator
            a, b = int(sqrt(a2)), int(sqrt(b2))
            if a * a == a2 and b * b == b2 and a < b < 36:
                t = s + Fraction(a, b)
                fset.add((t.numerator, t.denominator))

            # n = -2
            z = x3 + y3
            a2, b2 = z.numerator, z.denominator
            a, b = int(sqrt(a2)), int(sqrt(b2))
            if a * a == a2 and b * b == b2 and b < a < 36:
                t = s + Fraction(b, a)
                fset.add((t.numerator, t.denominator))

    c = 0
    for a, b in fset:
        c += Fraction(a, b)
        
    # answer: 285196020571078987
    return c.numerator + c.denominator