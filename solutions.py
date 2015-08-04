# -*- coding: utf-8 -*-

"""
solutions.py

Some interesting solutions for problems in Project Euler

@author: Jasper Wu
"""

from math import sqrt
from . prime import isprime, atkin_sieve


def hexCoord(k):
    """
    Compute the specific coordinates in the hexagon table of ep128
    """

    i = int((3 + sqrt(12*k-15)) / 6)
    r = k - 3*i*(i-1) - 2
    j, k = divmod(r, i)
    return i, j, k

def hexNeighbors(limit):
    """
    Compute neighbors of numbers in the hexagon table of ep128
    """

    cdict = {(0,0,0): 1, (1,0,0): 2, (1,1,0): 3, (1,2,0): 4, 
                         (1,3,0): 5, (1,4,0): 6, (1,5,0): 7,}
    ndict = {1: (0,0,0), 2: (1,0,0), 3: (1,1,0), 4: (1,2,0),
                         5: (1,3,0), 6: (1,4,0), 7: (1,5,0),}
    neighbour = {1: [2, 3, 4, 5, 6, 7],
                2: [1, 3, 7, 8, 9, 19],
                3: [1, 2, 4, 9, 10, 11],
                4: [1, 3, 5, 11, 12, 13],
                5: [1, 4, 6, 13, 14, 15],
                6: [1, 5, 7, 15, 16, 17],
                7: [1, 2, 6, 17, 18, 19],}
    if limit < 7:
        return neighbour
        
    for n in range(8, limit+1):
        c = hexCoord(n)
        cdict[c], ndict[n] = n, c
    
    for n in range(8, limit+1):
        i, j, k = ndict[n]
        if j == k == 0:
            try:
                neighbour[n] = [n+1, 3*(i-1)*(i-2)+2, 3*(i+1)*(i+2)+1, 
                                3*i*(i+1)+1, 3*i*(i+1)+2, 3*i*(i+1)+3]
            except:
                pass
        elif j == 5 and k == i - 1:
            try:
                neighbour[n] = [n+1, 3*(i-1)*(i-2)+2, 3*i*(i-1)+2, 
                                cdict[(i-1,5,i-2)], cdict[(i+1,5,i-1)], cdict[(i+1,5,i)]]
            except:
                pass
        elif k == 0:
            try:
                x = cdict[(i+1,j,0)]
                neighbour[n] = [cdict[(i-1,j,0)], n-1, n+1, x-1, x, x+1]
            except:
                pass
        else:
            try:
                x = cdict[(i+1,j,0)]
                neighbour[n] = [cdict[(i-1,j,k-1)], n-1, n+1, 
                                cdict[(i-1,(j+(k==i-1))%6,(k!=i-1)*k)], 
                                cdict[(i+1,j,k)], cdict[(i+1,j,k+1)]]
            except:
                pass
    return neighbour

def pe128(index):
    """
    Using patterns found in hexNeighbors()
    """
    
    x = i = 2
    while 1:
        if isprime(6*i-1):
            if isprime(6*i+1) and isprime(12*i+5):
                x += 1
                if x == index:
                    return 3*i*(i-1)+2
            if isprime(6*i+5) and isprime(12*i-7):
                x += 1
                if x == index:
                    return 3*i*(i+1)+1
        i += 1

def pe136(n=50000000):
    """
    The problem can be reducted to following form:
    
    Find how many positive integers n less than 50m such that n = k * (4*d - k) has only one set of solution (k, d) satisfying k <= 2*d and k, d are positive integers.
    
    By detailed analysis, n that satisfies all conditions can only be:
    
    1) prime p when and only when p % 4 == 3
    2) 4 and 4*p, where p is ODD prime
    3) 16*p, where p is prime (2 included)
    """

    plist = iter(atkin_sieve(50000000))
    c = 1    # n = 4
    for p in plist:
        if p == 2:
            c += 1    # n = 32
        else:
            if p % 4 == 3:
                c += 1
            if 4 * p < 50000000:
                c += 1
            if 16 * p < 50000000:
                c += 1
    return c

def pe147(n=47, m=43):
    """
    It's possible to find a recursion formula from n*m to n*(m+1).
    """

    def nx1(n):
        """rectangles in n*1"""
        
        return n * (n+1) / 2 + n - 1

    def nxmp1(n, m):
        """extra rectangles when n*m -> n*(m+1)"""

        # extra horizontal and vertical rectangles
        c = nx1(n) + m * n * (n+1) / 2 - n + 1

        # extra diagonal rectangles
        for i in range(n):
            for x in xrange(1, 2*i+2):
                # diag involving 1 diag line cell
                c += max(min(2*(n-i)-1, 2*m-x+1), 0)

                # diag involving 1 diag base cell
                if i and x < 2*i+1:
                    value = max(min(2*(n-i), 2*m-x+2), 0)
                    if value > 0:
                        c += value
                    else:
                        break
        return c 
    
    c, cdict = 1, {(1, 1): 1}
    for j in range(2, n+1):
        cdict[(j, 1)] = nx1(j)
        c += (1 + (j <= m)) * cdict[(j, 1)]
        
        for i in range(2, min(j+1, m+1)):
            cdict[(j, i)] = cdict[(j, i-1)] + nxmp1(j, i-1)
            c += (1 + (j <= m) * (i != j)) * cdict[(j, i)]
    return c

def pe151():
    """
    It's possible to find a delicate data structure to model the situation.
    """

    from copy import deepcopy
    from fractions import Fraction
    
    SINGLE = {1000, 100, 10, 1}
    CHANGE = {0: 889, 1: 89, 2: 9, 3: 1}

    week_old = {(1111, 0): Fraction(1, 1)}
    for _ in range(13):
        week_new = {}
        for (envlope, n), p in week_old.iteritems():
            sheets = envlope/1000, (envlope%1000)/100, (envlope%100)/10, envlope%10
            s = sum(sheets)
            for i,a in enumerate(sheets):
                if a:
                    envnew = envlope - CHANGE[i]
                    m = n + (envnew in SINGLE)
                    if (envnew, m) in week_new:
                        week_new[(envnew, m)] += p * Fraction(a, s)
                    else:
                        week_new[(envnew, m)] = p * Fraction(a, s)
        week_old = deepcopy(week_new)

    expection = sum(p * n for (_, n), p in week_new.items())
    return float(expection)

def pe154(N=200000):
    """
    semi-brute-force, using Legendre Theorem to accelerate cycling
    """

    from . formulas import padic

    def ndivisor(n, d):
        x = 0
        while n:
            n /= d
            x += n
        return x
    
    seq2 = [ndivisor(n, 2) for n in xrange(N)]
    n2 = ndivisor(N, 2) - 12
    
    d5 = {}
    s5 = [sum(map(int, padic(i, 5, 'l'))) for i in xrange(N+1)]
    for i in range(1, 30):
        d5[i] = [n for n,s in enumerate(s5) if s==i]
    
    n5 = 12 * 4 + s5[N]
    c = 0
    for i in xrange(1, N/3+1):
        jk5 = n5 - s5[i]
        for dj in xrange(max(1, jk5-29), 30):
            for j in d5[dj]:
                if j < i:
                    continue
                if j > (N - i) / 2:
                    break
                
                k = N - i - j
                if s5[k] + dj < jk5:
                    continue
                
                if seq2[i] + seq2[j] + seq2[k] <= n2:
                    c += 3 if (i==j or j==k) else 6

        # if i % 1000 == 0:
        #     print i, c
    return c

def pe155(N=18):
    """
    brute-force
    """

    from fractions import Fraction

    def sep(n):
        for i in range(1, n/2+1):
            yield i, n-i
    
    fcircuits = [{Fraction(1, 1),}, {Fraction(1, 2), Fraction(2, 1)}]
    for n in range(3, N+1):
        new = set([])
        for i,j in sep(n):
            for f1 in fcircuits[i-1]:
                for f2 in fcircuits[j-1]:
                    new.add(f1+f2)
                    new.add(1 / (1 / f1 + 1 / f2))
        fcircuits.append(new)
        # print n
        
    c = fcircuits[-1]
    for i in range(N-1):
        c |= fcircuits[i]
    
    # answer: 3857447
    return c