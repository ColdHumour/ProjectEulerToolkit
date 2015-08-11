# -*- coding: utf-8 -*-

"""
solutions(0-160).py

Some interesting solutions for problems 0-160 in Project Euler

@author: Jasper Wu
"""

from math import sqrt
import ProjectEuler.combinatorics as pec
import ProjectEuler.prime as pep


def pe128(index):
    """
    Using patterns found in hexNeighbors()
    """
    
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

    x = i = 2
    while 1:
        if pep.isprime(6*i-1):
            if pep.isprime(6*i+1) and pep.isprime(12*i+5):
                x += 1
                if x == index:
                    return 3*i*(i-1)+2
            if pep.isprime(6*i+5) and pep.isprime(12*i-7):
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

    plist = iter(pep.atkin_sieve(50000000))
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

    # answer: 0.464399
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

    # answer: 479742450
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
    return len(c)

def pe156():
    """
    (1) there's only 1641 turning points in range [0, 100000], end points included
    (2) the maximum and the minimum point are always two end points of range [b*100000, (b+1)*100000]
    """

    G = [0] + [d*10**(d-1) for d in range(1, 12)]

    def f(n, b, minus_n=True):
        s = map(int, str(n))
        d, count = len(s), 0
        for i,a in enumerate(s):
            tb = 10**(d-i-1)
            if i + 1 < d:
                count += a * G[d-i-1] + (a > b) * tb + (a == b) * (n % tb + 1)
            else:
                count += (a >= b)
        return count - n * minus_n

    def nbSearch(d, n, detail=True):
        sols = [n]
        for i in range(1, 10):
            xl = f(n-i, d)
            if xl == 0:
                sols.append(n-i)
            else:
                l = n - i
                break
        for i in range(1, 10):
            xr = f(n+i, d)
            if xr == 0:
                sols.append(n+i)
            else:
                r = n + i
                break

        if detail:
            return sols, (l, r, xl, xr)
        else:
            return sols

    def biSearch(d, a, b):
        sols = []
        
        xa, xb = f(a, d), f(b, d)
        if xa == 0:
            nbsols, (_, a, _, xa) = nbSearch(d, a)
            sols += nbsols
        if xb == 0:
            nbsols, (b, _, xb, _) = nbSearch(d, b)
            sols += nbsols

        if xa * xb > 0:
            return sols

        while b - a > 1:
            c = (a + b) / 2
            xc = f(c, d)
            if xc == 0:
                sols += nbSearch(d, c, False)
                return sols

            xa, xb = f(a, d), f(b, d)
            if xa < 0 and xb > 0:
                if xc < 0:
                    a, b = c, b
                else:
                    a, b = a, c
            else:
                if xc < 0:
                    a, b = a, c
                else:
                    a, b = c, b
        return []

    def sumSolutions(d, batch = 100000):
        sd = str(d)
        gb = G[len(str(batch)) - 1] - batch + 1

        def gplus(i, fi):
            # compute [f(i*batch, d) - i*batch] using f(i, d)
            x = str(i).count(sd)
            return i * gb + (fi - x) * batch + x - i 

        sols, tps, prg = set([]), [0], []
        trend, count, xn1, xb1 = 0, 0, 0, 0

        for i in range(1, batch):
            # sequentially compute f(i, d)
            count += str(i).count(sd)

            # record all solutions within batch
            xn2 = count - i
            if xn2 == 0:
                sols.add(i)

            # record all turning points
            if xn1 < xn2:
                if trend == -1:
                    tps.append(i-1)
                trend = 1
            elif xn1 > xn2:
                if trend == 1:
                    tps.append(i-1)
                trend = -1
            xn1 = xn2

            # record all possible ranges
            xb2 = gplus(i, count)
            if i > 1 and xb1 * xb2 <= 0:
                prg.append(((i-1)*batch, i*batch))
            xb1 = xb2
        tps += [batch]

        for a, b in prg:
            sols = sols.union(biSearch(d, a, b))

        # get sols in expand search region
        if d > 1:
            for i in range(batch, d*batch): 
                count += str(i).count(sd)
                xb2 = gplus(i, count)
                if i > 1 and xb1 * xb2 <= 0:
                    sols = sols.union(biSearch(d, (i-1)*batch, i*batch))
                xb1 = xb2

        return sum(list(sols))
    
    res = 0
    for d in range(1, 10):
        x = sumSolutions(d)
        print "Sum of all solutions of f(n, {0}) = n: {1}".format(d, x)
        res += x
    print "\nFinal sum:", res, '\n'

    # Sum of all solutions of f(n, 1) = n: 22786974071
    # Sum of all solutions of f(n, 2) = n: 73737982962
    # Sum of all solutions of f(n, 3) = n: 372647999625
    # Sum of all solutions of f(n, 4) = n: 741999999540
    # Sum of all solutions of f(n, 5) = n: 100000000000
    # Sum of all solutions of f(n, 6) = n: 2434703999430
    # Sum of all solutions of f(n, 7) = n: 1876917059570
    # Sum of all solutions of f(n, 8) = n: 15312327487352
    # Sum of all solutions of f(n, 9) = n: 360000000000

    # Final sum: 21295121502550 

    return res

def pe157():
    """
    It can be shown that a and b can be only written as:
        a = 2**x1 * 5**y1 * g
        b = 2**x2 * 5**y2 * g
    where g is a positive number not contain 2 or 5 as divisors.
    Then we can discuss different relationships between x1, x2, y1, y2:
    (1) x1 <= x2 & y1 <= y2
    (2) (x1 - x2) * (y1 - y2) < 0 & a < b
    and finally get the answer.
    """

    c, gdict1, gdict2 = 0, {}, {}
    for n in range(1, 10):
        cn = 0
        for s in range(n+1):
            for t in range(n+1):
                # case1: x1 <= x2 & y1 <= y2
                # need to find all divisors of (2**s * 5**t + 1)
                # it may have 2 or 5 in divisors
                if (s, t) in gdict1:
                    g, p2, p5 = gdict1[(s, t)]
                else:
                    x = 2**s * 5**t + 1
                    p2 = 0
                    while x%2 == 0:
                        x /= 2
                        p2 += 1
                    p5 = 0
                    while x%5 == 0:
                        x /= 5
                        p5 += 1
                    g = len(pep.divisorDecomp(x))
                    gdict1[(s, t)] = (g, p2, p5)
                
                g, p2, p5 = gdict1[(s, t)]
                cn += g * (n - s + p2 + 1) * (n - t + p5 + 1)

                # case2: (x1 - x2) * (y1 - y2) < 0
                # need to find all divisors of (2**s + 5**t)
                # it won't have 2 or 5 in divisors
                # it only happens when s >= 1 and t >= 1
                if s and t:
                    if (s, t) in gdict2:
                        g = gdict2[(s, t)]
                    else:
                        g = len(pep.divisorDecomp(2**s + 5**t))
                        gdict2[(s, t)] = g
                    
                    cn += g * (n - s + 1) * (n - t + 1)
        print 'Number of sols when n = {0}: {1}'.format(n, cn)
        c += cn
    print '\nFinal sum: {}\n'.format(c)
    
    # Number of sols when n = 1: 20
    # Number of sols when n = 2: 102
    # Number of sols when n = 3: 356
    # Number of sols when n = 4: 958
    # Number of sols when n = 5: 2192
    # Number of sols when n = 6: 4456
    # Number of sols when n = 7: 8260
    # Number of sols when n = 8: 14088
    # Number of sols when n = 9: 23058

    # Final sum: 53490
    
    return c

def pe158():
    """
    Proven close form p(n) = C(26, n) * (2**n - n - 1)
    So easy. 
    """
    
    return max(pec.C(26, n) * (2**n - n - 1) for n in range(1, 27))