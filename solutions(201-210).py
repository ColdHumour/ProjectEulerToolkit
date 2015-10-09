# -*- coding: utf-8 -*-

"""
solutions(201-210).py

Some interesting solutions for problems 201-210 in Project Euler

@author: Jasper Wu
"""

import math

import ProjectEuler.combinatorics as pec
import ProjectEuler.formulas as pef
import ProjectEuler.prime as pep


def pe201():
    """
    states[i] := all sums using i numbers among nlist[:i]
    For states[i], if m is unique, then states[i][m] = 1, else >1.
    """
    
    nlist = [n*n for n in range(1, 101)]
    d = len(nlist) / 2
    states = [{0: 1}] + [{} for _ in range(d)]

    for i,n in enumerate(nlist[:-1]):
        for k in range(min(d-1, i), max(0, i-d)-1, -1):
            for m0, v in states[k].iteritems():
                m1 = m0 + n
                if v >= 2:
                    states[k+1][m1] = 2
                else:    
                    states[k+1][m1] = states[k+1].get(m1, 0) + 1

    n = nlist[-1]
    for m0,v in states[d-1].iteritems():
        m1 = m0 + n
        if v >= 2:
            states[d][m1] = 2
        else:    
            states[d][m1] = states[d].get(m1, 0) + 1

    # answer: 115039000
    return sum(k for k,v in states[-1].iteritems() if v==1)

def pe202():
    """
     | \| \| \| \| \| \|   Expand triangle mirrors infinitely, then every 
     +--+--+--+--+--+--+-  possible path is a C-C line from bottomleft to
    C|\B|\A|\C|\B|\A|\C|\  somewhere upright and doesn't cross any vertex.
     | \| \| \| \| \| \|   
     +--+--+--+--+--+--+-  From the table can easity see that C-C lines are
    B|\A|\C|\B|\A|\C|\B|\  symmetric and lie on y = x + 3k. Without loss of
     | \| \| \| \| \| \|   generality, we just need to count all (a, b) with
     +--+--+--+--+--+--+-  a > b.
    A|\C|\B|\A|\C|\B|\A|\  
     | \| \| \| \| \| \|   Also, from (0, 0) to (a, b), the amount of cross
     +--+--+--+--+--+--+-  points with horizontal, vertical and diagonal
    C  B  A  C  B  A  C    lines is (a-1) + (b-1) + (a+b-1) = 2a + 2b - 3.
    
    Therefore to bounce N times, we need to find all (a, b) that:
        (1) a + b = (N + 3) / 2
        (2) gcd(a, b) = 1
        (3) a - b = 3k
    
    Recording to reference articles PE202, the number of solutions of (1), (2)
    and (3) is:
    
        sum m(d) * f(n/d) for all d satisfying d|n
    
    where m(x) if mobius function, f(x) = floor(x/2) - floor(x/3)
    """
    
    def f(x):
        return x / 2 - x / 3
    
    N = 12017639147
    N = (N + 3) / 2 # 6008819575
    
    # 6008819575 = 5*5 * 11 * 17 * 23 * 29 * 41 * 47
    pfactors = (5, 11, 17, 23, 29, 41, 47)
    count = f(N)
    for i in range(1, len(pfactors)+1):
        for pcomb in pec.combinations(pfactors, i):
            d = pef.cprod(pcomb)
            if i % 2:
                count -= f(N / d)
            else:
                count += f(N / d)
    
    # answer: 1209002624
    return count * 2

def pe203(N=51):
    """
    Use prime divisors to decide whether a combination number is square-free.
    """
    
    def oprfac(n1, n2, direction=1):
        for n,a in n2.iteritems():
            if n in n1:
                n1[n] += direction * a
            else:
                n1[n] = a
        return n1

    NFACS = {n: dict(pep.primeDivisorDecomp(n)) for n in range(2, N)}

    def squarefreeC(n):
        c, f, sfc = n, NFACS[n], set()
        if all(v==1 for v in f.values()):
            sfc.add(c)

        for i in range(2, n/2+1):
            c = c * (n+1-i) / i
            f = oprfac(f, NFACS[n+1-i], 1)
            f = oprfac(f, NFACS[i], -1)
            if all(v==1 for v in f.values()):
                sfc.add(c)
        return sfc

    distinct_sf = {1, 2}
    for i in xrange(3, N):
        distinct_sf |= squarefreeC(i)

    # answer: 34029210557338
    return sum(list(distinct_sf))

def pe204(N=10**9):
    """
    Recursive brute-froce.
    """

    p100 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 
            43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]

    def f(n, plist):
        
        if len(plist) == 1:
            return int(math.log(N / n, plist[0]))
        else:
            c = 0
            for i,p in enumerate(plist):
                q = p
                while q * n <= N:
                    c += 1 + f(q * n, plist[i+1:])
                    q *= p
            return c
    
    # answer: 2944730
    return f(1, p100) + 1

def pe205():
    def ndicedist(d, n, iscum=True):
        """
        probability (or cumulative) distribution of n d-sided dices.
        """

        dice = range(1, d+1)
        state = {a:1 for a in dice}
        for _ in range(n-1):
            new_stat = {}
            for a in state.keys():
                for b in dice:
                    if a+b in new_stat:
                        new_stat[a+b] += state[a]
                    else:
                        new_stat[a+b] = state[a]
            state = new_stat

        total = float(sum(state.values()))
        for m in xrange(min(state), max(state)+1):
            state[m] = state[m] / total
            if iscum:
                state[m] += state.get(m-1, 0)
        return state

    # peter beats colin
    peter = ndicedist(4, 9, False)
    colin = ndicedist(6, 6, True)
    p = sum(colin.get(i-1, 1) * peter.get(i, 0) \
            for i in xrange(min(colin)+1, max(peter)+1))

    # answer: 0.5731441
    return round(p, 7)

def pe206():
    """
    Brute-force with little optimization
    """
    
    check = '123456789'
    choices = [list('3210')] + \
              [list(check[::-1] + '0') for _ in range(6)] + \
              [['7', '3']]

    for p in pec.limitedCombinations(choices):
        n = int('1' + ''.join(p))
        n = str(n*n)
        if n[::2] == check:

            # answer: 1389019170
            # 1389019170 ** 2 = 1929374254627488900

            return int(math.sqrt(int(n))) * 10

def pe207(N=12345):
    """
    4^t = 2^t + k  =>  2^t = (1 + sqrt(1+4k)) / 2
    Let k = n(n+1), then 2^t = n + 1, beginning with n = 1.
    
    Therefore:
    (1) when 2^x - 1 <= n <= 2^(x+1) - 2, there are n partitions
        and x perfect partitions.
    (2) when n(n+1) <= m <= (n+1)(n+2) - 1, P(m) = x / n, where x,n
        satisfying (1).
    (3) from (1) and (2), when 2^x - 1 <= n <= 2^(x+1) - 2, P(m)
        monotonically decreases from x / (2^x-1) to x / (2^(x+1)-2).
    (4) x / (2^(x+1)-2) < (x+1) / (2^(x+1)-1)
    
    The properties above is enough to decide a search algorithm:
        find x that 2^x - 1 <= x * N <= 2^(x+1) - 2
        find n that n-1 <= x * N < n, i.e. x * N + 1
        return n * (n+1)
    """

    i = 1
    while 2**(i+1) - 2 <= i * N:
        i += 1
    n = int(i * N) + 1
    
    # answer: 44043947822
    return n*(n+1)

def pe208(N=70):
    """
    Record Cartesian coordinates for every step and merge same states.
    The most difficult part is to deal with precision.
    Can be optimized from various ways (float caching for instance).
    There are lots of valuable things in the discussion thread.
    Definitely worth to review again and again.
    """
    
    # always facing to the direaction whose angle between y-axis is 2pi/5 * k
    # and directions of next steps are (alpha + pi/5) and (alpha - pi/5)
    VECMAP = {0: [(4, -0.588,  0.809), (1,  0.588,  0.809)],
              1: [(0,  0.588,  0.809), (2,  0.951, -0.309)],
              2: [(1,  0.951, -0.309), (3,  0.000, -1.000)],
              3: [(2,  0.000, -1.000), (4, -0.951, -0.309)],
              4: [(3, -0.951, -0.309), (0, -0.588,  0.809)]}
    
    # final move must end with facing north
    FINOPT = {1: (0.588,  0.809), 4: (-0.588,  0.809)}

    ARCLIB = {(4, -0.588,  0.809): 1, (1,  0.588,  0.809): 1}
    for i in range(2, N):
        tmplib = {}
        for (d,x,y),n in ARCLIB.iteritems():
            for dd, vx, vy in VECMAP[d]:
                state = (dd, round(x+vx, 3), round(y+vy, 3))        
                if state in tmplib:
                    tmplib[state] += n
                else:
                    tmplib[state] = n
        ARCLIB = tmplib

    c = 0
    for (d,x,y),n in ARCLIB.iteritems():
        if d in (1, 4):
            vx, vy = FINOPT[d]
            if abs(x+vx) + abs(y+vy) < 0.01:
                c += n
    
    # answer: 331951449665644800
    return c

def pe209():
    """
    Take input as integer in 0~63. Then it's equivalent to how many binary
    functions F(x) with 0 <= x <= 63 subjected to F(a) * F(b) = 0 for some
    a, b which b is determined by a.
    There are some cycles in the mapping. So we just need to find valid
    ways of each cycle, which becomes a easy combinatoric problem.
    """

    def valid_permutations(n):
        """
        Count valid permutations of chain with length n. Each element of chain 
        can only be 0 or 1. No consecutive 1s. Head and tail is regarded as
        consecutive too.
        Therefore it would be 0~floor(n/2) possible 1s.
        For a chain of A 0s and B 1s, there's C(A-1, B-1) chains with head of 1
        and C(A, B) chains with head of 0.
        """
        
        x = 1
        for i in range(1, n/2+1):
            x += pec.C(n-i-1, i-1) + pec.C(n-i, i)
        return x

    pairs = {}
    for n in range(64):
        s = bin(n)[2:]
        s = '0' * (6-len(s)) + s
        a, b, c = map(int, s[:3])
        s = s[1:] + str(a ^ (b & c))
        pairs[n] = int(s, 2)

    # There are 6 cycles in pairs, started with 0, 1, 3, 5, 9, 21
    starts, lchains = (0, 1, 3, 5, 9, 21), []
    for i in starts:
        z = set([i])
        n = pairs[i]
        while n > i:
            z.add(n)
            n = pairs[n]
        lchains.append(len(z))
    
    # answer: 15964587728784
    return pef.cprod(map(valid_permutations, lchains))