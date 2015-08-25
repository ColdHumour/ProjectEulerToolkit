# -*- coding: utf-8 -*-

"""
solutions(161-180).py

Some interesting solutions for problems 161-180 in Project Euler

@author: Jasper Wu
"""

import math

import ProjectEuler.combinatorics as pec
import ProjectEuler.formulas as pef


def pe161():
    """
    Dynamic programming. 
    Number 0-107 from upleft to bottomright row by row. Then a subset of (0~107) is a possible state. Define possible states at i as all subsets containing i. 
    For all possible states at (i-1), try to add 6 triominoes to them to cover the ith grid, then record the new state as possible states at i, and accumulate counters.
    """

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
            if not state & tm:
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

def pe162():
    """
    Easy combinatorics problem. Work out a close-form formula.
    Amount of strings with length n and 0, 1, A appears at least once: 
    F(n) = 15*16**(n-1) - 43*15**(n-1) + 41*14**(n-1) - 13**n, n >= 3
    """

    c = 0
    for n in range(3, 17):
        c += 15*16**(n-1) - 43*15**(n-1) + 41*14**(n-1) - 13**n

    slist = '0123456789ABCDEF'
    s = ''
    while c:
        c, q = divmod(c, 16)
        s = slist[q] + s

    # answer: 3D58725572C62302
    return s

def pe163(n=36):
    """
    Divide and conquer for all kinds of triangles.
    Can summarize out some close-form formula. See reference articles.
    """

    def backdiff3_120(n):
        return 3 * ( 1 + [0,1,1,1][n % 4] + 2*[0,1][n % 2] )

    def backdiff3_90(n):
        return 6 * ( 1 + [1,2][n % 2] + [0,1,1][n % 3] + [0,2,1,1,2][n % 5] )

    def backdiff3_60(n):
        return 1 + [0,1][n % 2] + 6 * [0,1,1][n % 3] + 2 * (n>1)

    d0 = d1 = d2 = 0
    for i in range(n):
        d3 = backdiff3_120(i) + backdiff3_60(i) + backdiff3_90(i)
        d2 = d2 + d3
        d1 = d1 + d2
        d0 = d0 + d1
    
    # answer: 343047
    return d0

def pe164(limit=20):
    """
    Dynamic programming. 
    Define possible states at i as all left two digits when the number is i-digit. 
    For all possible states at (i-1), try to expand to i-digit number while satisfying constraints, then record the new state as possible states at i, and accumulate counters.
    """

    situation = {}
    for i in range(10):
        for j in range(i, 10-i):
            situation[(i, j)] = situation[(j, i)] = 1

    for _ in range(3, limit):
        temp = {}
        for (i, j), n in situation.iteritems():
            for k in range(10-i-j):
                if (j, k) in temp:
                    temp[(j, k)] += n
                else:
                    temp[(j, k)] = n
        situation = temp

    c = 0
    for (i, j), n in situation.iteritems():  
        if i + j < 9: 
            for k in range(1, 10-i-j):
                c += n

    # answer: 378158756814587
    return c

def pe165():
    """
    Brute-force 
    """

    def intersection(line1, line2):
        x1, y1, x2, y2 = line1
        x3, y3, x4, y4 = line2

        a, b = x2 - x1, y2 - y1
        c, d = x3 - x1, y3 - y1
        e, f = x4 - x1, y4 - y1

        t1 = c * b - a * d
        t2 = e * b - a * f

        if t1 * t2 >= 0:
            return 0

        g, h = x4 - x2, y4 - y2
        i, j = x4 - x3, y4 - y3

        t1 = e * j - i * f
        t2 = g * j - i * h

        if t1 * t2 >= 0:
            return 0
        
        denom = a * j - b * i
        ynum = b * c * j + a * j * y1 - b * i * y3
        xnum = a * (c * j - d * i) + x1 * denom
        g = pef.ggcd((denom, ynum, xnum))
        return (ynum/g, xnum/g, denom/g)
    
    lines, t, c = [], [], set()
    s = 290797
    for _ in xrange(20000):
        s = (s*s) % 50515093
        t.append(s % 500)
        if len(t) == 4:
            for l in lines:
                c.add(intersection(t, l))
            lines.append(tuple(t))
            t = []
    
    # answer: 2868868
    return len(c) - 1

def pe166():
    """
    Using Gauss-Jordan eliminate algorithm to simplify to 7 variables.
    Enumerating all possible 1e7 permutations, checking possible magic sum, and counting.
    """
    
    # v and (9-n for n in v) are symmetric.
    # so look up to 4999999 is enough.
    def vecgen(n):
        v = [0] * (n-1) + [1]
        while v < [5] + [0] * (n-1):
            yield v

            i = 1
            while v[-i] == 9:
                i += 1
            v = v[:-i] + [v[-i]+1] + [0]*(i-1)

    def check(v, (a, b, c), nlist):
        # using the reduced matrix and magic sum to 
        # check the validation of vector
        
        # VARMAT = [[-1.,  0.,  0., -1., -1., -1., -1.],
        #           [-1.,  1., -1., -1.,  0., -1., -2.],
        #           [ 1., -1.,  1.,  1.,  1.,  2.,  2.],
        #           [ 1.,  0.,  0.,  1.,  0.,  0.,  1.],
        #           [ 1., -1., -1.,  0.,  0.,  0.,  0.],
        #           [ 1.,  0.,  1.,  1.,  1.,  1.,  2.],
        #           [-1.,  1.,  0., -1., -1., -1., -2.],
        #           [ 0.,  1.,  1.,  1.,  0.,  0.,  0.],
        #           [ 0.,  0.,  0.,  0.,  1.,  1.,  1.]]
        # RHSCOEF = [ 1.,  1., -2., -1., 0., -2.,  1., -1., -1.] 
        
        count = 0
        for n in nlist:
            if not 0 <= n - a <= 9:
                continue
            if not 0 <= n - b <= 9:
                continue
            if not 0 <= n - c <= 9:
                continue
            if not 0 <= v[0] + v[3] + b - n <= 9:
                continue
            if not 0 <= b + c - v[1] - n <= 9:
                continue
            if not 0 <= 2*n - b - c - v[2] <= 9:
                continue
            if not 0 <= 2*n - b - c + v[1] - v[2] - v[5] <= 9:
                continue
            if not 0 <= c - v[1] + v[2] + v[5] + v[6] - n <= 9:
                continue
            count += 1
        return count
    
    # precalculate restricted magic sums by vmin and vmax
    magicsum = {}
    for a in range(37):
        bmax = min(a+10, 37)
        for b in range(a, bmax):
            magicsum[(a, b)] = range(b, bmax)

    # there's only 1 sol for all-0s and all-9s, count from 2
    c = 2
    for v in vecgen(7):
        if not 0 <= v[1] + v[2] - v[0] <= 9:
            continue
        
        vsum = (v[1] + v[2] + v[3], 
                v[4] + v[5] + v[6], 
                v[0] + v[3] + v[6])
        vmin, vmax = min(vsum), max(vsum)
        if vmax - vmin > 9:
            continue

        c += 2 * check(v, vsum, magicsum[(vmin, vmax)])
    
    # answer: 7130034
    return c

def pe167():
    """
    See reference articles.
    """
    
    def descUlam(v):
        """
        To get useful information for 1-addtive (Ulam) sequence U(2, v).
        Output: uncycled part, cycling differences
        """

        a = (v+7) / 2
        b = 2 * v + 2

        # initialize Ulam sequence till 2*v+5
        s = [2, v, v+2]
        while s[-1] < 2*v - 3:
            s.append(s[-1] + 2)
        s += [2*v - 1, 2*v + 1, 2*v + 2, 2*v + 3, 2*v + 5]

        head, diff, cycle, c = s[:a], [2], [], 0

        # just keep a reasonable length to get next U(2, v)
        s = [n for n in s if n%2]
        while 1:
            if s[0] + b == s[-1] + 2:
                n = s[1] + b
                d = n - s[-1]
                s[:2] = []
            else:
                n = s[-1] + 2
                d = 2
            s.append(n)

            # try to find shortest cycling part
            diff.append(d)
            if len(diff) > v:
                while c and d != diff[c]:
                    cycle[:1] = []
                    c -= 1

                if d == diff[c]:
                    cycle.append(d)
                    c += 1

                if c * 2 == len(diff):
                    return head, cycle

    def nUlam(v, n):
        """
        Using periodicity to find nth U(2, v).
        """
        
        head, cycle = descUlam(v)
        d1, d2, csum = len(head), len(cycle), sum(cycle)
        nc = (n - d1) / d2
        nr = n - d1 - nc * d2
        return 2 * v + 3 + csum * nc + sum(cycle[:nr-1])

    c = 0
    for v in range(5, 22, 2):
        c += nUlam(v, 10**11)
        
    # answer: 3916160068885
    return c

def pe168():
    """
    easy brute-force with a, c, k, n satisfying:
    c*10**(n-1) + a*10**(n-2) + b = k * (a*10**(n-1) + b*10 + c)
    where b must be integer less and equal than n-2 digit.
    """
    
    csum = sum(range(11, 100, 11))
    for n in range(3, 101):
        for k in range(1, 10):
            for a in range(1, 10):
                for c in range(a, 10):
                    s = 10 * k - 1
                    t = (10**(n-1) - k) * c - 10**(n-2) * s * a
                    if t % s == 0:
                        b = t / s
                        if len(str(b)) <= n - 2:
                            csum += (a*10**(n-1) + b*10 + c) % 100000
    
    # answer: 59206
    return csum % 100000

def pe169(N=10**25):
    """
    After represent N in binary, then we can get a list of powers of 2
    [pn, p(n-1), ... , p1], where N = 2**an + ... + 2**a1
    Then for [pn, p(n-1), ... , p1, 0], final amount of expression
    E(N) = a0 + b0, where a, b can be calculated recursivly:
    an = 0, bn = 1
    a(i-1) = ai + bi, b(i-1) = ai * (pi-p(i-1)-1) + bi * (pi-p(i-1))
    """
    
    s = bin(N)[2:]
    m = len(s) - 1
    p2 = [m-i for i,n in enumerate(s) if n=='1'] + [0]

    a, b = 0, 1
    for i,n in enumerate(p2[1:]):
        k = p2[i] - n
        a, b = a + b, a * (k - 1) + b * k
    
    # answer: 178653872807
    return a + b

def pe170():
    """
    Properties of the first integer
    1) contains at most 2 digits
    2) must be less than 49
    3) 2, 5, 10*k, 11*k can be eliminated by hand
    4) 9 gives out a naive lower bound 9768352140, which can be used to eliminate 7 and 8
    """

    res = '1234567890'
    
    # test 9 as the first number
    for n in pec.permutations('2345678', 7):
        n = int(''.join(n))
        c = str(n * 9)
        if len(c) == 8 and set(c) == set('12345678'):
            res = max(res, '9{}0'.format(c))
    
    # 9 * 0 = 0
    # 9 * 1 = 9
    # 9 * 8537246 = 76835214
    # This gives out a lower bound 9768352140
    
    # 20+ is enough, so I'm lazy to test 10+.
    # 21, 23 do not need to test
    for n in range(24, 30):
        sp = [s for s in '0456789' if s != str(n%10)]
        for m in pec.permutations(sp):
            for i in range(5):
                a = int('1' + ''.join(m[:i]))
                b = int('3' + ''.join(m[i:]))
                c = str(a * n) + str(b * n)
                if len(c) <= 10:
                    if set(c) == set('0123456789'):
                        res = max(res, c, str(b * n) + str(a * n))

                a = int('1' + ''.join(m[i:]))
                b = int('3' + ''.join(m[:i]))
                c = str(a * n) + str(b * n)
                if len(c) <= 10:
                    if set(c) == set('0123456789'):
                        res = max(res, c, str(b * n) + str(a * n))
    
    # 31, 32, 41, 42 do not need to test
    for n in range(34, 40) + range(43, 50):
        if n == 44:
            continue

        t = '4' if n < 40 else '3'
        sp = [s for s in '056789'+t if s != str(n%10)]
        for m in pec.permutations(sp):
            for i in range(5):
                a = int('1' + ''.join(m[:i]))
                b = int('2' + ''.join(m[i:]))
                c = str(a * n) + str(b * n)
                if len(c) <= 10:
                    if set(c) == set('0123456789'):
                        res = max(res, c, str(b * n) + str(a * n))

                a = int('1' + ''.join(m[i:]))
                b = int('2' + ''.join(m[:i]))
                c = str(a * n) + str(b * n)
                if len(c) <= 10:
                    if set(c) == set('0123456789'):
                        res = max(res, c, str(b * n) + str(a * n))

    # 27 * 36508 = 985716
    # 27 * 149   = 4023 
    # answer: 9857164023
    return res

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