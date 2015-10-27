# -*- coding: utf-8 -*-

"""
formula.py

Functions implementing formulas via fast algorithms.
Function list:
    sqrt, is_square, gcd, ggcd 
    factorial, cprod
    sum_mod, pow_mod, iter_associate
    legendre_symbol
    padic, max_subarray

    pythag_triple_tree
    co_prime_tree
    stern_brocot_tree

    rational_continous_frac, 
    irrational_continous_frac, 
    continous_frac_convergent
    
@author: Jasper Wu
"""


from fractions import Fraction
from collections import deque

try:
    from gmpy2 import is_square, factorial, sqrt, gcd
    from gmpy2 import powmod as pow_mod
except:
    from math import sqrt
    from fractions import gcd
    pow_mod = pow
    is_square = _is_square
    factorial = _factorial

try:
    from ProjectEulerToolkit.ext._formula import _c_sum_mod
    sum_mod = _c_sum_mod
except:
    sum_mod = _sum_mod


# Supplementry Implementations 
def _is_square(n):
    """return whether n is a perfect square"""

    s = int(sqrt(n))
    return s * s == n

def _factorial(n):
    """return n!"""

    if n < 0:
        raise ValueError("n in n! must be positive!")
    if n==0 or n==1: 
        return 1

    output = 1
    for i in xrange(2, n+1):
        output *= i
    return output


# Useful Functions
def cprod(seq):
    """return seq[0] * seq[1] * ... * seq[-1]"""

    output = 1
    for i in iter(seq):
        output *= i
    return output

def ggcd(seq):
    """
    return the greatest common divisor (gcd) for n integers, 
    where n can larger than 2
    """

    if len(seq) < 2:
        raise ValueError("There should be at least 2 integers!")
    elif len(seq) == 2:
        return gcd(seq[0], seq[1])
    else:
        g = gcd(seq[-2], seq[-1])
        if g == 1: return g
        for n in seq[:-2]:
            g = gcd(g, n)
            if g == 1: return 1
        return g

def padic(n, p, ntype='s'):
    """change integer n from base 10 to base p"""

    base = '0123456789' + ''.join(map(chr, range(65, 92)))
    snp = ''
    while True:
        if n == 0: break
        n, r = divmod(n, p)
        snp += base[r]

    snp = snp[::-1]
    if ntype == 's':
        return snp
    elif ntype == 'n':
        return int(snp)
    elif ntype == 'l':
        return list(snp)

def max_subarray(array):
    """return max sum of any continous subarray of an array"""

    max_so_far = max_ending_here = 0
    for x in array:
        max_ending_here = max(0, max_ending_here + x)
        max_so_far = max(max_so_far, max_ending_here)
    return max_so_far


# Modulo Functions
def _sum_mod(n):
    """return n%2 + n%3 + ... + n%(n-1)"""
    
    from itertools import takewhile, count

    sm = i = 0
    for i in takewhile(lambda x: n/x - n/(x+1) > 4, count(1)):
        a = n % (n/(i+1) + 1)
        b = n % (n/i) if i > 1 else 1
        c = (a-b) / i + 1
        sm += b*c + i*(c - 1)*c / 2
    sm += sum(n%j for j in xrange(2, n/(i+1) + 1))
    return sm

def power_mod(a, b, n):
    """
    return (a ** b) % n
    """
    
    r = a % n
    if r in (0, 1):
        return r
    
    a, b, n, r = int(a), int(b), int(n), 1
    while b:
        if b % 2:
            r = (r * a) % n
        b /= 2
        a = (a * a) % n
    return r % n

def iter_associate(f, x, n):
    """
    f is a bivariable function following associate law, namely f(a, f(b, c)) = f(f(a, b), c)
    return doing f iteratively for n times with identical input x, namely f(f(...f(x, x)..., x), x)
    """

    n, r = n - 1, x
    while n:
        if n % 2:
            r = f(x, r)
        n /= 2
        x = f(x, x)
    return r

def legendre_symbol(a, p):
    """
    Legendre symbol
    Define if a is a quadratic residue of odd prime p
    Details see: http://en.wikipedia.org/wiki/Legendre_symbol
    """

    ls = pow(a, (p - 1)/2, p)
    if ls == p - 1:
        return -1
    return ls


# Tree Generation
def pythag_triple_tree(triple=(3, 4, 5), forward=True, trust=False):
    """
    Primitive Pythagorean Triple (PPT) (a, b, c) is integer triple satisfying 
        a**2 + b**2 = c**2
        gcd(a, b) = gcd(b, c) = gcd(a, c) = 1
    Given a PPT, it can generate three more PPTs. And if recursively applying this branching function from (3, 4, 5), all PPTs in form (odd, even, odd) can be generated in a trinary tree, which covers the entire set of PPTs completely and uniquely.
    When forward is True, return all three children of current PPT.
    When forward is False, return its parent in the PPT tree.
    """

    a, b, c = triple
    if not trust:
        if a**2 + b**2 != c**2:
            raise ValueError("Invalid Primitive Pythagorean Triple")
        if gcd(a, b) * gcd(a, c) * gcd(b, c) != 1:
            raise ValueError("Invalid Primitive Pythagorean Triple")

    if forward:
        return (( a-2*b+2*c,  2*a-b+2*c,  2*a-2*b+3*c),
                ( a+2*b+2*c,  2*a+b+2*c,  2*a+2*b+3*c),
                (-a+2*b+2*c, -2*a+b+2*c, -2*a+2*b+3*c))
    else:
        if triple == (3, 4, 5):
            return triple
        else:
            return (abs(-a-2*b+2*c), abs(-2*a-b+2*c), -2*a-2*b+3*c)

def co_prime_tree(pair=(0, 0), trust=False):
    """
    All co-prime pairs can be generated from (2, 1) (for (odd, even) and (even, odd) pairs) and (3, 1) (for (odd, odd) pairs). 
    It follows a trinary tree, from co-prime pair (a, b), we get: (2*a - b, a), (2*a + b, a), (a + 2*b, b)
    It can be shown that the co-prime pairs in the tree are disjoint complete.
    """

    if pair == (0, 0):
        return ((2, 1), (3, 1))

    a, b = pair
    if not trust:
        if gcd(a, b) != 1:
            raise ValueError("Invalid co-prime pair!")            
    return ((2*a - b, a), (2*a + b, a), (a + 2*b, b))

def stern_brocot_tree():
    """
    Stern-Brocot Tree, an infinite complete binary tree in which the vertices correspond one-for-one to the positive rational numbers, whose values are ordered from left to right as in a search tree. It related to Farey series closely.
    https://en.wikipedia.org/wiki/Stern%E2%80%93Brocot_tree
    """

    sbt = deque([1, 1])
    while True:
        sbt += [sbt[0] + sbt[1], sbt[1]]
        sbt += [sbt[1] + sbt[2], sbt[2]]
        yield (sbt.popleft(), sbt.popleft())


# Continuous Fraction Functions
def rational_continous_frac(p, q=1, limit=100):
    """
    return continued fraction expansion of rational number p / q
    """

    p, q, cfrac = int(p), int(q), []
    while p%q:
        print p, q
        cfrac.append(p/q)
        p, q = q, p%q
    if q:
        cfrac.append(p)
    return tuple(cfrac)

def irrational_continous_frac(d, p=0, q=1, limit=100):
    """
    return continued fraction expansion of quadratic irrational number (p + sqrt(d)) / q
    about quadratic irrational number: https://en.wikipedia.org/wiki/Quadratic_irrational
    """
    
    d, p, q = int(d), int(p), int(q)
    
    sd = sqrt(d)
    if int(sd) * int(sd) == d:
        raise ValueError("D is perfect square!")
    
    a = int((p + sd) / q)
    repetend, pairspq = [a], []
    for _ in xrange(limit):
        p = Fraction(a*q - p, 1)
        q = Fraction(d - p*p, q)
        a = int((p + sd) / q)
        
        if (p, q) in pairspq:
            i = pairspq.index((p, q))
            return tuple(repetend[:i+1] + [tuple(repetend[i+1:])])
        
        pairspq.append((p, q))
        repetend.append(a)
    raise ValueError("Repetend is longer than {0:d}, please try higher limit!".format(limit))

def continous_frac_convergent(cfrac):
    """
    given continued fraction, return series of convergents 
    """

    if len(cfrac) < 2:
        raise ValueError("Continued fraction must longer than 2!")
    
    a0, a1 = cfrac[:2]
    p0, p1, q0, q1 = a0, a0*a1+1, 1, a1
    cvg = [Fraction(p0, q0), Fraction(p1, q1)]
    for a in cfrac[2:]:
        p0, p1 = p1, a*p1 + p0
        q0, q1 = q1, a*q1 + q0
        cvg.append(Fraction(p1, q1))
    return cvg