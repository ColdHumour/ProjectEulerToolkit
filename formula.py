# -*- coding: utf-8 -*-

"""
formula.py

Functions implementing formulas via fast algorithms.
Function list:
    sqrt, is_square, isqrt, iroot,
    gcd, ggcd, extended_gcd, lcm, llcm,
    sum_floor, generate_integer_quotients,
    legendre_symbol,
    padic, max_subarray,
    mex,

    pythag_triple_tree
    co_prime_tree
    stern_brocot_tree

    rational_continous_frac
    irrational_continous_frac
    continous_frac_convergent

    best_rational_approx
    best_rational_approx_for_log
    find_closest_lattice_point_to_line

@author: Jasper Wu
"""

from math import gcd, sqrt
from collections import deque

try:
    from gmpy2 import is_square, iroot
    isqrt = lambda x: int(_isqrt(int(x)))
except:
    is_square = None
    isqrt = None
    iroot = None


# Supplementry Implementations
def _is_square(n):
    """return whether n is a perfect square"""

    s = int(sqrt(n))
    return s * s == n

if is_square is None:
    is_square = _is_square


def _isqrt(n):
    """return integer square root of n"""

    return int(n**0.5)

if isqrt is None:
    isqrt = _isqrt


def _iroot(n, m):
    """return integer m-th root of n, and whether n is a perfect power"""

    r = int(n**(1./m))
    check = (r + 1)**m
    if check == n:
        return r + 1, True
    elif check < n:
        return r + 1, False
    else:
        return r, r**m == n

if iroot is None:
    iroot = _iroot


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
        if g == 1:
            return g
        for n in seq[:-2]:
            g = gcd(g, n)
            if g == 1:
                return 1
        return g


def extended_gcd(a, b):
    """
    return gcd(a, b), x, y that a*x + b*y = gcd(a, b)
    using Extended Euclid Algorithm, where a, b > 0
    """

    assert a >= 0 and b >= 0

    x = v = 0
    y = u = 1
    while a:
        q = b // a
        r = b - q * a
        m = x - u * q
        n = y - v * q

        b = a
        a = r
        x = u
        y = v
        u = m
        v = n

    gcd = b
    return gcd, x, y


def lcm(a, b):
    """return the least common multiple of a and b"""

    return a // gcd(a, b) * b


def llcm(seq):
    """
    return the greatest common divisor (gcd) for n integers,
    where n can larger than 2
    """

    if len(seq) < 2:
        raise ValueError("There should be at least 2 integers!")
    elif len(seq) == 2:
        return lcm(seq[0], seq[1])
    else:
        l = lcm(seq[-2], seq[-1])
        for n in seq[:-2]:
            if l % n:
                l = lcm(l, n)
        return l


def padic(n, p, ntype='s'):
    """change integer n from base 10 to base p"""

    base = '0123456789' + ''.join(map(chr, range(65, 92)))
    snp = ''
    while True:
        if n == 0:
            break
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


def mex(gvalues):
    """return mex(g1, g2, ...), gvalues is list or set, set is better"""

    n = 0
    while n in gvalues:
        n += 1
    return n


def sum_floor(n, xmin, xmax):
    """sum up n//x from x = xmin to xmax"""

    nrt = isqrt(n)
    res = 0
    if xmin > n:
        return 0
    if xmax > n:
        xmax = n

    if xmax <= nrt:
        for x in range(xmin, xmax+1):
            res += n // x
    elif xmin >= nrt:
        real_xmin = n // xmax
        real_xmax = n // xmin
        a0 = 0
        a1 = n // real_xmin
        for x in range(real_xmin, real_xmax+1):
            a0 = a1
            a1 = n // (x+1)
            ub = a0 if a0 < xmax else xmax
            lb = a1 if a1 >= xmin-1 else xmin-1
            res += (ub - lb) * x
    else:
        real_xmin = n // xmax
        if real_xmin > xmin:
            real_xmin = xmin

        a0 = 0
        a1 = n // real_xmin
        for x in range(real_xmin, nrt+1):
            a0 = a1
            a1 = n // (x+1)

            if x >= xmin:
                res += a0

            if a1 < xmax:
                ub = a0 if a0 < xmax else xmax
                res += (ub - a1) * x

        if x == n // x:
            res -= x
    return res


def generate_integer_quotients(n):
    """return list of all n//x, sorted descendingly"""

    iqs = []
    for m in range(1, isqrt(n)+1):
        iqs.append(n // m)

    if m != n // m:
        iqs.append(m)

    iqs += list(range(m-1, 0, -1))
    return iqs


def legendre_symbol(a, p):
    """
    Legendre symbol
    Define if a is a quadratic residue of prime p
    Details see: http://en.wikipedia.org/wiki/Legendre_symbol
    """

    if p == 2 or a % p == 0:
        return 1

    ls = pow(a % p, (p - 1) >> 1, p)
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
        return ((a-2*b+2*c,   2*a-b+2*c,  2*a-2*b+3*c),
                (a+2*b+2*c,   2*a+b+2*c,  2*a+2*b+3*c),
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
def rational_continous_frac(p, q=1):
    """
    return continued fraction expansion of rational number p / q
    """

    p, q, cfrac = int(p), int(q), []
    while p % q:
        cfrac.append(p//q)
        p, q = q, p % q
    if q:
        cfrac.append(p)
    return tuple(cfrac)


def irrational_continous_frac(d, p=0, q=1, limit=100):
    """
    return continued fraction expansion of quadratic irrational number (p + sqrt(d)) / q
    about quadratic irrational number: https://en.wikipedia.org/wiki/Quadratic_irrational
    """

    from fractions import Fraction

    d, p, q = int(d), int(p), int(q)

    sd = sqrt(d)
    if int(sd) * int(sd) == d:
        raise ValueError("D is perfect square!")

    a = int((p + sd) / q)
    repetend, pairspq = [a], []
    for _ in range(limit):
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

    from fractions import Fraction

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


# Approximation
def best_rational_approx(D, N):
    """
    get best lower and upper rational approximation of sqrt(D) and denominator no larger than N
    use features of Stern-Brocot Tree and Farey Sequence
    return (a, b, c, d), where a/b < sqrt{D} < c/d 
    """

    a, b, c, d = 0, 1, 1, 0
    while a + c <= N:
        if (b + d) * (b + d) * D > (a + c) * (a + c):
            a += c
            b += d
        else:
            c += a
            d += b
    return (a, b, c, d)


def best_rational_approx_for_log(a, b, n=10):
    """
    generate first n best rational approximations for log_b(a)
    according to the semiconvergents of its continued fraction expression
    """

    a = float(a)
    b = float(b)

    p0 = 1
    q0 = 0
    i = 0
    while a >= b:
        a /= b
        i += 1
    a, b = b, a
    p1 = i
    q1 = 1
    output = [(p1, q1)]

    while b > 1 and len(output) < n:
        i = 0
        while a >= b:
            a /= b
            i += 1
            output.append((p1*i+p0, q1*i+q0))

        a, b = b, a
        p0, q0 = p1, q1
        p1, q1 = output[-1]

    return output


def find_closest_lattice_point_to_line(A, B, C, x_min, x_max):
    """
    find the closest lattice point (x, y) that closest to Ax + By + C = 0
    which is equivalent to minimize |Ax + By + C|
    return x
    """    
    
    if x_min >= x_max or A == 0:
        return x_min
    if B == 0:
        x = -(-(-2*C // A) // 2)
        if x < xmin:
            x = xmin
        elif x > xmax:
            x = xmax
        return x
    if A < 0:
        A, C = -A, -C
    if B < 0:
        B = -B
    if A >= B:
        A %= B
        if A == 0:
            return x_min

    y_min = -(A * (2*x_max + 1) + 2*C) // (2*B) + 1
    y_max = -(A * (2*x_min - 1) + 2*C) // (2*B)
    y = find_closest_lattice_point_to_line(B, A, C, y_min, y_max)
    x = -(-(-(2*B*y + 2*C) // A) // 2)
    if x < x_min:
        x = x_min
    elif x > x_max:
        x = xmax
    
    d = abs(A*x + B*y + C)

    dd = abs(A*x_min + B*(y_max+1) + C)
    if dd < d:
        dd = d
        x = x_min

    dd = abs(A*x_max + B*(y_min-1) + C)
    if dd < d:
        x = x_max
    return x
