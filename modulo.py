# -*- coding: utf-8 -*-

"""
modulo.py

Modulo arithmetic functions implementing formulas via fast algorithms.
Function list:
    add_mod, cprod, mul_mod, pow_mod,
    fac, fac_mod, inv_mod,
    sum_over_mod, sum_power_series_mod,
    tabulate_inv_mod,
    tabulate_fac_mod, tabulate_fac_inv,
    tabulate_bernoulli_mod,
    faulhaber_mod_coefs,

@author: Jasper Wu
"""

try:
    from gmpy2 import isqrt as _isqrt, fac as _fac, powmod as _powmod
    fac = lambda x: int(_fac(int(x)))
    pow_mod = lambda x, y, m: int(_powmod(int(x), int(y), int(m)))
except:
    fac = None
    pow_mod = pow

try:
    from . ext.c_formula_int64 import c_sum_over_mod_int64
    sum_over_mod = c_sum_over_mod_int64
except:
    sum_over_mod = None


def add_mod(MOD, *args):
    """return sum(args) % MOD"""

    s = 0
    for n in args:
        s = (s + n) % MOD
    return s


def cprod(seq):
    """return seq[0] * seq[1] * ... * seq[-1]"""

    output = 1
    for i in iter(seq):
        output *= i
    return output


def mul_mod(MOD, *args):
    """return cprod(args) % MOD"""

    s = 1
    for n in args:
        s = (s * n) % MOD
    return s


def _factorial(n):
    """return n!"""

    if n < 0:
        raise ValueError("n in n! must be positive!")
    if n == 0 or n == 1:
        return 1

    output = 1
    for i in range(2, n+1):
        output *= i
    return output

if fac is None:
    fac = _factorial


def fac_mod(n, m):
    """return n! % m"""

    if n < 0:
        raise ValueError("n in n! must be positive!")
    if n == 0 or n == 1:
        return 1

    output = 1
    for i in range(2, n+1):
        output *= i
        output %= m
    return output


def inv_mod(n, m):
    """return n^(-1) mod m using Extended Euclid Algorithm"""

    n %= m
    if n == 0 or m <= 0:
        return 0
    
    m0, x0, x1 = m, 0, 1
    while n != 0:
        x0, x1 = x1, x0 - m // n * x1
        m, n = n, m % n

    if m == 1:
        return x0 % m0
    else:
        return 0


def _sum_over_mod(n):
    """return n % 2 + n % 3 + ... + n % (n-1)"""

    from itertools import takewhile, count

    sm = i = 0
    for i in takewhile(lambda x: n//x - n//(x+1) > 4, count(1)):
        a = n % (n//(i+1) + 1)
        b = n % (n//i) if i > 1 else 1
        c = (a-b) // i + 1
        sm += b*c + i*(c - 1)*c // 2
    sm += sum(n % j for j in range(2, n//(i+1) + 1))
    return sm

if sum_over_mod is None:
    sum_over_mod = _sum_over_mod


def sum_power_series_mod(i, n, m):
    """sum of x^i mod m from i=1 to n, for i = 0, 1, 2, 3"""

    if i == 0:
        return n
    elif i == 1:
        n %= 2 * m
        return ((n*(n+1)) >> 1) % m
    elif i == 2:
        n %= 6 * m
        m3 = 3 * m
        res = ((n*(n+1)) >> 1) % m3
        return ((res * ((2*n+1) % m3)) % m3) // 3
    elif i == 3:
        n %= 2 * m
        res = ((n*(n+1)) >> 1) % m
        return (res * res) % m
    else:
        return 0


def tabulate_inv_mod(n, MOD):
    """tabulate i^(-1) % MOD for 0 <= i <= n"""

    invs = [1] * (n+1)
    for i in range(2, n+1):
        invs[i] = inv_mod(i, MOD)
    return invs


def tabulate_fac_mod(n, MOD):
    """tabulate i! % MOD for 0 <= i <= n"""

    facs = [1] * (n+1)
    for i in range(2, n+1):
        facs[i] = facs[i-1] * i % MOD
    return facs


def tabulate_fac_inv(facs, MOD):
    """tabulate i!^(-1) % MOD for 0 <= i <= n"""

    facinvs = [1] * len(facs)
    for i in range(2, len(facs)):
        facinvs[i] = inv_mod(facs[i], MOD)
    return facinvs


def tabulate_bernoulli_mod(n, MOD, invs):
    """tabulate Bernoulli numbers % MOD for 0 <= i <= n"""

    bernoulli = [0] * (n+1)
    bernoulli[0] = 1
    bernoulli[1] = MOD - invs[2]
    for i in range(2, n+1, 2):
        B = 0
        for k in range(1, i+1):
            b = 0
            c = k
            for r in range(1, k+1):
                if r & 1:
                    b -= c * pow(r, i, MOD) % MOD
                    if b < 0:
                        b += MOD
                else:
                    b += c * pow(r, i, MOD) % MOD
                    if b >= MOD:
                        b -= MOD

                if r < k:
                    c = mul_mod(MOD, c, k-r, invs[r+1])

            B += b * invs[k+1] % MOD
            if B >= MOD:
                B -= MOD
        bernoulli[i] = B
    return bernoulli


def faulhaber_mod_coefs(k, MOD, invs, bernoulli):
    """
    generate coefficients of faulhaber's formula for 1^k + ... + n^k
    invs is modular inverse from 0 to at least k+1
    bernoulli is bernoulli numbers modulo MOD from 0 to k+1
    """

    if k == 1:
        return [0, invs[2], invs[2]]

    faulhaber = [0] * (k+2)
    faulhaber[k] = MOD - bernoulli[1]

    c = k * invs[2] % MOD if k & 1 else 1
    for i in range(1+(k&1), k, 2):
        faulhaber[i] = c * bernoulli[k-i+1] % MOD
        c = mul_mod(MOD, c, (k+1-i)*(k-i), invs[i+1], invs[i+2])
    faulhaber[i+2] = c

    return faulhaber
