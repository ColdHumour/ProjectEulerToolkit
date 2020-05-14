# -*- coding: utf-8 -*-

"""
algolibrary.py

Interesting algorithms, but not optimal
Function list: 
    atkin_sieve
    is_coprime
    power_mod
    
@author: Jasper Wu
"""


try:
    from gmpy2 import sqrt, gcd
except:
    from math import sqrt
    from fractions import gcd


def atkin_sieve(limit=1000000):
    """
    Sieve of Atkin, which is a much stronger version than sieve of Eratosthenes
    https://en.wikipedia.org/wiki/Sieve_of_Atkin
    """

    plist = [0] * limit

    # n = 3x^2 + y^2 section
    x = 3
    for i in range(0, 12*int(sqrt((limit-1)/3)), 24):
        x += i
        y_limit = int(12*sqrt(limit-x)-36)
        n = x + 16
        for j in range(-12, y_limit+1, 72):
            n += j
            plist[n] = not plist[n]

        n = x + 4
        for j in range( 12, y_limit+1, 72):
            n += j
            plist[n] = not plist[n]

    # n = 4x^2 + y^2 section
    x = 0
    for i in range(4, 8*int(sqrt((limit-1)/4))+4, 8):
        x += i
        n = x + 1
        if x % 3:
            for j in range(0, 4*int(sqrt(limit-x))-3, 8):
                n += j
                plist[n] = not plist[n]
        else:
            y_limit = 12 * int(sqrt(limit-x)) - 36
            
            n = x + 25
            for j in range(-24, y_limit+1, 72):
                n += j
                plist[n] = not plist[n]

            n = x + 1
            for j in range( 24, y_limit+1, 72):
                n += j
                plist[n] = not plist[n]

    # n = 3x^2 - y^2 section
    x = 1
    for i in range(3, int(sqrt(limit/2))+1, 2):
        x += 4 * i - 4
        n = 3 * x
        if n > limit:
            y = (int(sqrt(n-limit)) >> 2) << 2
            n -= y * y
            s = 4 * y + 4
        else:
            s = 4

        for j in range(s, 4*i, 8):
            n -= j
            if n <= limit and n % 12 == 11:
                plist[n] = not plist[n]

    x = 0
    for i in range(2, int(sqrt(limit/2))+1, 2):
        x += 4 * i - 4
        n = 3 * x

        if n > limit:
            y = ((int(sqrt(n-limit)) >> 2) << 2) - 1
            n -= y * y
            s = 4 * y + 4
        else:
            n -= 1
            s = 0

        for j in range(s, 4*i, 8):
            n -= j
            if n <= limit and n % 12 == 11:
                plist[n] = not plist[n]

    # eliminate squares        
    for n in range(5, int(sqrt(limit))+1, 2):
        if plist[n]:
            for k in range(n*n, limit, n*n):
                plist[k] = False

    return [2,3] + filter(plist.__getitem__, range(5,limit,2))


def is_coprime(a, b):
    """return whether a and b are coprime"""

    return gcd(a, b) == 1


def power_mod(a, b, n):
    """return (a ** b) % n"""
    
    r = a % n
    if r in (0, 1):
        return r
    
    r = 1
    while b:
        if b % 2:
            r = (r * a) % n
        b /= 2
        a = (a * a) % n
    return r % n
