# -*- coding: utf-8 -*-

"""
datastructure.py

Classes representing different data structure.

Function list:
- IntQuadIrr(D)
- DisjointSet(data)
- FenwickTree(data, MOD)

@author: Jasper Wu
"""

import numpy as np


class _IntQuadIrrBase:
    """class of a + b * sqrt(D)"""

    D = d = None

    def __init__(self, a, b=0):
        self.a = a
        self.b = b
        self.hash = hash((a, b))

    def eval(self):
        return self.a + self.b * self.d

    def __neg__(self):
        return self.__class__(-self.a, -self.b)
    
    def __add__(self, other):
        if isinstance(other, int):
            return self.__class__(self.a+other, self.b)
        elif isinstance(other, self.__class__):
            return self.__class__(self.a+other.a, self.b+other.b)
        else:
            raise NotImplementedError

    def __iadd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, int):
            return self.__class__(self.a-other, self.b)
        elif isinstance(other, self.__class__):
            return self.__class__(self.a-other.a, self.b-other.b)
        else:
            raise NotImplementedError

    def __rsub__(self, other):
        return -self.__sub__(other)

    def __mul__(self, other):
        if isinstance(other, int):
            return self.__class__(other*self.a, other*self.b)
        elif isinstance(other, self.__class__):
            return self.__class__(self.a*other.a + self.b*other.b*self.D, self.a*other.b + self.b*other.a)
        else:
            raise NotImplementedError

    def __rmul__(self, other):
        return self.__mul__(other)

    def __eq__(self, other):
        if isinstance(other, int):
            return self.b == 0 and self.a == other
        elif isinstance(other, self.__class__):
            return self.b == other.b and self.a == other.a
        else:
            raise NotImplementedError

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        if isinstance(other, int):
            other = self.__class__(other)

        if self.a <= other.a and self.b <= other.b and not (self.a == other.a and self.b == other.b):
            return True
        elif self.a <= other.a and self.b >= other.b:
            return (other.a - self.a)**2 > (self.b - other.b)**2 * self.D
        elif self.a >= other.a and self.b <= other.b:
            return (other.a - self.a)**2 < (self.b - other.b)**2 * self.D
        else:
            return False

    def __le__(self, other):
        return self.__eq__(other) or self.__lt__(other)

    def __gt__(self, other):
        return not self.__eq__(other) and not self.__lt__(other)

    def __ge__(self, other):
        return not self.__lt__(other)
    
    def __repr__(self):
        return "{} + {}√{}".format(self.a, self.b, self.D)
    
    def __hash__(self):
        return self.hash

def IntQuadIrr(D):
    """class of a + b * sqrt(D)"""

    return type("IntQuadIrr", (_IntQuadIrrBase, ), {"D": D, "d": np.sqrt(D)})


class DisjointSet:
    """disjoin set or union set"""

    def __init__(self, data):
        self.data = data
        self.parent = {k: k for k in data}
        self.group = {k: {k, } for k in data}

    def find(self, k):
        if self.parent[k] != k:
            self.parent[k] = self.find(self.parent[k])
        return self.parent[k]

    def union(self, a, b):
        x = self.find(a)
        y = self.find(b)
        if x != y:
            if len(self.group[x]) > len(self.group[y]):
                self.parent[y] = x
                self.group[x] |= self.group[y]
                del self.group[y]
            else:
                self.parent[x] = y
                self.group[y] |= self.group[x]
                del self.group[x]


class FenwickTree:
    def __init__(self, data, MOD=0):
        self.data = data
        self.N = len(data)
        self.M = MOD

    def add(self, i, val):
        """sum from index 1 to i based on Fenwick tree"""
        while i <= self.N: 
            self.data[i] += val
            if self.M and self.data[i] >= self.M:
                self.data[i] -= self.M
            i += i & -i
        return 1

    def sum_range(self, i):
        s = 0
        while i:
            s += self.data[i]
            if self.M and s >= self.M:
                s -= self.M
            i -= i & -i
        return s
