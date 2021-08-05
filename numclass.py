# -*- coding: utf-8 -*-

"""
numclass.py

Classes representing special numbers.

Function list:
- IntQuadIrr(D)

@author: Jasper Wu
"""

class IntQuadIrrBase:
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
    return type("IntQuadIrr", (IntQuadIrrBase, ), {"D": D, "d": np.sqrt(D)})