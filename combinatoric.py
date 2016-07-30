# -*- coding: utf-8 -*-

"""
combinatoric.py

Functions using to dealing with combinatorics problems.
Function list:
    C, C_mod, MP
    multiset_permutations, limited_combinations
    all_subsets, all_partitions, seq_partitions

@author: Jasper Wu
"""

from . prime import primes_list, euler_phi
from . formula import fac_mod


def C(n, k):
    if k > n:
        return 0
    if k > n/2:
        k = n - k
    if k == 0:
        return 1
    if k == 1:
        return n

    output = n
    for i in range(n-1, n-k, -1):
        output *= i
    for i in range(2, k+1):
        output //= i
    return output


def C_mod(n, k, m):
    """Return C(n, k) % m"""

    if k > n:
        return 0
    if k > n//2:
        k = n - k
    if k == 0:
        return 1
    if k == 1:
        return n % m

    output = fac_mod(n, m)

    x = (fac_mod(k, m) * fac_mod(n - k, m)) % m
    if x:
        t = euler_phi(m)
        output *= pow(x, t-1, m)
        output %= m
    else:
        plist = primes_list(n+1)
        output = 1
        for p in plist:
            x = 0
            q = p
            while q <= n:
                x += n // q
                q *= p

            if p <= k:
                q = p
                while q <= k:
                    x -= k // q
                    q *= p

            if p <= n - k:
                q = p
                while q <= n - k:
                    x -= (n - k) // q
                    q *= p

            output *= pow(p, x, m)
            output %= m
    return output


def MP(amounts):
    """
    Calculate number of permutations of multiset, which defined by
    amounts = [amount1, amount2, ...]
    """

    from ProjectEulerToolkit.formula import fac

    s, p = 0, 1
    for v in amounts:
        s += v
        p *= fac(v)
    return fac(s) // p


def multiset_permutations(multiset):
    """
    List out all permutations of a multiset [a, a, ..., b, b, ..., c, c, ...]
    """

    class Node(object):
        def __init__(self, v, to=None):
            self.val = v
            self.to = to

    def visit(n):
        vlist = [n.val]
        while n.to:
            n = n.to
            vlist.append(n.val)
        return vlist

    if len(multiset) == 1:
        yield multiset
    elif len(multiset) == 2:
        if multiset[0] == multiset[1]:
            yield multiset
        else:
            yield multiset
            yield multiset[::-1]
    else:
        E = [Node(n) for n in sorted(multiset, reverse=True)]
        for i, n in enumerate(E[:-1]):
            n.to = E[i+1]

        head, i, afteri = E[0], E[-2], E[-1]
        yield visit(head)

        while afteri.to or afteri.val < head.val:
            if afteri.to and i.val >= afteri.to.val:
                beforek = afteri
            else:
                beforek = i
            k = beforek.to
            beforek.to = k.to
            k.to = head
            if k.val < head.val:
                i = k
            head, afteri = k, i.to
            yield visit(head)


def limited_combinations(choices):
    """
    Generate all combinations [x1, x2, ..., xn] which subjected to
    limited choices [[possible choices for xi] for i in 1..n]

    e.g. limited_combinations([[1, 2], [3, 4]]) == [[1, 3], [1, 4], [2, 3], [2, 4]]
    """

    if len(choices) == 1:
        for x in choices[0]:
            yield [x]
    else:
        for x in choices[0]:
            for remains in limited_combinations(choices[1:]):
                yield [x] + remains


def all_subsets(fullset, xmin=1, xmax=None):
    """
    return all subsets of the fullset, minimum and maximum set size can be specified
    xmin: minimum subset size
    xmax: maximum subset size

    e.g. all_subsets([1, 2, 3], 1, None) = [[1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]]
    """

    from itertools import combinations

    if len(fullset) < xmin:
        raise ValueError("Minimum subset size too large!")

    if xmax is None:
        xmax = len(fullset)

    for i in range(xmin, xmax+1):
        for subset in combinations(fullset, i):
            yield subset


def all_partitions(n, s, xmin=1, xmax=None):
    """
    make n partitions of s goods, return all possible partitions
    xmin: minimum partition size
    xmax: maximum partition size

    e.g. all_partition(3, 5) == [[1, 1, 3], [1, 2, 2]]
    """

    if xmax is None:
        if n == 1:
            yield [s]
        else:
            for i in range(xmin, s // n + 1):
                for result in all_partitions(n-1, s-i, i, xmax):
                    yield [i] + result
    else:
        if s > n * xmax:
            yield None
        elif n == 1:
            yield [s]
        else:
            for i in range(max(xmin, s-(n-1)*xmax), min(s//n, xmax)+1):
                for result in all_partitions(n-1, s-i, i, xmax):
                    if result is not None:
                        yield [i] + result


def seq_partitions(sequence, p):
    """
    list all permutations of the sequence satisfying given partition p

    e.g. seq_partition([1, 2, 3], [1, 2]) == [[[1], [2, 3]], [[2], [1, 3]], [[3], [1, 2]]]
    """

    from itertools import combinations

    if len(sequence) != sum(p):
        raise ValueError("The length of sequence doesn't match given partition!")

    if len(p) == 1:
        output = []
        for subp in combinations(sequence, p[0]):
            output += [[list(subp)]]
        return output
    else:
        output = []
        for subp in combinations(sequence, p[0]):
            newseq = [ele for ele in sequence if ele not in subp]
            output += [[list(subp)] + s for s in seq_partitions(newseq, p[1:])]
        return output
