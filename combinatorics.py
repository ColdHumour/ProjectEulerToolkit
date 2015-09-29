# -*- coding: utf-8 -*-

"""
combinatorics.py

Functions using to dealing with combinatorics problems.
Function list: 
    permutations, multisetPermutations
    C, combinations, combinations_with_replacement, limitedCombinations
    allPartitions, seqPartitions

@author: Jasper Wu
"""

from itertools import permutations, combinations, combinations_with_replacement

from . formulas import factorial


def C(n, k):
    if k == 0: return 1
    if k == 1: return n
    if k > n/2: k = n - k

    x = n
    for i in range(n-1, n-k, -1):
        x *= i
    for i in range(2, k+1):
        x /= i
    return x

def MP(amounts):
    """
    Calculate number of permutations of multiset, which defined by
    amounts = [amount1, amount2, ...]
    """

    s, p = 0, 1
    for v in amounts:
        s += v
        p *= factorial(v)
    return factorial(s) / p

def multisetPermutations(multiset):
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
        multiset = sorted(multiset, reverse=True)
        d = len(multiset) - 1
        
        E = [Node(n) for n in multiset]
        for i,n in enumerate(E[:-1]):
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

def limitedCombinations(choices):
    """
    Generate all combinations [x1, x2, ..., xn] which subjected to
    limited choices [[possible choices for xi] for i in 1..n]
    
    e.g. limitedCombinations([[1, 2], [3, 4]]) == [[1, 3], [1, 4], [2, 3], [2, 4]]
    """
    
    if len(choices) == 1:
        for x in choices[0]:
            yield [x]
    else:
        for x in choices[0]:
            for remains in limitedCombinations(choices[1:]):
                yield [x] + remains

def allPartitions(n, s, init=1):
    """
    make n partitions of s goods, output all possible partitions

    e.g allPartition(3, 5) == [[1, 1, 3], [1, 2, 2]]
    """

    if n == 1:
        yield [s]
    else:
        for i in range(init, s / n + 1):
            for result in allPartitions(n-1, s-i, i):
                yield [i] + result

def seqPartitions(sequence, p):
    """
    list all permutations of the sequence satisfying given partition p

    e.g seqPartition([1, 2, 3], [1, 2]) == [[[1], [2, 3]], [[2], [1, 3]], [[3], [1, 2]]]
    """

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
            output += [[list(subp)] + s for s in seqPartitions(newseq, p[1:])]
        return output