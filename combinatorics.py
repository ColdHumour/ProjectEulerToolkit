# -*- coding: utf-8 -*-

"""
combinatorics.py

Functions using to dealing with combinatorics problems.
Function list: 
    permutations, combinations, combinations_with_replacement
    allPartitions, seqPartitions

@author: Jasper Wu
"""

from itertools import permutations, combinations, combinations_with_replacement


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

def allPartitions(n, s, init=1):
    """
    make n partitions of s goods, output all possible partitions

    e.g partition(3, 5) == [[1, 1, 3], [1, 2, 2]]
    """

    if n == 1:
        return [[s]]
    else:
        output = []
        for i in range(init, s / n + 1):
            output += [[i] + result for result in allPartitions(n-1, s-i, i)]
        return output

def seqPartitions(sequence, p):
    """
    list all permutations of the sequence satisfying given partition p

    e.g partition(3, 5) == [[1, 1, 3], [1, 2, 2]]
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