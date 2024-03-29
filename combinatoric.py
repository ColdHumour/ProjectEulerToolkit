# -*- coding: utf-8 -*-

"""
combinatoric.py

Functions using to dealing with combinatorics problems.
Function list:
    C, C_mod, MP
    multiset_permutations, limited_combinations
    all_subsets, all_partitions, seq_partitions
    composite_perm, inverse_perm, rank_perm, unrank_perm
    cycle_index_mod_p, merge_cycle_index_mod_p

@author: Jasper Wu
"""

import numpy as np

from . prime import primes_list, euler_phi, all_divisors
from . formula import gcd
from . modulo import inv_mod, fac_mod, tabulate_fac_mod, tabulate_fac_inv


def C(n, k):
    if n < 0 or k < 0 or k > n:
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

    if n < 0 or k < 0 or k > n:
        return 0
    if k > n//2:
        k = n - k
    if k == 0:
        return 1
    if k == 1:
        return n % m

    output = fac_mod(n, m)

    x = (fac_mod(k, m) * fac_mod(n - k, m)) % m
    if x and gcd(x, m) == 1:
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


def composite_perm(perm1, perm2, *other_perms):
    """composite two or more perms"""

    perm = perm2[perm1]
    for other_perm in other_perms:
        perm = other_perm[perm]
    return perm.copy()


def inverse_perm(perm):
    """
    find inverse of a permuation in bijection form, return np.array
    e.g. inverse([1, 2, 3, 0]) = [3, 0, 1, 2]
    """

    n = len(perm)
    inv = np.zeros(n, dtype=np.short)
    for i in range(n):
        inv[perm[i]] = i
    return inv


def rank_perm(perm):
    """map perm to a unique integer between 0 and n!-1"""

    inv = inverse_perm(perm)
    perm = perm.copy()

    slist = []
    for i in range(len(perm)-1, 0, -1):
        s = perm[i]
        perm[i], perm[inv[i]] = perm[inv[i]], perm[i]
        inv[s], inv[i] = inv[i], inv[s]
        slist.append(s)

    r = 0
    for i, s in enumerate(slist[::-1]):
        r = int(s) + (i+2)*r
    return r


def unrank_perm(r, n):
    """resolve proper perm from 0 to n-1 with rank r"""

    perm = np.arange(n, dtype=np.short)
    for i in range(n-1, -1, -1):
        r, q = divmod(r, i+1)
        perm[i], perm[q] = perm[q], perm[i]
    return perm


def cycle_index_mod_p(gtype, n, MOD):
    """
    get cycle index of permutation/cyclic/dihedral groups, with parameters mod p
    return [(param, ((cycle length, exponent), ...)), ...]

    gtype (S/C/D): permutation/cyclic/dihedral group
    n (int): order
    MOD (int): prime number to modulo
    """

    if gtype == "S":
        facs = tabulate_fac_mod(n, MOD)
        facinv = tabulate_fac_inv(facs, MOD)
        def binom(n, i):
            return (facs[n] * facinv[i] % MOD) * facinv[n-i] % MOD

        states = [[[], n, 1]]  # [term, remaining length, coeff]
        cycle_index = []  # (term, coeff)
        for i in range(n, 0, -1):
            states_new = states[:]
            for term, l, c in states:
                l_new = l
                c0_new = c
                e = 0
                while l_new >= i:
                    e += 1
                    term_new = term + [(i, e)]
                    c0_new = c0_new * binom(l_new, i) % MOD
                    l_new -= i
                    c_new = (c0_new * facinv[e] % MOD) * pow(facs[i-1], e, MOD) % MOD
                    if l_new:
                        states_new.append([term_new, l_new, c_new])
                    else:
                        cycle_index.append((tuple(sorted(term_new)), c_new * facinv[n] % MOD))
            states = states_new
    elif gtype == "C":
        cycle_index = []  # (term, coeff)
        ninv = inv_mod(n, MOD)
        for d in all_divisors(n):
            cycle_index.append((((d, n // d),), euler_phi(d) * ninv % MOD))
    elif gtype == "D":
        if n <= 2:
            return cycle_index_mod_p("C", n, MOD)

        inv2 = inv_mod(2, MOD)
        cycle_index = {}
        for term, coeff in cycle_index_mod_p("C", n, MOD):
            cycle_index[term] = coeff * inv2 % MOD

        if n & 1:
            cycle_index[((1, 1), (2, (n-1)//2))] = inv2
        else:
            cycle_index[((1, 2), (2, (n-2)//2))] = inv2 * inv2 % MOD
            cycle_index[((2, n//2),)] += inv2 * inv2 % MOD
            if cycle_index[((2, n//2),)] >= MOD:
                cycle_index[((2, n//2),)] -= MOD
        cycle_index = [(term, coeff) for term, coeff in cycle_index.items()]

    cycle_index.sort(key=lambda x: sum(xx[1] for xx in x[0]), reverse=True)
    return cycle_index


def merge_cycle_index_mod_p(cycle_index_1, cycle_index_2, MOD):
    """
    merge cycle index of Z(G1) and Z(G2) to get Z(G1 * G2), with parameters mod p
    return [(param, ((cycle length, exponent), ...)), ...]
    """

    cycle_index_3 = {}
    for term1, c1 in cycle_index_1:
        for term2, c2 in cycle_index_2:
            term = {}
            for l1, e1 in term1:
                for l2, e2 in term2:
                    g = gcd(l1, l2)
                    l = l1 * l2 // g
                    e = e1 * e2 * g
                    if l in term:
                        term[l] += e
                    else:
                        term[l] = e

            term = tuple((l, term[l]) for l in sorted(term))
            c = c1 * c2 % MOD
            if term in cycle_index_3:
                cycle_index_3[term] += c
                if cycle_index_3[term] >= MOD:
                    cycle_index_3[term] -= MOD
            else:
                cycle_index_3[term] = c

    cycle_index_3 = [(term, coeff) for term, coeff in cycle_index_3.items()]
    cycle_index_3.sort(key=lambda x: sum(xx[1] for xx in x[0]), reverse=True)
    return cycle_index_3
