# -*- coding: utf-8 -*-

"""
Project Euler Toolkit

Dependencies:
    shutil
    fractions
    itertools
    collections
    cython, Cython

    numpy, scipy: scientific computation
    sympy: symbolic computation
    gmpy2: number theory
    ortools: programming and optimization

@author: Jasper Wu
"""


# dependency test
# default: numpy, scipy, sympy
DEPENDENCIES = ['gmpy2']  # 'cython', 'ortools'

for pkg in DEPENDENCIES:
    try:
        __import__(pkg)
    except:
        print("WARNING: Fail to import {}!".format(pkg))


# cython extensions
from . import ext


# useful functions in current package
from . combinatoric import (
    C, C_mod, MP,

    multiset_permutations,
    limited_combinations,

    all_subsets,
    all_partitions,
    seq_partitions,

    composite_perm,
    inverse_perm,
    rank_perm,
    unrank_perm,
)

from . formula import (
    sqrt, is_square, isqrt, iroot,
    gcd, ggcd, extended_gcd, lcm, llcm,
    sum_floor, legendre_symbol,
    padic, max_subarray,

    pythag_triple_tree,
    co_prime_tree,
    stern_brocot_tree,

    rational_continous_frac,
    irrational_continous_frac,
    continous_frac_convergent,

    best_rational_approx,
    find_closest_lattice_point_to_line,
)

from . linalg import (
    mat_pow_mod,
    gauss_jordan_elimination,
)

from . modulo import (
    add_mod, cprod, mul_mod, pow_mod,
    fac, fac_mod, inv_mod,
    sum_over_mod, sum_power_series_mod,
    tabulate_inv_mod,
    tabulate_fac_mod, tabulate_fac_inv,
    tabulate_bernoulli_mod,
    faulhaber_mod_coefs,
)

from . prime import (
    primes_list,
    is_prime,
    prime_divisor_decomp,
    all_divisors,
    euler_phi,
    mobius,
    mobius_list,
    factor_sieve,
    primepi,
)

from . utils import (
    timepast,
    memoize,
    clear_cython_cache,
    # find_solution,
)

from . import equation
from . import linalg
from . import polynomial


__all__ = [
    'combinatoric',
    'equation',
    'formula',
    'generator',
    'linalg',
    'polynomial',
    'prime',
    'utils',
]
