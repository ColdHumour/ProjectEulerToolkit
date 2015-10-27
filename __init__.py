# -*- coding: utf-8 -*-

"""
Project Euler Toolkit

Dependencies:
    shutil
    fractions
    itertools
    collections
    cython, Cython
    
    numpy: arrays and matrices
    pulp: linear and non-linear optimization
    gmpy2: integers, rationals, etc.

@author: Jasper Wu
"""


# dependency test
DEPENDENCIES = ['cython', 'numpy', 'pulp', 'gmpy2']

for pkg in DEPENDENCIES:
    try:
        __import__(pkg)
    except:
        print("WARNING: Fail to import {}!".format(pkg))


# cython extensions
import ProjectEulerToolkit.ext


# useful functions in dependent packages
from itertools import permutations, combinations, combinations_with_replacement


# useful functions in current package
from ProjectEulerToolkit.combinatoric import (
    C, MP,
    
    multiset_permutations, 
    limited_combinations,
    
    all_partitions, 
    seq_partitions,
)

from ProjectEulerToolkit.formula import (
    sqrt, is_square, gcd, ggcd,
    factorial, cprod,
    sum_mod, pow_mod, legendre_symbol,
    padic, max_subarray,

    pythag_triple_tree,
    co_prime_tree,
    stern_brocot_tree,

    rational_continous_frac, 
    irrational_continous_frac, 
    continous_frac_convergent,
)

from ProjectEulerToolkit.prime import (
    primes_list,
    is_prime,
    is_coprime,
    prime_divisor_decomp,
    all_divisors,
)


import ProjectEulerToolkit.combinatoric
import ProjectEulerToolkit.equation
import ProjectEulerToolkit.formula
import ProjectEulerToolkit.generator
import ProjectEulerToolkit.prime
import ProjectEulerToolkit.utils


__all__ = [
    'combinatoric',
    'equation',
    'formula',
    'generator',
    'prime',
    'utils',
]