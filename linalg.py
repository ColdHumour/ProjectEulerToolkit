# -*- coding: utf-8 -*-

"""
linalg.py

Functions implementing maths in linear algebra.
Function list:
    mat_pow_mod
    gauss_jordan_elimination
    gauss_jordan_elimination_with_unknown_RHS
    get_integer_matrix_inverse_as_list
    get_integer_matrix_inverse_as_numpy_array

@author: Jasper Wu
"""

from copy import deepcopy

import numpy as np
from sympy import Symbol, Rational

from . formula import gcd


def mat_pow_mod(mat, n, m=0):
    """return (mat^n) % m"""

    if n < 0:
        raise ValueError("power must be positive!")

    d = len(mat)
    res = np.eye(d, dtype=np.int64)
    while n:
        if n & 1:
            if m:
                res = np.mod(np.dot(res, mat), m)
            else:
                res = np.dot(res, mat)
        if m:
            mat = np.mod(np.dot(mat, mat), m)
        else:
            mat = np.dot(mat, mat)
        n >>= 1
    return res



def gauss_jordan_elimination(coeffs):
    """
    Gauss-Jordan elimination algorithm, can only be used when there are more variables than equations
    coeffs: 2D list, all elements float
    """

    w, d = len(coeffs[0]), len(coeffs)
    coefmat = np.matrix(coeffs)

    for i in range(d):
        flag = 1
        j = i
        while flag and j < d:
            if abs(coefmat[i, j]) < 0.001:
                for k in range(i+1, d):
                    if abs(coefmat[k, j]) > 0.001:
                        flag = 0
                        coefmat[k], coefmat[i] = deepcopy(coefmat[i]), deepcopy(coefmat[k])
                        break
                if flag:
                    j += 1
            else:
                flag = 0

        if j == d:
            break

        if coefmat[i, j] != 1:
            coefmat[i] /= coefmat[i, j]

        for k in range(i+1, d):
            if coefmat[k, j]:
                coefmat[k] = coefmat[k] - coefmat[k, j] * coefmat[i]

    for i in range(1, d):
        for j in range(w):
            if abs(coefmat[i, j]) > 0.001:
                break
        for k in range(i):
            if coefmat[k, j]:
                coefmat[k] = coefmat[k] - coefmat[k, j] * coefmat[i]

    return coefmat


def gauss_jordan_elimination_with_unknown_RHS(coeffs):
    """
    Gauss-Jordan elimination algorithm, can only be used when there are more variables than equations
    Allow RHS to be sympy.Symbol
    coeffs: 2D list, all elements sympy.Rational or sympy.Symbol
    """

    w, d = len(coeffs[0]), len(coeffs)
    coefmat = deepcopy(coeffs)
    for i in range(d):
        flag = 1
        j = i
        while flag and j < d:
            if isinstance(coefmat[i][j], Rational) and coefmat[i][j] == 0:
                for k in range(i+1, d):
                    if isinstance(coefmat[i][j], Rational) and coefmat[k][j] != 0:
                        flag = 0
                        coefmat[k], coefmat[i] = deepcopy(coefmat[i]), deepcopy(coefmat[k])
                        break
                if flag:
                    j += 1
            else:
                flag = 0

        if j == d:
            break

        if coefmat[i][j] != 1:
            coefmat[i] = [n / coefmat[i][j] for n in coefmat[i]]

        for k in range(i+1, d):
            if coefmat[k][j] != 0:
                coefmat[k] = [coefmat[k][x] - coefmat[k][j] * coefmat[i][x] for x in range(w)]

    for i in range(1, d):
        for j in range(w-1):
            if coefmat[i][j] != 0:
                break
        for k in range(i):
            if coefmat[k][j]:
                coefmat[k] = [coefmat[k][x] - coefmat[k][j] * coefmat[i][x] for x in range(w)]

    return coefmat


def get_integer_matrix_inverse_as_list(matrix):
    """
    get inverse of matrix by elementry row transformation
    the matrix is defined by list of list, to avoid overflow in numpy
    """

    L = len(matrix)
    matrix = [row + [0]*L for row in matrix]
    for r in range(L):
        matrix[r][r+L] = 1

    # handle every row
    for r in range(L):
        r2 = r
        while r2 < L and matrix[r2][r] == 0:
            r2 += 1

        if r2 == L:
            raise ValueError("Singular matrix!")

        if r2 != r:
            matrix[r], matrix[r2] = matrix[r2][:], matrix[r][:]
        if matrix[r][r] < 0:
            matrix[r] = [-x for x in matrix[r]]

        g = matrix[r][r]
        for c in range(r+1, 2*L):
            x = matrix[r][c]
            if x:
                g = gcd(g, abs(x))
        matrix[r] = [x // g for x in matrix[r]]

        pivot = matrix[r][r]
        for r2 in range(L):
            if r2 != r and matrix[r2][r]:
                x = matrix[r2][r]
                if x:
                    g = gcd(abs(x), pivot)
                    matrix[r2] = [pivot // g * matrix[r2][c] - x // g * matrix[r][c] for c in range(2*L)]

    # reduce every row
    det = 1
    for r in range(L):
        if matrix[r][r] < 0:
            matrix[r] = [-x for x in matrix[r]]

        g = matrix[r][r]
        for c in range(L, 2*L):
            x = matrix[r][c]
            if x:
                g = gcd(g, abs(x))
        matrix[r] = [x // g for x in matrix[r]]

        g = matrix[r][r]
        det = g // gcd(g, det) * det

    # handle diagonal
    for r in range(L):
        g = det // matrix[r][r]
        matrix[r] = [g * x for x in matrix[r]]

    # get inverse
    matrix = [row[L:] for row in matrix]

    return det, matrix


def get_integer_matrix_inverse_as_numpy_array(matrix):
    """
    get inverse of matrix by elementry row transformation
    the matrix is defined by numpy array with int64, but may still overflow when values are large
    """

    L = len(matrix)
    I = np.eye(L, dtype=np.int64)
    matrix = np.concatenate((matrix, I), axis=1)

    # handle every row
    for r in range(L):
        r2 = r
        while r2 < L and matrix[r2, r] == 0:
            r2 += 1

        if r2 == L:
            raise ValueError("Singular matrix!")

        if r2 != r:
            matrix[r], matrix[r2] = matrix[r2].copy(), matrix[r].copy()
        if matrix[r, r] < 0:
            matrix[r] *= -1

        g = matrix[r, r]
        for c in range(r+1, 2*L):
            x = matrix[r, c]
            if x:
                g = gcd(g, abs(x))
        matrix[r] //= g

        pivot = matrix[r, r]
        for r2 in range(L):
            if r2 != r and matrix[r2, r]:
                x = matrix[r2, r]
                if x:
                    g = gcd(abs(x), pivot)
                    matrix[r2] = pivot // g * matrix[r2] - x // g * matrix[r]

    # reduce every row
    det = 1
    for r in range(L):
        if matrix[r, r] < 0:
            matrix[r] *= -1

        g = matrix[r, r]
        for c in range(L, 2*L):
            x = matrix[r, c]
            if x:
                g = gcd(g, abs(x))
        matrix[r] //= g

        g = int(matrix[r, r])
        det = g // gcd(g, det) * det

    # handle diagonal
    for r in range(L):
        matrix[r] *= det // matrix[r, r]

    return det, matrix[:, L:]
