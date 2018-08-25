# -*- coding: utf-8 -*-

"""
c_linalg_int64.pyx

Cython extension of functions implementing linear algebra functions via fast algorithms.
Function list:
    c_dot_mod_int64
    c_mat_pow_mod_int64
    c_mat_sum_pow_mod_int64

@author: Jasper Wu
"""

import numpy as np
cimport numpy as np

ctypedef long long int64
ctypedef np.ndarray arr

cpdef arr[int64, ndim=2] c_dot_mod_int64(arr[int64, ndim=2] A, arr[int64, ndim=2] B, int64 MOD):
    cdef:
        arr[int64, ndim=2] C
        long w, m, h, i, j, k
        int64 v

    w = len(A)
    m = len(B)
    h = len(B[0])
    C = np.zeros((w, h), dtype=np.int64)
    for i in range(w):
        for j in range(h):
            for k in range(m):
                C[i, j] = (C[i, j] + A[i, k] * B[k, j]) % MOD
    return C

cpdef arr[int64, ndim=2] c_mat_pow_mod_int64(arr[int64, ndim=2] mat, int64 n, int64 m=0):
    """return (mat^n) % m"""

    cdef:
        short d = len(mat)
        arr[int64, ndim=2] res = np.eye(d, dtype=np.int64)

    while n:
        if n & 1:
            if m:
                res = c_dot_mod_int64(res, mat, m)
            else:
                res = np.dot(res, mat)
        if m:
            mat = c_dot_mod_int64(mat, mat, m)
        else:
            mat = np.dot(mat, mat)
        n >>= 1
    return res

cpdef arr[int64, ndim=2] c_mat_sum_pow_mod_int64(arr[int64, ndim=2] A0, arr[int64, ndim=2] Q, int64 n, int64 m=0):
    """return (A0 + Q A0 + Q^2 A0 + ... + Q^n A0) % m"""

    cdef:
        short d = len(Q)
        arr[int64, ndim=2] O, I, Q_ext, Q_ext_pow, I2, res
    
    if n == 0:
        return A0

    if m:
        A0 = A0 % m
        Q = Q % m

    O = np.zeros((d, d), dtype=np.int64)
    I = np.eye(d, dtype=np.int64)
    Q_ext = np.concatenate([np.concatenate([Q, I], axis=1), np.concatenate([O, I], axis=1)])
    Q_ext_pow = c_mat_pow_mod_int64(Q_ext, n, m)
    I2 = np.concatenate([I, I], axis=0)
    if m:
        res = c_dot_mod_int64(c_dot_mod_int64(Q_ext_pow, I2, m), A0, m)
    else:
        res = np.dot(np.dot(Q_ext_pow, I2), A0) 

    return res[:d]
