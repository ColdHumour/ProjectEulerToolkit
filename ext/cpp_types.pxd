# -*- coding: utf-8 -*-
# distutils: language = c++

"""
cpp_types.pxd

Declaration file for c++ types accessed by other pyx

@author: Jasper Wu
"""

from libcpp cimport bool
from libcpp.map cimport map as cmap
from libcpp.pair cimport pair
from libcpp.set cimport set as cset
from libcpp.vector cimport vector

ctypedef long long int64
ctypedef pair[int64, int64] lpair
ctypedef cset[int64] lset
ctypedef cset[int64].iterator lset_it
ctypedef vector[bool] bvec
ctypedef vector[int64] lvec

ctypedef lpair sbi  # simple big integer
