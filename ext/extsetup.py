# -*- coding: utf-8 -*-

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

import numpy as np

extensions = []

EXT_FILES = ['c_formula_int64', 'cpp_formula_int64']
for f in EXT_FILES:
    extensions.append(Extension(f, [f + '.pyx']))

EXT_FILES_INC = ['c_linalg_int64', 'c_prime_int64']
for f in EXT_FILES_INC:
    extensions.append(Extension(f, [f + '.pyx'],
                                include_dirs=[np.get_include()]))

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=cythonize(extensions),
)
