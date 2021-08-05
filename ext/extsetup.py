# -*- coding: utf-8 -*-

from importlib import import_module
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

import numpy as np

extensions = []

EXT_FILES = ['c_formula_int64', 'cpp_formula_int64', 'cpp_polynomial_int64', 'cpp_prime_int64', 'simple_bigint']
for f in EXT_FILES:
    extensions.append(Extension("ProjectEulerToolkit.ext.{}".format(f),
                                ["ProjectEulerToolkit/ext/{}.pyx".format(f)]))

EXT_FILES_INC = ['c_linalg_int64', 'c_prime_int64']
for f in EXT_FILES_INC:
    extensions.append(Extension("ProjectEulerToolkit.ext.{}".format(f),
                                ["ProjectEulerToolkit/ext/{}.pyx".format(f)],
                                include_dirs=[np.get_include()],
                                define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]))

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=cythonize(extensions),
)
