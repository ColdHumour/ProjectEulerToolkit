# -*- coding: utf-8 -*-

"""
Cython extensions
File list:
    cpp_types.pxd
    c_formula_int64.pyx/pxd -> c_formula_int64.pyd
    c_linalg_int64.pyx/pxd -> c_linalg_int64.pyd
    c_prime_int64.pyx/pxd -> c_prime_int64.pyd
    cpp_formula_int64.pyx/pxd -> cpp_formula_int64.pyd
    cpp_prime_int64.pyx/pxd -> cpp_prime_int64.pyd
    simple_bigint.pyx/pxd -> simple_bigint.pyd

@author: Jasper Wu
"""

import os
import shutil

try:
    from . import (
        c_formula_int64,
        c_linalg_int64,
        c_prime_int64,
        cpp_formula_int64,
        cpp_prime_int64,
        simple_bigint,
    )
except:
    CUR_DIR = os.getcwd()
    EXT_DIR = os.path.dirname(os.path.abspath(__file__))

    # remove .pyd files since building process is for all files
    for f in os.listdir(EXT_DIR):
        faffix = f.split('.')
        if faffix[-1] == 'pyd':
            os.remove(os.path.join(EXT_DIR, f))

    # build cython extensions
    # the build command is executed at CUR_DIR, thus extsetup.py path must be hardcoded
    # ~/build directory will also created at CUR_DIR, but .c and .cpp files are created
    # at the same level of extsetup.py
    state = os.system('python ProjectEulerToolkit/ext/extsetup.py build_ext -i --compiler=msvc')

    if state == 0:
        # clean intermediate files
        build = os.path.join(CUR_DIR, 'build')
        if os.path.exists(build):
            shutil.rmtree(build)
        for f in os.listdir(EXT_DIR):
            faffix = f.split('.')
            if faffix[-1] in 'c|cpp':
                os.remove(os.path.join(EXT_DIR, f))

        from . import (
            c_formula_int64,
            c_linalg_int64,
            c_prime_int64,
            cpp_formula_int64,
            cpp_prime_int64,
            simple_bigint,
        )
    else:
        print("WARNING: Fail to build Cython extensions!")
