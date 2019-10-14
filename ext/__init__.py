# -*- coding: utf-8 -*-

"""
Cython extensions
File list:
    c_formula_int64.pyx/pxd -> c_formula_int64.pyd
    c_linalg_int64.pyx/pxd -> c_linalg_int64.pyd
    c_prime_int64.pyx/pxd -> c_prime_int64.pyd
    cpp_formula_int64.pyx/pxd -> cpp_formula_int64.pyd

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
    )
except:
    CUR_DIR = os.getcwd()
    EXT_DIR = os.path.dirname(os.path.abspath(__file__))

    # build cython extensions
    os.chdir(EXT_DIR)
    state = os.system('python extsetup.py build_ext -i --compiler=msvc')

    if state == 0:
        # clean intermediate files
        build = os.path.join(EXT_DIR, 'build')
        if os.path.exists(build):
            shutil.rmtree(build)
        for f in os.listdir(EXT_DIR):
            faffix = f.split('.')
            if faffix[-1] in 'c|cpp':
                os.remove(os.path.join(EXT_DIR, f))

        # change back to workdir
        os.chdir(CUR_DIR)

        from . import (
            c_formula_int64,
            c_linalg_int64,
            c_prime_int64,
            cpp_formula_int64,
        )
    else:
        # change back to workdir
        os.chdir(CUR_DIR)

        print("WARNING: Fail to build Cython extensions!")
