# -*- coding: utf-8 -*-

"""
Cython extensions
File list:
    c_formula_int64.pyx/pxd -> c_formula_int64.pyd
    c_linalg_int64.pyx/pxd -> c_linalg_int64.pyd
    c_prime_int64.pyx/pxd -> c_prime_int64.pyd

@author: Jasper Wu
"""

import os
import shutil

try:
    from . import (
        c_formula_int64,
        c_linalg_int64,
        c_prime_int64,
    )
except:
    CUR_DIR = os.getcwd()
    EXT_DIR = os.path.dirname(os.path.abspath(__file__))

    # build cython extensions
    os.chdir(EXT_DIR)
    state = os.system('python extsetup.py build_ext -i --compiler=msvc')

    if state == 0:
        # clean intermediate files
        for folder in ("build", "__pycache__"):
            folder = os.path.join(EXT_DIR, folder)
            if os.path.exists(folder):
                shutil.rmtree(folder)
        for f in os.listdir(EXT_DIR):
            faffix = f.split('.')
            if faffix[-1] == 'c':
                os.remove(os.path.join(EXT_DIR, f))

        # change back to workdir
        os.chdir(CUR_DIR)

        from . import (
            c_formula_int64,
            c_prime_int64,
        )
    else:
        # change back to workdir
        os.chdir(CUR_DIR)

        print("WARNING: Fail to build Cython extensions!")
