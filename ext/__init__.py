# -*- coding: utf-8 -*-

"""
Cython extensions
File list:
    _formula.pxd
    _formula.pyx
    _prime.pxd
    _prime.pyx

@author: Jasper Wu
"""

import os
import shutil

try:
    import _formula
    import _prime
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
            if faffix[-1] == 'c':
                os.remove(os.path.join(EXT_DIR, f))

        # change back to workdir
        os.chdir(CUR_DIR)

        import _formula
        import _prime
    else:
        print "WARNING: Fail to build Cython extensions"