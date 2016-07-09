# -*- coding: utf-8 -*-

"""
utils.py

Utility functions using to analysis and evaluate solutions
Function list: 
    timepast
    clear_cython_cache

@author: Jasper Wu
"""

import os
import shutil


def timepast(func):
    import time
    def _deco(*args, **kwargs):
        t = time.time()
        ret = func(*args, **kwargs)
        print("Time consumed by {0}(): {1}s".format(func.__name__, round(time.time() - t, 2)))
        return ret
    return _deco


def clear_cython_cache(url="C:\\Users\\yudi.wu\\.ipython\\cython"):
    if os.path.exists(url):
        for f in os.listdir(url):
            filepath = os.path.join(url, f)
            if '.' in f:
                os.remove(filepath)
            else:
                shutil.rmtree(filepath)
    else:
        raise ValueError("Bad cython cache URL!")