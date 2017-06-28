# -*- coding: utf-8 -*-

"""
utils.py

Utility functions using to analysis and evaluate solutions
Function list:
    timepast
    memoize
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


def memoize(cache=None, key=lambda x: x):
    if cache is None:
        raise ValueError("cache must be an existed dict!")

    def _deco1(func):
        def _deco2(*args, **kw):
            idx = key(args)
            if idx not in cache:
                cache[idx] = func(*args, **kw)
            return cache[idx]
        return _deco2
    return _deco1


def clear_cython_cache(url="C:\\Users\\wuyd\\.ipython\\cython"):
    if os.path.exists(url):
        for f in os.listdir(url):
            filepath = os.path.join(url, f)
            if '.' in f:
                os.remove(filepath)
            else:
                shutil.rmtree(filepath)
    else:
        raise ValueError("Bad cython cache URL!")
