# -*- coding: utf-8 -*-

"""
utils.py

Utility functions using to analysis and evaluate solutions
Function list:
    timepast
    memoize
    clear_cython_cache
    find_solution

@author: Jasper Wu
"""

import os
import functools
import shutil
import time
# import webbrowser


def timepast(func):
    @functools.wraps(func)
    def _deco(*args, **kwargs):
        t = time.time()
        print("Func {0}() begins at: {1}".format(func.__name__, time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(t))))
        ret = func(*args, **kwargs)
        print("Time consumed by {0}(): {1:.2f}s".format(func.__name__, time.time() - t))
        return ret
    return _deco


def memoize(cache={}, key=lambda x: x):
    @functools.wraps(func)
    def _deco(func):
        def __deco(*args, **kwargs):
            idx = key(args)
            if idx not in cache:
                cache[idx] = func(*args, **kwargs)
            return cache[idx]
        return __deco
    return _deco


def clear_cython_cache():
    url = os.path.expanduser("~\\.ipython\\cython")
    if os.path.exists(url):
        for f in os.listdir(url):
            filepath = os.path.join(url, f)
            if '.' in f:
                os.remove(filepath)
            else:
                shutil.rmtree(filepath)
    else:
        raise ValueError("Cython cache not found: {}".format(url))


# def find_solution(id):
#     if id % 10 == 0:
#         lid, rid = id-9, id
#     else:
#         lid, rid = id // 10 * 10 + 1, id // 10 * 10 + 10
#     webbrowser.open('http://htmlpreview.github.io/?https://github.com/ColdHumour/ProjectEulerSolutions/blob/master/Solutions%20{}-{}.html#{}'.format(lid, rid, id))
