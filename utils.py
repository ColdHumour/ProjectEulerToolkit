# -*- coding: utf-8 -*-

"""
utils.py

Utility functions using to analysis and evaluate solutions
Function list: 
    timepast

@author: Jasper Wu
"""


def timepast(func):
    import time
    def _deco(*args, **kwargs):
        t = time.time()
        ret = func(*args, **kwargs)
        print "Time consumed by {0}(): {1}s".format(func.__name__, round(time.time() - t, 2))
        return ret
    return _deco