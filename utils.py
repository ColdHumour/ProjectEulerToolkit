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
        print "Time (s) consumed by {0}(): {1}".format(func.__name__, time.time() - t)
        return ret
    return _deco