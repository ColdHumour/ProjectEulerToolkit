# -*- coding: utf-8 -*-

"""
utils.py

Utility functions using to analysis and evaluate solutions
Function list:
    timepast
    memoize
    clear_cython_cache
    find_solution
    generate_catalogue_draft
    find_keyword_in_solutions

@author: Jasper Wu
"""

import os
import functools
import shutil
import time
import webbrowser
import inspect
from inspect import getmembers, signature


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


def generate_catalogue_draft(module):
    def print_func_info(name, func):
        try:
            sig = str(signature(func))
        except:
            sig = "(<unknown>)"

        name += sig
        doc = func.__doc__
        if not doc:
            doc = name + "\n"
        elif sig == "(<unknown>)":
            doc = "{}\n".format(doc.replace("\n\n", "\n").replace("\n", "\n    "))
        elif not doc.startswith("\n"):
            doc = "{}\n    {}\n".format(name, doc)
        else:
            doc = name + doc
        print(doc)

    blocklist = ["Rational", "Symbol", "deepcopy"]
    basic_pool = dict(getmembers(module, inspect.isfunction))
    basic_pool.update(dict(getmembers(module, inspect.isbuiltin)))
    basic_pool.update(dict(getmembers(module, inspect.isclass)))
    
    print("{}:\n".format(module.__name__))
    for name in sorted(basic_pool):
        print_func_info(name, basic_pool[name])

    for submodule_name in ["equation", "linalg", "numclass", "polynomial"]:
        submodule = getattr(module, submodule_name)
        print("\n{}:\n".format(submodule.__name__))

        pool = getmembers(submodule, inspect.isfunction) + \
               getmembers(submodule, inspect.isbuiltin) + \
               getmembers(submodule, inspect.isclass)
        for name, func in sorted(pool):
            if name not in basic_pool and name not in blocklist and not name.startswith("_"):
                print_func_info(name, func)


def find_keyword_in_solutions(keyword):
    keyword = keyword.lower()

    folder = r"E:\Project Euler\ProjectEulerSolutions"
    solutions = [file for file in os.listdir(folder) if file.endswith("ipynb")]

    check = []
    for file in solutions[::-1]:
        with open(os.path.join(folder, file), 'r', encoding="utf-8") as data:
            data = data.readlines()

        for i, line in enumerate(data):
            if keyword in line.lower():
                check.append(file.replace("ipynb", "html"))
                break

    for file in check:
        print(file)


def find_solution(id):
    if id % 10 == 0:
        lid, rid = id-9, id
    else:
        lid, rid = id // 10 * 10 + 1, id // 10 * 10 + 10
    # webbrowser.open('http://htmlpreview.github.io/?https://github.com/ColdHumour/ProjectEulerSolutions/blob/master/Solutions%20{}-{}.html#{}'.format(lid, rid, id))
    webbrowser.open('file:///E:/Project Euler/ProjectEulerSolutions/Solutions%20{}-{}.html#{}'.format(lid, rid, id))
