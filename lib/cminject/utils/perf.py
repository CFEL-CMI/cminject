#!/usr/bin/env python3
# -*- coding: utf-8; fill-column: 100; truncate-lines: t -*-
#
# This file is part of CMInject
#
# Copyright (C) 2018,2020 CFEL Controlled Molecule Imaging group
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public 
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any 
# later version. 
#
# If you use this program for scientific work, you should correctly reference it; see the LICENSE.md file for details. 
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. 
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

"""
Utility functions for improving performance.
"""
from _thread import RLock
from functools import wraps, lru_cache

import numpy as np


_NOT_FOUND = object()


class cached_property:
    """
    Copied verbatim from the Python 3.8 source code. See the documentation for functools.cached_property there.
    The reason it lives here is that this is the only part in CMInject currently requiring Python >3.6, and by copying
    over this self-contained implementation, we can keep the version requirement at >3.6.
    """
    def __init__(self, func):
        self.func = func
        self.attrname = None
        self.__doc__ = func.__doc__
        self.lock = RLock()

    def __set_name__(self, owner, name):
        if self.attrname is None:
            self.attrname = name
        elif name != self.attrname:
            raise TypeError(
                "Cannot assign the same cached_property to two different names "
                f"({self.attrname!r} and {name!r})."
            )

    def __get__(self, instance, owner=None):
        if instance is None:
            return self
        if self.attrname is None:
            raise TypeError(
                "Cannot use cached_property instance without calling __set_name__ on it.")
        try:
            cache = instance.__dict__
        except AttributeError:  # not all objects have __dict__ (e.g. class defines slots)
            msg = (
                f"No '__dict__' attribute on {type(instance).__name__!r} "
                f"instance to cache {self.attrname!r} property."
            )
            raise TypeError(msg) from None
        val = cache.get(self.attrname, _NOT_FOUND)
        if val is _NOT_FOUND:
            with self.lock:
                # check if another thread filled cache while we awaited lock
                val = cache.get(self.attrname, _NOT_FOUND)
                if val is _NOT_FOUND:
                    val = self.func(instance)
                    try:
                        cache[self.attrname] = val
                    except TypeError:
                        msg = (
                            f"The '__dict__' attribute on {type(instance).__name__!r} instance "
                            f"does not support item assignment for caching {self.attrname!r} property."
                        )
                        raise TypeError(msg) from None
        return val


def numpy_method_cache(*args, **kwargs):
    """
    LRU cache (functools.lru_cache) based caching implementation for functions whose SECOND parameter is a numpy array.
    Mostly useful for instance methods (whose first parameter is `self`) with one np.array as their second argument.

    Based on https://gist.github.com/Susensio/61f4fee01150caaac1e10fc5f005eb75, this is a more specialised and thus
    simpler and faster implementation for 1D arrays that can be converted into tuples without recursing down.
    More general ways might need to be derived depending on the needs.
    """

    def decorator(function):
        @wraps(function)
        def wrapper(_self, np_array, *wrapped_args, **wrapped_kwargs):
            hashable_array = tuple(np_array)
            return cached_wrapper(_self, hashable_array, *wrapped_args, **wrapped_kwargs)

        @lru_cache(*args, **kwargs)
        def cached_wrapper(_self, hashable_array, *wrapped_args, **wrapped_kwargs):
            array = np.array(hashable_array)
            return function(_self, array, *wrapped_args, **wrapped_kwargs)

        # copy refs for lru_cache attributes over
        wrapper.cache_info = cached_wrapper.cache_info
        wrapper.cache_clear = cached_wrapper.cache_clear
        return wrapper

    return decorator
