#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is part of CMInject
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# If you use this program for scientific work, you should correctly reference it; see LICENSE file for details.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

"""
Utility functions for improving performance.
"""

from functools import wraps, lru_cache

import numpy as np


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
        def wrapper(_self, np_array, *args, **kwargs):
            hashable_array = tuple(np_array)
            return cached_wrapper(_self, hashable_array, *args, **kwargs)

        @lru_cache(*args, **kwargs)
        def cached_wrapper(_self, hashable_array, *args, **kwargs):
            array = np.array(hashable_array)
            return function(_self, array, *args, **kwargs)

        # copy refs for lru_cache attributes over
        wrapper.cache_info = cached_wrapper.cache_info
        wrapper.cache_clear = cached_wrapper.cache_clear
        return wrapper

    return decorator


### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
