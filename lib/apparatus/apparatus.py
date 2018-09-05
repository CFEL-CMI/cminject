# -*- coding: utf-8 -*-
#
# This file is part of CMInject


class apparatus:
    """This class describe the general feature of any experimental apparatus. There
       will be several classes inherited from this class that describe the
       details of the apparatus, such as the ALS and ADL classes
    """
    def __init__(self, x, y, length, width):
        """ For any device the initial coordinates and the length and width has to be defined"""
        self.x = x
        self.y = y
        self.length = length
        self.width = width



### Local Variables:
### fill-column: 80
### truncate-lines: t
### End:
