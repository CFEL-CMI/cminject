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

from scipy import interpolate
import numpy as np

class field:
    def __init__(self, Type='read', fname=None):
        if Type=='read':
            f = open(fname)
            EE=[]
            Tt=[]
            for i in f:
                token = i.split()
                Tt.append(float(token[0]))
                EE.append(float(token[1]))
            Tt = np.array(Tt)
            EE = np.array(EE)
            self.F = interpolate.interp1d(Tt,EE)
    
    def E(self, t, a):
        try:
            an = a * self.F(t)
        except:
            an = 0.
        return an

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
