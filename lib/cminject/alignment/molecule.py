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

import numpy as np 
from pyquaternion import Quaternion

class molecule(object):
        """Class to carry the data of the molecule"""
        def __init__(self, I, P, T=0.5, Q=None, W=None):
                self.I = I
                self.P = P
                if Q is None:
                    self.Q = Quaternion.random() 
                else: 
                    self.Q = Q
                if W is None:
                   k = 1.38042e-23
                   mu = np.sqrt(2*k*T/I)
                   sig = np.sqrt(k*T/I)
                   self.W = np.random.normal(mu, sig)
                self.PA = self.get_PA()

        def update_w(self, E, a, dt):
                """Calculate the torque then update 
                the angular velocities"""
                dipole = np.dot(self.P, E)
                torque = np.cross(a*E, dipole)
                dw = self.dw(self.W, torque)
                w1 = self.W + dw * dt/2.
                dw = self.dw(w1, torque)
                self.W = self.W + dw * dt

        def dw(self, w, T):
                """Calculate the acceleration given the torque"""
                dw = np.zeros((3))
                dw[0] = (T[0] -(self.I[2] - self.I[1])*w[1]*w[2])/self.I[0]
                dw[1] = (T[1] -(self.I[0] - self.I[2])*w[0]*w[2])/self.I[1]
                dw[2] = (T[2] -(self.I[1] - self.I[0])*w[0]*w[1])/self.I[2]
                return dw

        def update_q(self, dt):
                """the derivatives of the quatrionions given omega"""
                dq = self.Q.derivative(self.W)
                self.Q = self.Q + dq * dt
                self.Q = self.Q/self.Q.norm
        def get_PA(self):
            autoval, autovect = np.linalg.eig(self.P)
            auto_ord = np.argsort(autoval)
            ord_autoval = autoval[auto_ord]
            ord_autovect_complete = autovect[:, auto_ord].T
            return ord_autovect_complete[0]

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
