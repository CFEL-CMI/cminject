#!/usr/bin/env python3
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

from integrate import *
from field import field
I = np.array([8.83613997e-43, 2.48769068e-39, 2.48799524e-39]) # Inertia
P = np.array([[231.29e-30,   4.23e-30, 0.43e-30],  
              [  4.23e-30, 216.56e-30, 0.29e-30],
              [  0.43e-30, 0.29e-30, 201.17e-30]])  # Polarizability Tensor
T = 0.4 #K
n_mol = 100 # number of molecules
startTime = 1.e-10
endTime = 3.e-10
f = field(Type='read', fname='laser_pulse.txt') # read the field in Tera Watt per cm^2
dt = abs((np.pi/320.)/np.max(speed(T, I)))  # calculate the time step
costhetaP, costhetaI, time = integrate_all_particles(n_mol, I, P, f, dt, startTime, endTime) # propagate
import matplotlib.pyplot as plt
plt.plot(time, costhetaP)
plt.show()

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
