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

import sys
sys.path.insert(0, '../lib')
from field.LBMD2Q9 import LBM
import threading

def setup_viz():
    import matplotlib
    import matplotlib.figure
    import matplotlib.backends.backend_tkagg
    import Tkinter

    viz = {}

    # Initialise Tk ...
    figure = matplotlib.figure.Figure()
    figure.set_size_inches((8, 6))

    viz['axes1'] = figure.add_subplot(311)
    viz['axes2'] = figure.add_subplot(312)
    viz['axes3'] = figure.add_subplot(313)
    viz['tk'] = Tkinter.Tk()
    viz['canvas'] = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(
        figure,
        master=viz['tk']
    )
    viz['canvas'].get_tk_widget().pack(expand = True, fill = Tkinter.BOTH)

    return viz


Xdim = 3 
Ydim = 1
dx = 0.02
Umax = 0.01
dt = 1
time = 1000
Re = 10
eta = 0.01
OutR = 0.10
als = [  (0.5, 0.35) , (1, 0.3), (1.5, 0.25), (2, 0.20), (2.5, 0.18)  ] 

viz = setup_viz()
thread = threading.Thread(target = LBM, args=( Xdim, Ydim, dx, Umax, dt, time, Re, eta, als, OutR, viz))
thread.start()
viz['tk'].mainloop()

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
