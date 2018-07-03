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
