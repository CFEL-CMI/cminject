import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation



fig, ax = plt.subplots()
line, = ax.plot([], [], 'o-', lw=0)
ax.grid()

class visualizer:
  def __init__(self, exp):
    self.exp = exp
    self.visualize()

  def data_gen(self):
    for i in range(len(self.exp.source.particles[0].trajectory)):
        datax=[]
        datay=[]
        for j in self.exp.source.particles:
           datax.append(j.trajectory[i][0])
           datay.append(j.trajectory[i][1])
        yield datax, datay

  def init(self):
    ax.set_ylim( 0.005, 0.015 )
    ax.set_xlim( 0, 0.5 )
    return line,


  def run(self, data):
    # update the data
    xdata=[]
    ydata=[]
    x, y = data
    xdata.append(x)
    ydata.append(y)
    xmin, xmax = ax.get_xlim()

    line.set_data(xdata, ydata)

    return line,


  def visualize(self):
    ani = animation.FuncAnimation(fig, self.run, self.data_gen, blit=True, interval=3,
                              repeat=False, init_func=self.init)
    plt.show()

