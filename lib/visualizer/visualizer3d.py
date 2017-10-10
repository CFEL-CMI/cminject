import numpy as np
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class visualizer:
  def __init__(self, exp):
    self.exp = exp
    self.n = len(self.exp.source.particles[0].trajectory)
    self.n_particles = self.exp.source.particles
    self.plot()

  def data_gen(self, i):
        datax=[]
        datay=[]
        dataz=[]
        for j in self.n_particles:
           datax.append(j.trajectory[i][0])
           datay.append(j.trajectory[i][1])
           dataz.append(j.trajectory[i][2])
        return datax, datay, dataz

  def init(self):
    ax.set_ylim( 0.005, 0.015 )
    ax.set_xlim( 0, 0.5 )
    return line,


  def plot(self):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(0, 0.1)
    ax.set_ylim(0, 0.03)
    ax.set_zlim(0, 0)
    for i in range(self.n):
     if plt.fignum_exists(fig.number):
      x, y, z = self.data_gen(i)
      ax.set_ylim(min(y), max(y))
      p = ax.plot(x, y, z, 'o-', lw=0)
      plt.pause(0.01)
      p.pop(0).remove()
     else:
       plt.close()
