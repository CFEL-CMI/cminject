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

  def plot(self):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(-0.01, 0.01)
    ax.set_xlabel("X")
    ax.set_ylim(-0.01, 0.01)
    ax.set_ylabel("Y")
    ax.set_zlim(0, -0.05)
    for i in range(self.n):
     if plt.fignum_exists(fig.number):
      x, y, z = self.data_gen(i)
      p = ax.plot(x, y, z, 'o-', lw=0)
      plt.pause(0.01)
      p.pop(0).remove()
     else:
       plt.close()
