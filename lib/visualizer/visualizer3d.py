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

  def plot_device(self, ax):
    r1 = [ (-0.01, 0.01), (-0.001, 0.001), (-0.01, 0.01), (-0.001, 0.001)]
    z1 = [ ( 0.0, -0.025), (-0.025, -0.030), (-0.03, -0.045), (-0.045, -0.05)] 

    for i,j in zip(r1,z1):
      x = np.linspace(i[0], i[1], 100)
      z = np.linspace(j[0], j[1], 100)
      Xc, Zc=np.meshgrid(x, z)
      Yc = np.sqrt(i[0]*i[0]-Xc**2)

      # Draw parameters
      rstride = 20
      cstride = 10
      ax.plot_surface(Xc, Yc, Zc, alpha=0.2, rstride=rstride, cstride=cstride)
      ax.plot_surface(Xc, -Yc, Zc, alpha=0.2, rstride=rstride, cstride=cstride)


  def plot(self):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    self.plot_device(ax)
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
