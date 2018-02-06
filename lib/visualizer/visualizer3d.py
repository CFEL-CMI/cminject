from device.als import *
from device.bgc import *
import numpy as np
#import matplotlib
#matplotlib.use("Agg")
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class visualizer:
  def __init__(self, exp):
    self.exp = exp
    self.n = exp.LongestTrajectory()
    self.n_particles = self.exp.source.particles
    self.plot()

  def data_gen(self, i):
        datax=[]
        datay=[]
        dataz=[]
        for j in self.n_particles:
         if len(j.trajectory)>i:
           datax.append(j.trajectory[i][0])
           datay.append(j.trajectory[i][1])
           dataz.append(j.trajectory[i][2])
        return datax, datay, dataz

  def plot_device(self, ax):
    for i in self.exp.devices:
      if type(i) is AerodynamicsLensStack:
        print "Plot ALS"
        offset = i.position[2]
        for j in i.segments:
          x = np.linspace(i.position[0]+j[0], i.position[0]-j[0], 100)
          z = np.linspace(offset, offset+j[1], 100)
          Xc, Zc=np.meshgrid(x, z)
          Yc = np.sqrt(j[0]*j[0]-Xc**2)
          offset = offset + j[1]
      
      # Draw parameters
          rstride = 20
          cstride = 10
          ax.set_xlim(-0.01, 0.01)
          ax.set_xlabel("X")
          ax.set_ylim(-0.01, 0.01)
          ax.set_ylabel("Y")
          ax.set_zlim(0, -0.05)


          ax.plot_surface(Xc, Yc, Zc, alpha=0.2, rstride=rstride, cstride=cstride)
          ax.plot_surface(Xc, -Yc, Zc, alpha=0.2, rstride=rstride, cstride=cstride)
          ax.view_init(elev=5., azim=-55)

      if type(i) is BufferGasCell:
        yvec = np.arange(0., 0.035, 0.00005)
        zvec = np.arange(0., 0.035, 0.00005)
        Y, Z = np.meshgrid(yvec, zvec)
        X = np.sqrt((Z-0.015)**2 + (Y-0.015)**2)-0.015**2
        X[X < 0.002] = np.nan
        X[X > 0.012] = np.nan
        X1 = -(X-0.042)
        rstride = 20
        cstride = 10
        Xc = np.sqrt(0.015**2-(Y-0.015)**2) +0.015
        #ax.plot_wireframe(X, Y, Z)
        ax.plot_surface(X, Y, Z, alpha=0.2, rstride=rstride, cstride=cstride)
        ax.plot_surface(X1, Y, Z, alpha=0.2, rstride=rstride, cstride=cstride)
        Z[Z<0.012]=np.nan
        Z[Z>0.032]=np.nan
        ax.plot_surface(Z, Y, Xc, alpha=0.2, rstride=rstride, cstride=cstride)
        ax.plot_surface(Z, Y, 0.03-Xc, alpha=0.2, rstride=rstride, cstride=cstride)
        ax.view_init(elev=5., azim=-90)
       # ax.set_zlim(0,2)
        ax.set_xlim(-0.0, 0.045)
        ax.set_xlabel("X")
        ax.set_ylim(-0.0, 0.030)
        ax.set_ylabel("Y")
        ax.set_zlim(0, 0.03)

  def plot(self):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    self.plot_device(ax)


#    FFMpegWriter = animation.writers['ffmpeg']
#    metadata = dict(title='Movie Test', artist='Matplotlib',
#                comment='Movie support!')
#    writer = FFMpegWriter(fps=15, metadata=metadata)
#    with writer.saving(fig, "writer_test.mp4", 100):
    for i in range(self.n):
       if plt.fignum_exists(fig.number):
         x, y, z = self.data_gen(i)
         p = ax.plot(x, y, z, 'o-', lw=0)
#         writer.grab_frame()
         plt.pause(0.01)
         p.pop(0).remove()
   
       else:
         plt.close()
