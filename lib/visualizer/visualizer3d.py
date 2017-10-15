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

    offset = self.exp.als.position[2]
    for i in self.exp.als.segments:
      x = np.linspace(self.exp.als.position[0]+i[0], self.exp.als.position[0]-i[0], 100)
      z = np.linspace(offset, offset+i[1], 100)
      Xc, Zc=np.meshgrid(x, z)
      Yc = np.sqrt(i[0]*i[0]-Xc**2)
      offset = offset + i[1]
      
      # Draw parameters
      rstride = 20
      cstride = 10
      ax.plot_surface(Xc, Yc, Zc, alpha=0.2, rstride=rstride, cstride=cstride)
      ax.plot_surface(Xc, -Yc, Zc, alpha=0.2, rstride=rstride, cstride=cstride)


  def plot(self):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(elev=5., azim=-55)
    self.plot_device(ax)
    ax.set_xlim(-0.01, 0.01)
    ax.set_xlabel("X")
    ax.set_ylim(-0.01, 0.01)
    ax.set_ylabel("Y")
    ax.set_zlim(0, -0.05)



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
