import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm, colors
from math import pi


#fig, ax = plt.subplots()
#line, = ax.plot([], [], 'o-', lw=0)
#ax.grid()

class visualizer:
  def __init__(self, exp):
    self.exp = exp
    self.n = exp.LongestTrajectory()
    self.n_particles = self.exp.source.particles
    b = BeamVisualizer(self.exp.beam, None)
    x = np.arange(0, 0.00005, 0.000002)
    y = np.arange(0, 0.00005, 0.000002)
    z = np.arange(-0.005, 0.005, 0.000003)
    R = (x**2 + y**2)
    self.RR, self.ZZ = np.meshgrid(R, z, sparse=True)
    self.I=b.VortexIntensity(self.RR,self.ZZ)

    self.plot()
#    self.visualize()

  def data_gen(self,i):
        datax=[]
        datay=[]
        for j in self.exp.source.particles:
          if len(j.trajectory)>i:
           datax.append(j.trajectory[i][0])
           datay.append(j.trajectory[i][2])
        return datax, datay

  def init(self, ax):
    b = BeamVisualizer(self.exp.beam, None)
    ax.set_xlim( -0.00005, 0.00005 )
    ax.set_ylim( -0.005, 0.005 )
    p = ax.pcolor(np.sqrt(self.RR), self.ZZ, self.I, cmap=cm.hot, vmin=abs(self.I).min(), vmax=abs(self.I).max())
    p1 = ax.pcolor(-1*np.sqrt(self.RR), self.ZZ, self.I, cmap=cm.hot, vmin=abs(self.I).min(), vmax=abs(self.I).max())
#    cb = fig.colorbar(p)

#    return line,


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

  def plot(self):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    self.init(ax)

    FFMpegWriter = animation.writers['ffmpeg']
    metadata = dict(title='Movie Test', artist='Matplotlib',
                comment='Movie support!')
    writer = FFMpegWriter(fps=15, metadata=metadata)
    with writer.saving(fig, "optical_force_5w.mp4", 100):



     for i in range(self.n):
       if plt.fignum_exists(fig.number):
         x, y = self.data_gen(i)
         p = ax.plot(x, y, 'o-', lw=0)
         writer.grab_frame()
         plt.pause(0.005)
         p.pop(0).remove()
       else:
         plt.close()


  def visualize(self):
    ani = animation.FuncAnimation(fig, self.run, self.data_gen, blit=True, interval=10,
                              repeat=False, init_func=self.init)
    plt.show()

class BeamVisualizer:
   def __init__(self, beam, particle):
     self.beam = beam
     self.particle = particle
#     x = np.arange(0, 0.00004, 0.0000005)
#     y = np.arange(0, 0.00004, 0.0000005)
#     z = np.arange(-0.003, 0.003, 0.000005)
#     R = (x**2 + y**2)
#     RR, ZZ = np.meshgrid(R, z, sparse=True)
#     I=self.VortexIntensity(RR,ZZ)
#     rstride = 20
#     cstride = 10
#     fig = plt.figure()
#     ax = fig.add_subplot(111)
#     p = ax.pcolor(np.sqrt(RR), ZZ, I, cmap=cm.hot, vmin=abs(I).min(), vmax=abs(I).max())
#     p1 = ax.pcolor(-1*np.sqrt(RR), ZZ, I, cmap=cm.hot, vmin=abs(I).min(), vmax=abs(I).max())
#     cb = fig.colorbar(p)
#     circle=plt.Circle((np.sqrt(self.particle.position[0]**2 + self.particle.position[1]**2), self.particle.position[2]), radius=self.particle.radius, color='white')
#     ax.add_artist(circle)
#     plt.show()


   def VortexIntensity(self, r2, z):
      """Vortex intensity distribution propagating along the z axis"""
      z0 = pi * (self.beam.w0)**2 / self.beam.lamda    # Beyond the Rayleigh range of
      w =  (self.beam.w0) * np.sqrt( 1 + (z/z0)**2 )
      I = ( 4  / pi ) * r2 / w**4 * np.exp( -2 * r2 / w**2 )
      return I

