from particle import *
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

class Detector:
  def __init__(self, location):
    self.l = location
    self.x=[];self.y=[];self.z=[];
    self.vx=[];self.vy=[];self.vz=[]
    self.particles=[]
    self.after = abs(location)+abs(location)/20.0
    self.before = abs(location)-abs(location)/20.0

  def Projection(self):
    x = np.array(self.x); y = np.array(self.y); z=np.array(self.z)
    vx = np.array(self.vx); vy = np.array(self.vy); vz=np.array(self.vz)

    f0 = interp1d(z, x); f1 = interp1d(z, y); f2 = interp1d(z, vx)
    f3 = interp1d(z, vy); f4 = interp1d(z, vz)
    self.particles.append((f0(self.l), f1(self.l), f2(self.l), f3(self.l), f4(self.l)))
    return
  
  def PlotDetector(self, directory):
    f = open(directory+'detector.txt', 'w+')
    x=[]; y=[]; vx=[]; vy=[]; vz=[]
    for i in self.particles:
      x.append(i[0])
      y.append(i[1])
      vx.append(i[2])
      vy.append(i[3])
      vz.append(i[4])
      f.write(str(i[0])+' '+str(i[1])+' '+str(i[2])
           +' '+str(i[3])+' '+str(i[4])+'\n')
    f.close()
    fig, ax = plt.subplots( nrows=1, ncols=1 )  # create figure & 1 axis
    ax.scatter(x,y)
    fig.savefig(directory+'detector.png')
    plt.close(fig)
    return
    

