from particle import *
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp1d
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
from matplotlib import rc
import matplotlib.pyplot as plt
plt.ioff()
from scipy import optimize
from scipy.stats import norm

class Detector:
  """ This class implement detectors. It obtain the projection of the particles on
      the detectors supplied by the users"""

  def __init__(self, locations):
    self.l = locations        # locations of detectors provided by the user
    self.n = len(locations)   # number of detectors
    self.x=[];self.y=[];self.z=[];self.ID=0;self.Tt=[]  # saving positions and velocities needed 
    self.vx=[];self.vy=[];self.vz=[];self.T=[]  # for interpolation
    self.particles=[]         # save phase space for all particles
    locations = np.absolute(locations)
    maximum = max(locations)
    self.end = maximum + maximum/20.0  # This is to know where to stop the particles
                                       # the particles should stop at the last detector
  def Check(self, z):
    """ This function checks if the particle is near any of the detectors to save
        its position for the interpolation"""
    for i in self.l:
      after = abs(i)+abs(i)/20.0
      before = abs(i)-abs(i)/20.0
      if z > before and z < after:
        return True
    return False
 
  def Projection(self):
    """This function des the interpolation to obtain physical properties at the detector"""
    x = np.array(self.x); y = np.array(self.y); z=np.array(self.z)
    vx = np.array(self.vx); vy = np.array(self.vy); vz=np.array(self.vz)
    T = np.array(self.T); Tt = np.array(self.Tt)
    f0 = interp1d(z, x); f1 = interp1d(z, y); f2 = interp1d(z, vx)
    f3 = interp1d(z, vy); f4 = interp1d(z, vz); f5 = interp1d(z, T); f6 = interp1d(z, Tt)
    self.particles.append((f0(self.l), f1(self.l), f2(self.l), f3(self.l), f4(self.l), f5(self.l), f6(self.l), self.ID))
    return

  def Gaussian2D(self, xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    """ This function is required for the 2D fittings of the detectors"""
    (x, y) = xdata_tuple
    x = np.array(x)
    y = np.array(y)
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)
                        + c*((y-yo)**2)))
    return g.ravel()

  def FitandPlot(self, x, y, location, phase, directory):
    x = np.array(x)
    y = np.array(y)
    (mux, sigmax) = norm.fit(x)
    (muy, sigmay) = norm.fit(y)
    n = 150
    xedges, yedges = np.linspace(min(x), max(x), n), np.linspace(min(y), max(y), n)
    H_rot, x, y = np.histogram2d(x,y, (xedges, yedges))
    X, Y = np.meshgrid(x[1:], y[1:])
    initial_guess = (np.amax(H_rot), mux, muy, sigmax, sigmay, 0., 0.)
    popt, pcov = optimize.curve_fit(self.Gaussian2D, (X, Y), H_rot.ravel(), p0=(initial_guess[0],
               initial_guess[1],initial_guess[2],initial_guess[3],initial_guess[4],initial_guess[5],
                 initial_guess[6]),maxfev=100000)
    data_fitted = self.Gaussian2D((X, Y), *popt)
    fig, ax = plt.subplots(1, 1)
    cax = ax.imshow(data_fitted.reshape(n-1,n-1), cmap=plt.cm.jet, origin='bottom',
                                                     extent=(X.min(), X.max(), Y.min(), Y.max()))
    ax.contour(X, Y, data_fitted.reshape(n-1,n-1), 8, colors='w')
    cbar = fig.colorbar(cax)
    ax.set_title(phase +" Detector at " + location + " Z" )
    txt ='X sigma '+str(round(popt[3],7)) +'Y sigma '+str(round(popt[4],7))+'Angle '+str(popt[5])
    rc('text', usetex=True)
    ax.set_ylabel('Y-axis')
    ax.set_xlabel(r'\begin{center}X-axis\\*\textit{\small{' + txt + r'}}\end{center}') 
#    fig.text(.5, .05, txt, ha='center')
    plt.tight_layout()
    fig.savefig(directory+'/images/'+location+'_'+phase+'_detector.png')
    plt.close(fig)

  def PlotDetector(self, directory, filename, source):
    """This function save the detectors' data to hdf5 file
       it also draw the detectors with the 2D fitting"""

    print len(self.particles), "Particles servived"
    f = h5py.File(directory+filename+".hdf5", "w")
    detectors = f.create_group("detectors")
    for j in range(self.n):
      x=[]; y=[]; vx=[]; vy=[]; vz=[]
      ID=[]; T=[]; Tt=[]
      for i in self.particles:
        x.append(i[0][j])
        y.append(i[1][j])
        vx.append(i[2][j])
        vy.append(i[3][j])
        vz.append(i[4][j])
        T.append(i[5][j])
        Tt.append(i[6][j])
        ID.append(i[7])
      data = np.array([x,y,vx,vy,vz,T, Tt, ID])
      detectors.create_dataset(str(self.l[j]), data=data)
      self.FitandPlot(x, y, str(self.l[j]), "Positions", directory)
      self.FitandPlot(vx, vy, str(self.l[j]), "Velocities", directory)
    x=[]; y=[]; z=[]; vx=[]; vy=[]; vz=[]; ID=[];
    for i in source.particles:
        x.append(i.position[0])
        y.append(i.position[1])
        z.append(i.position[2])
        vx.append(i.velocity[0])
        vy.append(i.velocity[1])
        vz.append(i.velocity[2])
        ID.append(i.ID)
    data = np.array([x,y,z,vx,vy,vz,ID])
    source = f.create_group("source")
    source.create_dataset("source", data=data)
    f.close()
    return
    

