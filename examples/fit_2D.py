import sys
import h5py
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from scipy import optimize
from scipy.stats import norm


def Gaussian2D(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    (x, y) = xdata_tuple
    x = np.array(x)
    y = np.array(y)
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g =  amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)
                        + c*((y-yo)**2)))
    return g.ravel()

def filter(x, y, vx, vy, vz, ID, X, Y, VX, VY, VZ):
  if x is not None and y is not None:
    mask = np.zeros(len(x), dtype=bool)
    for i in range(len(ID)):
      if (x[ID[i]]>X[0] and x[ID[i]]<X[1]
            and y[ID[i]]>Y[0] and y[ID[i]]<Y[1]
              and vx[ID[i]]>VX[0] and vx[ID[i]]<VX[1]
                and vy[ID[i]]>VY[0] and vy[ID[i]]<VY[1]
                  and vz[ID[i]]>VZ[0] and vz[ID[i]]<VZ[1]):
        mask[i] = True
    return mask



Xrange = (-1.0, 1.0)
Yrange = (-1.0, 1.0)
VXrange = (-80., 80.)
VYrange = (-80., 80.)
VZrange = (-80., 80.)

x = None
for i in range(3, len(sys.argv)):
  f = h5py.File(sys.argv[i], 'r')  
  d = f['detectors']
  s = f['source']['source']
  det = d[sys.argv[1]]
  xs = s[0]
  ys = s[1]
  vxs = s[3]
  vys = s[4]
  vzs = s[5] 
  ID = det[7]
  mask = filter(xs, ys, vxs, vys, vzs, ID, 
               Xrange, Yrange, VXrange, VYrange, VZrange)

  if sys.argv[2]=='p':
    if x is None:
      x = det[0][mask]
      y = det[1][mask]
    else:
      x = np.concatenate((x,det[0][mask]))
      y = np.concatenate((y,det[1][mask]))
  if sys.argv[2]=='v':
    if x is None:
      x = det[2][mask]
      y = det[3][mask]
    else:
      x = np.concatenate((x,det[2][mask]))
      y = np.concatenate((y,det[3][mask]))




print "Number of particles: ", len(x)
histy = np.array(y)
histx = np.array(x)
(mux, sigmax) = norm.fit(histx)
(muy, sigmay) = norm.fit(histy)
print 'number of particles', len(x)
print "======================Initial Guess================"
print 'mean X',mux, 'mean Y', muy, 'sig X', sigmax, 'sig Y', sigmay

n = 100
#plt.hist2d(histx, histy, bins=n)
xedges, yedges = np.linspace(min(histx), max(histx), n), np.linspace(min(histy), max(histy), n)
h, x, y = np.histogram2d(histx,histy, (xedges, yedges))
H_rot = h
X, Y = np.meshgrid(x[1:], y[1:])

initial_guess = (np.amax(h), mux,muy,sigmax,sigmay, 0.0, 0.)
#fitting is then done using scipy.curve_fit on the flattened array (numpy.ravel):
popt, pcov = optimize.curve_fit(Gaussian2D, (X, Y), H_rot.ravel(), p0=(initial_guess[0],
               initial_guess[1],initial_guess[2],initial_guess[3],initial_guess[4],initial_guess[5],
                 initial_guess[6]),maxfev=100000)

print "===================Fitting Results================="
print 'Amplitude', popt[0]
print 'X center', popt[1]
print 'Y center', popt[2]
print 'Major FWHM', max([abs(popt[3]), abs(popt[4])])*2.355
print 'Minor FWHM', min([abs(popt[3]), abs(popt[4])])*2.355
print 'Angle', popt[5]
print 'Offset', popt[6]

data_fitted = Gaussian2D((X, Y), *popt)
fig, ax = plt.subplots(1, 1)
cax = ax.imshow(data_fitted.reshape(n-1,n-1), cmap=plt.cm.jet, origin='bottom',
    extent=(X.min(), X.max(), Y.min(), Y.max()))
ax.contour(X, Y, data_fitted.reshape(n-1,n-1), 8, colors='w')
cbar = fig.colorbar(cax)
plt.tight_layout()
plt.show()

