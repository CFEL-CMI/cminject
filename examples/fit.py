import sys
import h5py
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from scipy.stats import norm

f = h5py.File(sys.argv[1], 'r')
x=[]; y=[]; z=[]; vx=[]; vy=[]; vz=[]
d = f['detectors']
det = d[sys.argv[2]]
if sys.argv[3]=='p':
  x= det[0]
  y= det[1]
if sys.argv[3]=='v':
  x= det[2]
  y= det[3]

def Gaussian2D(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
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

histy = np.array(y)
histx = np.array(x)
(mux, sigmax) = norm.fit(histx)
(muy, sigmay) = norm.fit(histy)
print 'number of particles', len(x)
print "Initial Guess:   ",
print 'mean X',mux, 'mean Y', muy, 'sig X', sigmax, 'sig Y', sigmay
#print max(vz)
#plt.hist2d(x, y, bins=100)
n = 150
xedges, yedges = np.linspace(min(histx), max(histx), n), np.linspace(min(histy), max(histy), n)
h, x, y = np.histogram2d(histx,histy, (xedges, yedges))
H_rot = h

X, Y = np.meshgrid(x[1:], y[1:])
initial_guess = (np.amax(h), mux,muy,sigmax,sigmay,0.,0.)
#fitting is then done using scipy.curve_fit on the flattened array (numpy.ravel):
popt, pcov = optimize.curve_fit(Gaussian2D, (X, Y), H_rot.ravel(), p0=(initial_guess[0],
               initial_guess[1],initial_guess[2],initial_guess[3],initial_guess[4],initial_guess[5],
                 initial_guess[6]),maxfev=100000)


data_fitted = Gaussian2D((X, Y), *popt)
fig, ax = plt.subplots(1, 1)
cax = ax.imshow(data_fitted.reshape(n-1,n-1), cmap=plt.cm.jet, origin='bottom',
    extent=(X.min(), X.max(), Y.min(), Y.max()))
ax.contour(X, Y, data_fitted.reshape(n-1,n-1), 8, colors='w')
cbar = fig.colorbar(cax)
#cbar.ax.set_yticklabels(['0', str(popt[0])])  # vertically oriented colorbar

print 'Amplitude', popt[0]
print 'X center', popt[1]
print 'Y center', popt[2]
print 'X sigma', popt[3]
print 'Y sigma', popt[4]
print 'Angle', popt[5]
print 'Offset', popt[6]

plt.show()
