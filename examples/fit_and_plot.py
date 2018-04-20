import sys
f = open(sys.argv[1])
x=[]; y=[]; z=[]; vx=[]; vy=[]; vz=[]
for i in f:
  token = i.split()
  if len(token)>1:
    x.append(float(token[0])); y.append(float(token[1]));
    vx.append(float(token[2])); vy.append(float(token[3])); vz.append(float(token[4]))

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
#hist = np.sqrt(np.array(y)**2 + np.array(z)**2)
#hist = np.sqrt(np.array(vx)**2+np.array(vx)**2)
histy = np.array(y)
histx = np.array(x)
(mux, sigmax) = norm.fit(histx)
(muy, sigmay) = norm.fit(histy)
print 'number of particles', len(x), 'mean',mux, 'standard deviation', (sigmax+sigmay)/2
#print max(vz)

n, bins, patches = plt.hist(histx, 100, normed=True) #,  100, normed=0.1, facecolor='green', alpha=0.70)
y = mlab.normpdf( bins, mux, sigmax)


#y = y/38.
l = plt.plot(bins, y, 'r--', linewidth=2)
sigmas= sigmax*1000000

plt.xlabel('X')
plt.ylabel('Number of Particles')
plt.title(r'$\mathrm{Histogram\ of\ Particle\ Distribution\ at\ 50sccm:}\ \mu=%.3f,\ \sigma=%.3f$ 10^-4' %(mux, sigmas))
plt.grid(True)


plt.show()

