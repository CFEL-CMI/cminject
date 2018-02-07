import sys

from scipy.interpolate import interp1d
import numpy as np

def LoadParticles(filename):
  f = open(filename)
  particles={}
  for i in f:
    token = i.split()
    if len(token)>1:
      if token[0] in particles.keys():
         particles[token[0]].append((float(token[1]), float(token[2]), float(token[3]), float(token[4]), float(token[5]), float(token[6])))
      else:
         particles[token[0]]=[(float(token[1]), float(token[2]), float(token[3]), float(token[4]), float(token[5]), float(token[6]))]

  return particles

def BuildArray(filename, position):
  ys=[]; zs=[]; vxs=[]; vys=[]; vzs=[]  
  p = LoadParticles(filename)
  for i in p.keys():
     x=[]; y=[]; z=[]; vx=[]; vy=[]; vz=[]
     for j in p[i]:
       x.append(j[0]); y.append(j[1]); z.append(j[2])
       vx.append(j[3]); vy.append(j[4]); vz.append(j[5])
     if x[-1]>=position:
       x = np.array(x); y = np.array(y); z=np.array(z) 
       vx = np.array(vx); vy = np.array(vy); vz=np.array(vz) 
       f0 = interp1d(x, y); f1 = interp1d(x, z); f2 = interp1d(x, vx)
       f3 = interp1d(x, vy); f4 = interp1d(x, vz)
       ys.append(float(f0(position)));
       zs.append(float(f1(position))); vxs.append(float(f2(position)));
       vys.append(float(f3(position))); vzs.append(float(f4(position)));
  return ys, zs, vxs, vys, vzs

y, z, vx, vy, vz = BuildArray(sys.argv[1], float(sys.argv[2]))

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
hist = None
if sys.argv[3]=='y':
  hist = np.array(y)
if sys.argv[3]=='z':
  hist = np.array(z)
if sys.argv[3]=='vx':
  hist = np.array(vx)
if sys.argv[3]=='vy':
  hist = np.array(vy)
if sys.argv[3]=='vz':
  hist = np.array(vz)
if hist is not None:
  plt.hist(hist)
  plt.show() 
else:
  print "Please enter the required input"
