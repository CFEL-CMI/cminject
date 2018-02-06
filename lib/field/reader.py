import vtk
import numpy as np

def ReadVTK(fileName):
  a = vtk.vtkXMLImageDataReader()
  a.SetFileName(fileName)
  a.Update()
  bound = a.GetOutput().GetBounds()
  # (0.75, 53.25, -0.25, 30.25, -0.25, 30.25)
  xbound = (bound[0], bound[1])
  ybound = (bound[2], bound[3])
  zbound = (bound[4], bound[5])
  xspace = a.GetOutput().GetSpacing()[0]
  yspace = a.GetOutput().GetSpacing()[1]
  zspace = a.GetOutput().GetSpacing()[2]
  nx, ny, nz = a.GetOutput().GetDimensions()
  c=0
  f = open('vtktotext', 'w+')
  X=[]; Y=[]; Z=[]; VX=[]; VY=[]; VZ=[]; P=[];
  for k in range(nz):
   for j in range(ny):
     for i in range(nx):
       x = i * xspace + xbound[0]
       y = j * yspace + ybound[0]
       z = k * zspace + zbound[0]
       vx=a.GetOutput().GetPointData().GetArray(0).GetTuple3(c)[0]
       vy=a.GetOutput().GetPointData().GetArray(0).GetTuple3(c)[1]
       vz=a.GetOutput().GetPointData().GetArray(0).GetTuple3(c)[2]
       p =a.GetOutput().GetPointData().GetArray(1).GetTuple1(c)
       X.append(x); Y.append(y); Z.append(z); VX.append(vx); VY.append(vy); VZ.append(vz); P.append(p)
  
  return X, Y, Z, VX, VY, VZ, P

def ReadFromFile(filename):
  f = open(filename)
  x=[]; y=[]; z=[]
  vx=[]; vy=[]; vz=[]
  p=[]
  for i in f:
    token = i.split()
    x.append(float(token[0]))
    y.append(float(token[1]))
    z.append(float(token[2]))

    if token[3]=='NaN':
      vx.append(0.0)
    else:
      vx.append(float(token[3]))

    if token[4]=='NaN':
      vy.append(0.0)
    else:
      vy.append(float(token[4]))

    if token[5]=='NaN':
      vz.append(0.0)
    else:
      vz.append(float(token[5]))

    if token[6]=='NaN':
      p.append(0.0)
    else:
      p.append(float(token[6]))

  x=sorted(set(x))
  y=sorted(set(y))
  z=sorted(set(z))
  n_x = len(x)
  n_y = len(y)
  n_z = len(z)
  c=0
  Vx = np.zeros((n_x, n_y, n_z))
  Vy = np.zeros((n_x, n_y, n_z))
  Vz = np.zeros((n_x, n_y, n_z))
  P = np.zeros((n_x, n_y, n_z))
  for i in range(n_z):
    for j in range(n_y):
      for k in range(n_x):
        Vx[k,j,i] = vx[c]
        Vy[k,j,i] = vy[c]
        Vz[k,j,i] = vz[c]
        P[k,j,i] = p[c]
        c+=1
  return x, y, z, Vx, Vy, Vz, P
