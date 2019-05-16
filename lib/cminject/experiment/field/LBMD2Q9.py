import numpy as np
import matplotlib


class LBM(object):
  """ This class implement Lattice Boltzmann method in 2D based on D2Q9 lattice,
      I used the following link as a guidance: 
      http://benchpress.readthedocs.io/autodoc_benchmarks/lattice_boltzmann_D2Q9.html"""

  def __init__(self, Lx, Ly, dx, uMax, dt, time, Re, nu, als, OutR, viz=None):
    self.nx = int(Lx/dx)   # number of LB cells in the x-direction 
    self.ny = int(Ly/dx)   # number of LB cells in the y-direction
    self.uMax = uMax  # maximum inlet speed in lattice units
    self.n = int(time/dt) # number of iterations
    self.Re = Re  # Renold number
    self.omega = omega  = 1. / (3*nu+1./2.)  # Relaxation parameter
    self.als = np.array(als)/dx
    self.OutR = OutR/dx
    self.t = [4/9., 1/9.,1/9.,1/9.,1/9., 1/36.,1/36.,1/36.,1/36.] # D2Q9 lattice constants     
    self.cx = [  0,   1,  0, -1,  0,    1,  -1,  -1,   1] # 9 lattice velocities vx
    self.cy = [  0,   0,  1,  0, -1,    1,   1,  -1,  -1] # 9 lattice velocities vy
    self.opp = [  0,   3,  4,  1,  2,    7,   8,   5,   6]  # for bounce back boundary condition
    t_3d    = np.asarray(self.t)[:, np.newaxis, np.newaxis] # make t (9,1,1) array to generate the population from it
    self.fin = t_3d.repeat(self.nx, axis = 1).repeat(self.ny, axis = 2) # Initial condition: (rho=1, u=0) ==> fIn[i] = t[i]
    self.PrepareGeometry()
    # Loop over time iterations 
    for itr in range(self.n):
      print("Iteration number", itr)
      
      self.ComputeRhoU()

      self.SetPoisseuilleProfile()

      self.SetDensityInlet()

      self.SetZeroGradientOutlet()

      self.Equilibrium()

      self.SetPopulationInletOutlet()

      self.Stream()
      
#      self.BounceBack(itr/50, 15, 2)
      if viz and not itr % 10:
        self.visualize(viz, itr) 
  
  
  def visualize(self, viz, itr):
        viz['axes1'].clear()
        viz['axes2'].clear()
        viz['axes3'].clear()
        im=viz['axes1'].imshow(self.fin.sum(axis=0).T.copy(), interpolation='nearest')
        im=viz['axes2'].imshow(self.ux.T.copy(), interpolation='nearest')
#        if itr==100:
#          fig= viz['axes2'].get_figure()
#          fig.colorbar(im, orientation='vertical')
        viz['axes3'].imshow(self.uy.T.copy(), interpolation='nearest')
        viz['canvas'].show()


  def PrepareGeometry(self):
    obst = np.zeros((self.nx,self.ny))
    obst[:,0] = obst[:,-1] = 1
    obst[self.nx-int(self.nx/20)-1:-1, 0:int((self.ny/2)-self.OutR)]=1
    obst[self.nx-int(self.nx/20)-1:-1, int((self.ny/2)+self.OutR):-1]=1
    for i in self.als:
      OutR = int(i[1])
      cord = int(i[0])
#      print(i[0], i[1], self.nx, self.ny)
      obst[cord:cord+1, 0:int((self.ny/2)-OutR)]=1
      obst[cord:cord+1, int((self.ny/2)+OutR):-1]=1
    self.bbRegion = np.asarray(obst, dtype = np.bool)

  def ComputeRhoU(self):
    """Calculates U and Rho given the population"""
    self.rho = np.sum(self.fin, axis=0) # the density for each cell is the summation of population in the cell
    print("Total Density = ", np.sum(self.fin)/(self.nx*self.ny))
    cx_3d   = np.asarray(self.cx)[:, np.newaxis, np.newaxis] # put cx in 3d format to multiply it with fin
    cy_3d   = np.asarray(self.cy)[:, np.newaxis, np.newaxis] # put cy in 3d format to multiply it with fin
    self.ux = np.sum(cx_3d * self.fin, axis = 0) / self.rho # the velocity is the lattice velocities multiplied by \
    self.uy = np.sum(cy_3d * self.fin, axis = 0) / self.rho # population devided by density

  def SetPoisseuilleProfile(self):
    col = np.arange(1, self.ny - 1)
    L = self.ny-2.0 # the first and last cell don't have a flow
    y = col-0.5  # take the coordinates as the center of the cell
    self.ux[0, 1:self.ny-1] = 4 * self.uMax / (L ** 2) * (y * L - y ** 2) # Poisseuille equation
    self.uy[0, 1:self.ny-1] = 0 # no flow in the y direction
  
  def SetDensityInlet(self):
    t1 = self.fin[[0, 2, 4], 0, 1:self.ny-1].sum(axis = 0) # The population of the left cells is unknown
    t2 = 2 * self.fin[[3, 6, 7], 0, 1:self.ny-1].sum(axis = 0) # however, it could be obtained by the known in the inlet velocity
    self.rho[0, 1:self.ny-1] = 1 / (1 - self.ux[0, 1:self.ny-1]) * (t1 + t2)

  def SetZeroGradientOutlet(self):
    """ Outlet: Zero gradient on rho/ux"""
    self.rho[self.nx - 1, 1:self.ny-1] = 1
    #self.rho[self.nx - 1, 1:self.ny-1] = 4.0 / 3 * self.rho[self.nx - 2, 1:self.ny-1] - \
    #                  1.0 / 3 * self.rho[self.nx - 3, 1:self.ny-1]
    #self.uy[self.nx - 1, 1:self.ny-1] = 0
    #self.ux[self.nx - 1, 1:self.ny-1] = 4.0 / 3 * self.ux[self.nx - 2, 1:self.ny-1] - \
    #                  1.0 / 3 * self.ux[self.nx - 3, 1:self.ny-1]
  
  def Equilibrium(self):
    self.fEq = np.zeros((9, self.nx, self.ny))
    self.fout = np.zeros((9, self.nx, self.ny))
    for i in range(0, 9):
      cu = 3 * (self.cx[i] * self.ux + self.cy[i] * self.uy)
      self.fEq[i] = self.rho * self.t[i] * (1 + cu + 0.5 * cu ** 2 - 1.5 * (self.ux ** 2 + self.uy ** 2))
      self.fout[i] = self.fin[i] - self.omega * (self.fin[i] - self.fEq[i])

  def SetPopulationInletOutlet(self):
    for i in range(0, 9):
      # Left boundary:
      self.fout[i, 0, 1:self.ny-1] = self.fEq[i,0,1:self.ny-1] + 18 * self.t[i] * self.cx[i] * self.cy[i] * \
             (self.fin[7,0,1:self.ny-1] - self.fin[6,0,1:self.ny-1] - self.fEq[7,0,1:self.ny-1] + self.fEq[6,0,1:self.ny-1])
      # Right boundary:
      self.fout[i,self.nx-1,1:self.ny-1] = self.fEq[i,self.nx-1,1:self.ny-1] + 18 * self.t[i] * self.cx[i] * self.cy[i] * \
             (self.fin[5,self.nx-1,1:self.ny-1] - self.fin[8,self.nx-1,1:self.ny-1] - \
              self.fEq[5,self.nx-1,1:self.ny-1] + self.fEq[8,self.nx-1,1:self.ny-1])
      # Bounce Back Wall
      self.fout[i,self.bbRegion] = self.fin[self.opp[i],self.bbRegion]
      

  def BounceBackParticle(self, obstx=0, obsty=0, obstr=0):
    y,x = np.meshgrid(np.arange(self.ny),np.arange(self.nx))
    obst = (x-obstx)**2 + (y-obsty)**2 <= obstr**2
 #  obst = np.zeros((self.nx,self.ny))
    obst[:,0] = obst[:,self.ny-1] = 1
    obst[self.nx-10:self.nx-1, 0:15]=1
    obst[self.nx/2:1+self.nx/2, 0:9]=1
    obst[self.nx/2:1+self.nx/2, self.ny-9:self.ny-1]=1
    obst[self.nx-10:self.nx-1, self.ny-15:self.ny-1]=1
    self.bbRegion = np.asarray(obst, dtype = np.bool)
    

  def Stream(self):
    for i in range(9):
      self.fin[i] = np.roll( np.roll(self.fout[i], self.cx[i], axis=0), self.cy[i], axis=1)

"""
viz = setup_viz()  
#LB = LBM( 1, 1, 0.02, 0.1, 0.1, 1, 10, 0.001, viz)
import threading
thread = threading.Thread(target = LBM, args=( 2, 1, 0.01, 0.025, 1, 10000, 10, 0.01, viz))
thread.start()
viz['tk'].mainloop() 
"""  
