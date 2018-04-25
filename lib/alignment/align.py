import sys
import numpy as np
from math import sin, cos
from numpy.linalg import inv
from scipy.integrate import ode

class atom:
  def __init__(self, coord, mass):
    self.coord = coord
    self.mass = mass

def read_pdb_xyz(pdb_name):
    """
    read xyz from pdb
    return:
    [[x1 y1 z1]
    [x2 y2 z2]
    [.. .. ..]
    [xn yn zn]]
    """
    atoms = []
    pdb_file = open(pdb_name, 'r')
    for line in pdb_file:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            if line[13] == "C":               # these atomic masses in atomic units
                 atoms.append(atom([x, y, z], 12.0107)) # if the pdb file contains other atoms then 
            if line[13] == "O":
                 atoms.append(atom([x, y, z], 15.9994)) # more elements need to be added
            if line[13] == "N":
                 atoms.append(atom([x, y, z], 14.0067))
            if line[13] == "H":
                 atoms.append(atom([x, y, z], 1.0079))

    pdb_file.close()
    return atoms


def matriz_inercia(atoms):
	'''
	DESCRIPTION

	The method calculates the mass center, the inertia tensor and the eigenvalues and eigenvectors
	for a given selection. Mostly taken from inertia_tensor.py
	'''

	totmass = 0.0
	x,y,z = 0,0,0
	for a in atoms:
		m = a.mass
		x += a.coord[0]*m
		y += a.coord[1]*m
		z += a.coord[2]*m
		totmass += m
	global cM
	cM = np.array([x/totmass, y/totmass, z/totmass])

	I = []
	for index in range(9):
		I.append(0)

	for a in atoms:
		temp_x, temp_y, temp_z = a.coord[0], a.coord[1], a.coord[2]
		temp_x -= x
		temp_y -= y
		temp_z -= z

		I[0] += a.mass * (temp_y**2 + temp_z**2)
		I[1] -= a.mass * temp_x * temp_y
		I[2] -= a.mass * temp_x * temp_z
		I[3] -= a.mass * temp_x * temp_y
		I[4] += a.mass * (temp_x**2 + temp_z**2)
		I[5] -= a.mass * temp_y * temp_z
		I[6] -= a.mass * temp_x * temp_z
		I[7] -= a.mass * temp_y * temp_z
		I[8] += a.mass * (temp_x**2 + temp_y**2)
	tensor = np.array([(I[0:3]), (I[3:6]), (I[6:9])]) * 1.66053886 * 1.e-27 * 1.e-20 #SI Units
        return tensor

def rotationMatrix(theta):
   tx,ty,tz = theta
   Rx = np.array([[1,0,0], [0, cos(tx), -sin(tx)], [0, sin(tx), cos(tx)]])
   Ry = np.array([[cos(ty), 0, -sin(ty)], [0, 1, 0], [sin(ty), 0, cos(ty)]])
   Rz = np.array([[cos(tz), -sin(tz), 0], [sin(tz), cos(tz), 0], [0,0,1]])
   return np.dot(Rx, np.dot(Ry, Rz))

def rotateTensor(angles, tensor):
   R = rotationMatrix(angles)
   return np.dot(R, np.dot(tensor, R.transpose()))

atoms = read_pdb_xyz(sys.argv[1])
I = matriz_inercia(atoms)
conversion =  (0.529177)**3 * (1.e-24  * 1.e-15) / 8.988
alpha = np.array([[1406.079,   -4.035, 5.050],
                   [-4.035, 1293.306, 3.469],
                   [5.050, 3.469, 1213.409]]) * conversion #1.11265e-16 * (5.29177e-9)**3

E = np.array([10., 0, 0])
dipole = np.dot(alpha, E)
torque = np.cross(dipole, E)
omega = np.array([0, 0, 0])
theta = (0., 0., 0.)
alpha1 = rotateTensor(theta, alpha)


dt=1.e-9
def get_v_and_a(t, p_and_v, alpha, E, I):
  theta = p_and_v[:3]
  omega = p_and_v[3:]
  alpha = rotateTensor(theta, alpha)
  dipole = np.dot(alpha, E)
  torque = np.cross(dipole, E)
  I = rotateTensor(theta, I)
  Iw = I.dot(omega)
  OIw = np.cross(omega, Iw)
  M = torque - OIw
  invI = inv(I)
  a = invI.dot(M)
  return  np.concatenate((omega, a)) 
  
integral = ode( get_v_and_a )
integral.set_integrator('lsoda') #, method='BDF',with_jacobian=False,atol=1e-8,rtol=1e-4,first_step=1e-5,nsteps=10000)
integral.set_initial_value( (np.array([0,0,0,0,0,0])), 0. ).set_f_params(alpha, E, I)
print("Calculate angular trajectories")

while integral.successful() and integral.t < 0.01:

  integral.integrate(integral.t + dt)

print integral.y[0], integral.y[1], integral.y[2]
print integral.y[3], integral.y[4], integral.y[5]
