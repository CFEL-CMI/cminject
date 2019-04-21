from pyquaternion import Quaternion
from scipy import interpolate
from math import pi, cos, sqrt
import numpy as np

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def readPulse(fname):
	""" This function just read the pulse from text file"""
	f = open(fname)
	EE=[]
	Tt=[]
	for i in f:
		token = i.split()
		Tt.append(float(token[0]))
		EE.append(float(token[1]))
	Tt = np.array(Tt)
	EE = np.array(EE)
	return Tt, EE

def timestep(w):
	"""This function calculates the time step 
	needed for the simulations based on the 
	maximum angular velocity"""
	return abs((np.pi/320.)/np.max(w))

def speed(T, I):
	"""This function calculates the intial angular velocity
	based on Temperature and inertia"""
	k = 1.38042e-23
	return 2.8*np.sqrt(2*k*T/I)

class molecule(object):
	"""Class to carry the data of the molecule"""
	def __init__(self, I=None, P=None):
		self.I = I
		self.P = P

	def update_w(self, w, E, a, dt):
		"""Calculate the torque then update 
		the angular velocities"""
		dipole = np.dot(self.P, E)
		torque = np.cross(a*E, dipole)
		dw = self.dw(w, torque)
		w1 = w + dw * dt/2.
		dw = self.dw(w1, torque)
		w = w + dw * dt
		return  w

	def dw(self, w, T):
		"""Calculate the acceleration given the torque"""
		dw = np.zeros((3))
		dw[0] = (T[0] -(self.I[2] - self.I[1])*w[1]*w[2])/self.I[0]
		dw[1] = (T[1] -(self.I[0] - self.I[2])*w[0]*w[2])/self.I[1]
		dw[2] = (T[2] -(self.I[1] - self.I[0])*w[0]*w[1])/self.I[2]
		return dw
			
	def update_q(self, q, w, dt):
		"""the derivatives of the quatrionions given omega"""
		dq = q.derivative(w)
		q = q + dq * dt
		return q/q.norm

Tt, EE = readPulse('laser_pulse.txt')
F = interpolate.interp1d(Tt,EE)
#q = Quaternion.random()  Turn this on if you want to try a new position everytime
q = Quaternion(0.772, 0.032, -0.077, -0.630)
E = np.array([0., 0., 1.0]) # field starts at Z direction
m = molecule(I=np.array([8.83613997e-43, 2.48769068e-39, 2.48799524e-39]), 
             P=np.array([[231.29e-30,   4.23e-30, 0.43e-30],
                         [  4.23e-30, 216.56e-30, 0.29e-30],
                         [  0.43e-30, 0.29e-30, 201.17e-30]]))

am = 4.0 * pi * 1.e16/(2.997929e8) # conversion factor from tera watt/cm2 to V/m
#Principle axis of polarizability for Trp cage, hard coded for bench mark
PA =np.array([-0.966, -0.258, -0.016])
# cos square theta
cosqr=[]
# time
time=[]
#start at 0. sec
t= 0. #1.e-12
# get the speed based on T 0.4 K and I
w = speed(0.4, m.I) 
# get the time step based on w
dt = timestep(w)
# just to plot the pulse
field=[]
# and theta between E and PA
theta=[]
# the Z axis in the lab frame
Z = np.array([0.,0.,1.0])
#start simulations
for i in range(300):
	try:
		a = am*F(t)
	except:
		a= 0.
	t = t + dt
	time.append(t)
	E = np.array([q.rotation_matrix[2][0], q.rotation_matrix[2][1], q.rotation_matrix[2][2]])
	theta.append(cos(angle_between(E, PA))**2)
	cosqr.append(cos(angle_between(E, Z))**2)
	w = m.update_w(w, E, a, dt)
	q = m.update_q(q,w,dt)
	field.append(a)

import matplotlib.pyplot as plt
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
ax1.plot(time,cosqr, label='ZZ')
ax1.set_xlabel('time')
ax1.set_ylabel('cos(theta)')
ax1.legend(loc='upper left', prop={'size': 6})
ax2.plot(time, theta)
ax2.set_xlabel('time')
ax2.set_ylabel('cos(theta) PA')
ax3.plot(time, field)
ax3.set_xlabel('time')
ax3.set_ylabel('I(TW/cm2)')
plt.show()


