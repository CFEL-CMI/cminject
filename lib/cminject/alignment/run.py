from integrate import *
from field import field
I = np.array([8.83613997e-43, 2.48769068e-39, 2.48799524e-39]) # Inertia
P = np.array([[231.29e-30,   4.23e-30, 0.43e-30],  
              [  4.23e-30, 216.56e-30, 0.29e-30],
              [  0.43e-30, 0.29e-30, 201.17e-30]])  # Polarizability Tensor
T = 0.4 #K
n_mol = 100 # number of molecules
startTime = 1.e-10
endTime = 3.e-10
f = field(Type='read', fname='laser_pulse.txt') # read the field in Tera Watt per cm^2
dt = abs((np.pi/320.)/np.max(speed(T, I)))  # calculate the time step
costhetaP, costhetaI, time = integrate_all_particles(n_mol, I, P, f, dt, startTime, endTime) # propagate
import matplotlib.pyplot as plt
plt.plot(time, costhetaP)
plt.show()

