from molecule import molecule
from field import field
import numpy as np
from math import pi

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::
    """
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    return np.dot(v1_u, v2_u)**2

def integrate_single_particle(j, a, dt):
    Z = np.array([0.,0.,1.0])
    E = np.array([j.Q.rotation_matrix[2][0],
    j.Q.rotation_matrix[2][1], j.Q.rotation_matrix[2][2]])
    costhetaP = angle_between(E, j.PA)
    costhetaI = angle_between(E, Z)
    j.update_w(E, a, dt)
    j.update_q(dt)
    return costhetaP, costhetaI

def integrate_all_particles(n_mol, I, P, field, dt, start, end):
    n = int((end-start)/dt)
    costhetaI = np.zeros(n)
    costhetaP = np.zeros(n)
    t = start
    time = []
    pulse = []
    am = 4.0 * pi * 1.e16/(2.997929e8)
    molecules = []
    for i in range(n_mol):
        molecules.append(molecule(I, P))
    for i in range(n):
        t = t + dt
        a = field.E(t,am)
        time.append(t)
        for j in molecules:
           costhetap, costhetai = integrate_single_particle(j, a, dt)
           costhetaP[i] += costhetap
           costhetaI[i] += costhetai
        pulse.append(a/am)
    costhetaP = costhetaP/n_mol
    costhetaI = costhetaI/n_mol
    return costhetaP, costhetaI, time

def speed(T, I):
        """This function calculates the intial angular velocity
        based on Temperature and inertia"""
        k = 1.38042e-23
        return 2.8*np.sqrt(2*k*T/I)

