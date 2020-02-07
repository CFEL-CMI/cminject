from cminject.definitions.fields.regular_grid_interpolation_field import RegularGridInterpolationField
from cminject.definitions.particles import Molecule
import h5py as hp
import numpy as np
import scipy as sc
from typing import Tuple
import sys

class B_StarkField(RegularGridInterpolationField):
    """
    Units of the variables:
    Estark: Joul. The field in cmistark: volt/cm
    Gradient of the electric field: volt/cm^2
    Mass: atomic unit
    Resulting acceleration is in cm/s^2
    """

    def __init__(self, f_filename: str):
        super().__init__(filename=f_filename)
        self.memory = {}

    def calculate_acceleration(self, particle: Molecule, time: float) -> np.array:

        au_kg = (1.660538782e-27)  #converting mass from atomic unit to kilogram
        icm_j = 1.98630e-23 # convert form inverse centimiters to joul
        icm_j_cm = 1.98630e-19 # joul=kg. m^(2)/s^(2). this is to express joul in centimeter squared
        mass = mass*au_kg

        voltage, grad = self.get_local_properties(particle.position[0:2])
        return np.concatenate([-1*(1/particle.mass)*icm_j_cm*self.energy_interpolate(particle.energy_filename, particle.q_n, voltage)*grad[0:2], np.array([0.])]) #I'm assuming mass is already in kilogram

    def energy_interpolate(self, energy_filename, particle_qn, voltage) -> float:
        # retrieve the values of the dictionary
        J = particle_qn['J']
        Ka = particle_qn['Ka']
        Kc = particle_qn['Kc']
        M = particle_qn['M']
        Iso = particle_qn['Isomer']
        path = '_' + str(J) + '/_' + str(Ka) + '/_' + str(Kc) + '/_' + str(M) + '/_' + str(Iso)
        key = energy_filename + path
        if key in self.memory.keys():
            mueff_interp = self.memory[key]
        else:
            with hp.File(energy_filename, 'r') as stark:
                # note that in cmi-stark the data is not stored as np arrays, so I need to edit this code a bit
                dc = stark.get(path + '/dcfield')
                dc = np.asarray(dc)
                dc = dc[0] # needed because data is stored as VLarray in cmi-stark output. If it was stored as numpy arrays this step is not needed.
                enr = stark.get(path + '/dcstarkenergy')
                enr = np.asarray(enr)
                enr = enr[0]
                mu_eff = np.gradient(enr, dc)
                mueff_interp = sc.interpolate.interp1d(dc,mu_eff, fill_value = 'extrapolate') # should I do it like this? Wathc interp object
            self.memory[key] = mueff_interp
        return mueff_interp(voltage)

    def field_interpolate(self, particle_position: np.array) -> np.array:
        return self._interpolator(particle_position)

    def get_local_properties(self, particle_position) -> Tuple[float, np.array]:
        data = self.field_interpolate(particle_position)
        grad = data[0:3]
        norm = data[3]          # we are only considering the norm of the field for now
        return norm, grad


if __name__ == '__main__':

    e_file = '../../../../../wd_nr_static_symmetrized'
    f_file = '../../../../../test'
    f2_file = '../../../../../5mmCone_10sccm.h5'

    d = {'J': 2, 'Ka': 2, 'Kc': 2, 'M': 2, 'Isomer':0}
    identifier= 0
    position = np.zeros((6,))
    start_time = 0.5
    pos = np.array([1,2,3])
    p = Molecule(34., d, e_file, identifier, start_time, position)
    f = B_StarkField(f_file)
    m = f.calculate_acceleration(p, start_time)
    #indic = [x[0] for x in np.nditer(a, order = 'F')]
    #print(indic)
    #assert f.get_particle_interpolator(p) = ...
    #assert f.calculate_acceleration(p) == something
