from cminject.definitions.fields.regular_grid_interpolation_field import RegularGridInterpolationField
from cminject.definitions.particles import Molecule
import h5py as hp
import numpy as np
import scipy as sc
from typing import Tuple
import sys

class B_StarkField(RegularGridInterpolationField):

    def __init__(self, f_filename, e_filename):
        super().__init__(filename=f_filename)
        self.energy_file = e_filename
        #self.memory = {}
        #...

    def calculate_acceleration(self, particle: Molecule, time: float) -> np.array:

        #voltage = ...
        #interp = self.get_particle_interpolator(particle.stark_file, particle.q, particle.j, ...)
        #E = interp(voltage)
        #...
        #pass
        # construct a path of the quantum state of the particle:
        _, grad = self.get_local_properties(particle_position)
        return -1*(1/particle.mass)*self.energy_interpolate(particle_position)*grad #I'm assuming mass is already in kilogram

    def energy_interpolate(self, particle_position: np.array) -> float:
        voltage, _ = self.get_local_properties(particle_position)
        print(particle.q_n['J'])
        path = '_' + str(particle.q_n['J']) + '/_' + str(particle.q_n['Ka']) + '/_' + str(particle.q_n['Kc']) + '/_' + str(particle.q_n['M']) + '/_' + str(particle.q_n['Isomer'])
        with hp.File(self.energy_file, 'r') as stark:
            dc = stark.get(path + '/dcfield')
            dc = np.asarray(dc)
            enr = stark.get(path + '/dcstarkenergy')
            enr = np.asarray(enr)
            print("we are reading data successfully")
            sys.exit()
            mu_eff = np.gradient(enr, dc)
            mueff_interp = sc.interpolate.interp1d(dc,mueff) # should I do it like this? Wathc interp object

        return mueff_interp(voltage)
        # calculate gradient

    def field_interpolate(self, particle_position: np.array) -> np.array:
        return self._interpolator(particle_position)

    def get_local_properties(self, particle_position: np.array) -> Tuple[float, np.array]:

        data = self.field_interpolate(particle_position)
        grad = data[:3]
        norm = data[3]          # we are only considering the norm of the field for now
        return norm, grad


if __name__ == '__main__':

    e_file = '../../../../../wd_nr_static_symmetrized'
    f_file = '../../../../../test'
    d = {'J': 2, 'Ka': 2, 'Kc': 2, 'M': 2, 'Isomer':0}
    pos = np.array([1,2,3])
    p = Molecule(34., d)
    f = B_StarkField(f_file, e_file)
    f.energy_interpolate(p)

    #assert f.get_particle_interpolator(p) = ...
    #assert f.calculate_acceleration(p) == something
