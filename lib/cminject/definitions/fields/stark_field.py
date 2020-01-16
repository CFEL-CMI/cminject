from .fields import RegularGridInterpolationField

class StarkField(RegularGridInterpolationField):
    def __init__(self, filename, ...):
        super().__init__(filename)
        self.memory = {}
        ...

    def calculate_acceleration(self, particle: StarkSphericalParticle, time: float) -> np.array:
        voltage = ...
        interp = self.get_particle_interpolator(particle.stark_file, particle.q, particle.j, ...)
        E = interp(voltage)
        ...
        pass



if __name__ == '__main__':
    from cminject.definitions.particles import StarkParticle
    p = StarkParticle(...)
    f = StarkField(...)

    assert f.get_particle_interpolator(p) = ...
    assert f.calculate_acceleration(p) == something
