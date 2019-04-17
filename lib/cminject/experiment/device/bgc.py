from cminject.experiment.device.cylinder import Cylinder
from cminject.experiment.device.cone import Cone

class BufferGasCell(object):
  """ This class implement the buffer gas cell device"""
  def __init__(self, position, cylinder1=None, cone1=None, cylinder2=None, cone2=None, cylinder3=None):
    self.position = position
    cyL1 = 0.002
    cyR1 = 0.002
    co1bRadius = 0.002
    cotRadius = 0.015
    coL1=0.010
    cyL2 = 0.020
    cyR2 = 0.015
    co2bRadius = 0.001
    cyR3 = 0.001
    cyL3 = 0.001
    self.maxZ = cyL1 + cyL3 + 2*coL1 + cyL2

    if cylinder1 is None:
      self.cylinder1 = Cylinder(cyR1, position, (position[0], position[1], position[2]+cyL1))

    if cone1 is None:
      self.cone1 = Cone((position[0], position[1], position[2]+cyL1), co1bRadius, cotRadius, coL1)

    if cylinder2 is None:
      self.cylinder2 = Cylinder(cyR2, (position[0], position[1], position[2]+cyL1+coL1), 
                                           (position[0], position[1], position[2]+cyL1+coL1+cyL2))
    if cone2 is None:
      self.cone2 = Cone((position[0], position[1], position[2]+cyL1+coL1+cyL2), cotRadius, co2bRadius, coL1)

    if cylinder3 is None:
      self.cylinder3 = Cylinder(cyR3, (position[0], position[1], position[2]+cyL1+(2*coL1)+cyL2), 
                                           (position[0], position[1], position[2]+(cyL1+cyL3)+(2*coL1)+cyL2))
   

  def ParticleInside(self, ParticlePosition):
    if ( self.cylinder1.ParticleInside(ParticlePosition) or
            self.cylinder2.ParticleInside(ParticlePosition) or
              self.cylinder3.ParticleInside(ParticlePosition) or
                        self.cone1.ParticleInside(ParticlePosition) or
                             self.cone2.ParticleInside(ParticlePosition) ):
      return True
    else:
      return False

