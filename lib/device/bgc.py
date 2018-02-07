class Cone:
   def __init__(self, position, bRadius, tRadius, height):
     self.position = position
     self.bRadius = bRadius
     self.tRadius = tRadius
     self.height = height

   def ConeInside(self, ParticlePosition):
          z = ParticlePosition[0] - self.position[0]
          permRadius = self.bRadius - (self.bRadius - self.tRadius) * float(z)/self.height
          pointRadius = ((ParticlePosition[1] - self.position[1]) ** 2 + (ParticlePosition[2] - self.position[2]) ** 2)
          if (z <= self.height and z >=0) and (pointRadius <= permRadius**2):
#           print "InSide Cone"
           return True
          else:
           return False
     
class Cylinder:
   def __init__(self, radius, start, end):
     self.radius = radius
     self.start = start
     self.end = end
   
   def CylinderInside(self, ParticlePosition):
     if (((ParticlePosition[1]-self.start[1])**2 + (ParticlePosition[2]-self.start[2])**2 <= self.radius**2) 
                              and (ParticlePosition[0] >= self.start[0] and ParticlePosition[0] <= self.end[0])):
#        print "InSide Cyliner"
        return True
     else:
        return False

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
    if cylinder1 is None:
      self.cylinder1 = Cylinder(cyR1, position, (position[0]+cyL1, position[1], position[2]))

    if cone1 is None:
      self.cone1 = Cone((position[0]+cyL1, position[1], position[2]), co1bRadius, cotRadius, coL1)

    if cylinder2 is None:
      self.cylinder2 = Cylinder(cyR2, (position[0]+cyL1+coL1, position[1], position[2]), 
                                           (position[0]+cyL1+coL1+cyL2, position[1], position[2]))
    if cone2 is None:
      self.cone2 = Cone((position[0]+cyL1+coL1+cyL2, position[1], position[2]), cotRadius, co2bRadius, coL1)

    if cylinder3 is None:
      self.cylinder3 = Cylinder(cyR3, (position[0]+cyL1+(2*coL1)+cyL2, position[1], position[2]), 
                                           (position[0]+(cyL1+cyL3)+(2*coL1)+cyL2, position[1], position[2]))
   

  def BufferGasCellInside(self, ParticlePosition):
    if ( self.cylinder1.CylinderInside(ParticlePosition) or
            self.cylinder2.CylinderInside(ParticlePosition) or
              self.cylinder3.CylinderInside(ParticlePosition) or
                        self.cone1.ConeInside(ParticlePosition) or
                             self.cone2.ConeInside(ParticlePosition)):
      return True
    else:
      return False

