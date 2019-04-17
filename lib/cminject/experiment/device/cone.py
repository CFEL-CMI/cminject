class Cone:
   def __init__(self, position, bRadius, tRadius, height):
     self.position = position
     self.bRadius = bRadius
     self.tRadius = tRadius
     self.height = height

   def ParticleInside(self, ParticlePosition):
          z = ParticlePosition[2] - self.position[2]
          permRadius = self.bRadius - (self.bRadius - self.tRadius) * float(z)/self.height
          pointRadius = ((ParticlePosition[0] - self.position[0]) ** 2 + (ParticlePosition[1] - self.position[1]) ** 2)
          if (z <= self.height and z >=0) and (pointRadius <= permRadius**2):
           return True
          else:
           return False

