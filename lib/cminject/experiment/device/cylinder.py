class Cylinder:
   def __init__(self, radius, start, end):
     self.radius = radius
     self.start = start
     self.end = end

   def ParticleInside(self, ParticlePosition):
     if (((ParticlePosition[0]-self.start[0])**2 + (ParticlePosition[1]-self.start[1])**2 <= self.radius**2)
                              and (ParticlePosition[2] >= self.start[2] and ParticlePosition[2] <= self.end[2])):
        return True
     else:
        return False

