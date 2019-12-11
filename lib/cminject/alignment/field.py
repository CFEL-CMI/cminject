from scipy import interpolate
import numpy as np

class field:
    def __init__(self, Type='read', fname=None):
        if Type=='read':
            f = open(fname)
            EE=[]
            Tt=[]
            for i in f:
                token = i.split()
                Tt.append(float(token[0]))
                EE.append(float(token[1]))
            Tt = np.array(Tt)
            EE = np.array(EE)
            self.F = interpolate.interp1d(Tt,EE)
    
    def E(self, t, a):
        try:
            an = a * self.F(t)
        except:
            an = 0.
        return an

