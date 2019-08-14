cimport numpy as np
import numpy as np
from libc.math cimport floor
from cython cimport boundscheck, wraparound, nonecheck, cdivision

np.import_array()

@cdivision(True)
@boundscheck(False)
@wraparound(False)
@nonecheck(False)
cpdef np.ndarray[np.float64_t, ndim=1] interp2D(np.float_t[:,:,::1] v, np.float_t x, np.float_t y,
                                                int V, int X, int Y):
    cdef:
        int i, x0, x1, y0, y1, dim
        np.float_t xd, yd, a0, a1
        np.float_t *v_a

    cdef np.ndarray[np.float64_t, ndim=1] result = np.full(V, np.nan, dtype=np.float64)

    x0 = <int>floor(x)
    x1 = x0 + 1
    y0 = <int>floor(y)
    y1 = y0 + 1

    xd = (x1-x)/(x1-x0)
    yd = (y1-y)/(y1-y0)

    if x0 >= 0 and y0 >= 0 and x1 <= X and y1 <= Y:
        for ai in range(V):
            v_a = &v[ai,0,0]

            a0 = xd*v_a[Y*x0+y0] + (1-xd)*v_a[Y*x1+y0]
            a1 = xd*v_a[Y*x0+y1] + (1-xd)*v_a[Y*x1+y1]

            result[ai] = yd*a0 + (1-yd)*a1
    else:
        raise ValueError("At least one of the coordinate components is outside of the grid!")
    return result


@cdivision(True)
@boundscheck(False)
@wraparound(False)
@nonecheck(False)
cpdef np.ndarray[np.float64_t, ndim=1] interp3D(np.float_t[:,:,:,::1] v, np.float_t x, np.float_t y, np.float_t z,
                                                int V, int X, int Y, int Z):
    cdef:
        int i, x0, x1, y0, y1, z0, z1, dim
        np.float_t xd, yd, zd, c00, c01, c10, c11, c0, c1, c
        np.float_t *v_c

    cdef np.ndarray[np.float64_t, ndim=1] result = np.full(V, np.nan, dtype=np.float64)

    x0 = <int>floor(x)
    x1 = x0 + 1
    y0 = <int>floor(y)
    y1 = y0 + 1
    z0 = <int>floor(z)
    z1 = z0 + 1

    xd = (x-x0)/(x1-x0)
    yd = (y-y0)/(y1-y0)
    zd = (z-z0)/(z1-z0)

    if x0 >= 0 and y0 >= 0 and z0 >= 0 and x1 <= X and y1 <= Y and z1 <= Z:
        for ci in range(V):
            v_c = &v[ci,0,0,0]

            c00 = v_c[Y*Z*x0+Z*y0+z0]*(1-xd) + v_c[Y*Z*x1+Z*y0+z0]*xd
            c01 = v_c[Y*Z*x0+Z*y0+z1]*(1-xd) + v_c[Y*Z*x1+Z*y0+z1]*xd
            c10 = v_c[Y*Z*x0+Z*y1+z0]*(1-xd) + v_c[Y*Z*x1+Z*y1+z0]*xd
            c11 = v_c[Y*Z*x0+Z*y1+z1]*(1-xd) + v_c[Y*Z*x1+Z*y1+z1]*xd

            c0 = c00*(1-yd) + c10*yd
            c1 = c01*(1-yd) + c11*yd

            result[ci] = c0*(1-zd) + c1*zd
    else:
        raise ValueError("At least one of the coordinate components is outside of the grid!")
    return result



### Local Variables:
### mode: Python
### fill-column: 100
### truncate-lines: t
### End:
