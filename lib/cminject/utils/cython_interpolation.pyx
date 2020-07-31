from libc.math cimport floor
from cython cimport boundscheck, wraparound, nonecheck, cdivision

@cdivision(True)
@boundscheck(False)
@wraparound(False)
@nonecheck(False)
cpdef int interp2d(double[:,:,::1] v, double x, double y,
                   int nd, int nx, int ny, double[::1] out) except -1:
    """
    Interpolates a n-dimensional vector field bilinearly based on a 2D regular data grid.

    :param v: The data array. (nx, ny, nd)-shaped, where nx and ny are the number of points in each
        dimension of the regular grid, and nd is the number of dimensions that are interpolated
        (the size of the output vector). MUST be C-contiguous.
    :param x: The x position to interpolate at.
    :param y: The y position to interpolate at.
    :param nd: As described in v.
    :param nx: As described in v.
    :param ny: As described in v.
    :param out: The interpolated output vector, which is an interpolated (nd,)-shaped np.array.

    :return: Nothing
    """
    cdef:
        int x0, x1, y0, y1
        double xd, yd, a0, a1
        double *v_a

    # Get the x0,y0 and x1,y1 corner indices of the rectangle that the point (x,y) lies within
    x0 = <int>floor(x)
    x1 = x0 + 1
    y0 = <int>floor(y)
    y1 = y0 + 1

    # These are the linear interpolation factors, xd for x and yd for y
    xd = (x1-x)/(x1-x0)
    yd = (y1-y)/(y1-y0)

    # All the calculated corners must be within the grid bounds...
    if x0 >= 0 and y0 >= 0 and x1 < nx and y1 < ny:
        # For each output dimension ai:
        for ai in range(nd):
            # Get the slice of the data array for this dimension
            v_a = &v[ai,0,0]

            # Interpolate bilinearly, getting the 1D index from the calculated 2D indices
            a0 = xd*v_a[ny*x0+y0] + (1-xd)*v_a[ny*x1+y0]
            a1 = xd*v_a[ny*x0+y1] + (1-xd)*v_a[ny*x1+y1]

            # Linearly interpolate and write the value to the result vector
            out[ai] = yd*a0 + (1-yd)*a1
    # ...otherwise raise an error
    else:
        raise ValueError("At least one of the coordinate components is outside of the grid!")


@cdivision(True)
@boundscheck(False)
@wraparound(False)
@nonecheck(False)
cpdef int interp3d(double[:,:,:,::1] v, double x, double y, double z,
                   int nd, int nx, int ny, int nz, double[::1] out) except -1:
    """
    Interpolates a n-dimensional vector field trilinearly based on a 3D regular data grid.

    :param v: The data array. (nx, ny, nz, nd)-shaped, where nx, ny, nz are the number of points in
      each dimension of the regular grid, and nd is the number of dimensions that are interpolated
      (the size of the output vector). MUST be C-contiguous.
    :param x: The x position to interpolate at.
    :param y: The y position to interpolate at.
    :param z: The z position to interpolate at.
    :param nd: As described in v.
    :param nx: As described in v.
    :param ny: As described in v.
    :param nz: As described in v.
    :param out: The interpolated output vector, which is an interpolated (nd,)-shaped np.array.
    :return: Nothing
    """
    cdef:
        int x0, x1, y0, y1, z0, z1
        double xd, yd, zd, c00, c01, c10, c11, c0, c1
        double *v_c

    # Get the indices for the 6 corners of the cube that (x,y,z) lies within
    x0 = <int>floor(x)
    x1 = x0 + 1
    y0 = <int>floor(y)
    y1 = y0 + 1
    z0 = <int>floor(z)
    z1 = z0 + 1

    # Calculate the interpolation factors for all 3 dimensions
    xd = (x-x0)/(x1-x0)
    yd = (y-y0)/(y1-y0)
    zd = (z-z0)/(z1-z0)

    # All the calculated corners must be within the grid bounds...
    if x0 >= 0 and y0 >= 0 and z0 >= 0 and x1 < nx and y1 < ny and z1 < nz:
        # For each output dimension ci:
        for ci in range(nd):
            # Get the slice of the data array for this dimension
            v_c = &v[ci,0,0,0]

            # Interpolate trilinearly, getting the 1D index from the calculated 3D indices
            c00 = v_c[ny*nz*x0+nz*y0+z0]*(1-xd) + v_c[ny*nz*x1+nz*y0+z0]*xd
            c01 = v_c[ny*nz*x0+nz*y0+z1]*(1-xd) + v_c[ny*nz*x1+nz*y0+z1]*xd
            c10 = v_c[ny*nz*x0+nz*y1+z0]*(1-xd) + v_c[ny*nz*x1+nz*y1+z0]*xd
            c11 = v_c[ny*nz*x0+nz*y1+z1]*(1-xd) + v_c[ny*nz*x1+nz*y1+z1]*xd

            # Bilinearly interpolate
            c0 = c00*(1-yd) + c10*yd
            c1 = c01*(1-yd) + c11*yd

            # Linearly interpolate and write the value to the result vector
            out[ci] = c0*(1-zd) + c1*zd
    # ...otherwise raise an error
    else:
        raise ValueError("At least one of the coordinate components is outside of the grid!")



### Local Variables:
### mode: Python
### fill-column: 100
### truncate-lines: t
### End:
