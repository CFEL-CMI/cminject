/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

/** \file
 * Helper functions for the implementation of LB dynamics. This file is all
 * about efficiency. The generic template code is specialized for commonly
 * used Lattices, so that a maximum performance can be taken out of each
 * case.
 */
#ifndef ADVECTION_DIFFUSION_LB_HELPERS_2D_H
#define ADVECTION_DIFFUSION_LB_HELPERS_2D_H


namespace olb {


template<typename T>
struct adLbDynamicsHelpers<T, descriptors::D2Q5DescriptorBase<T> > {
  /// equilibrium distribution
  static T equilibrium( int iPop, T rho, const T u[2] )
  {
    typedef descriptors::D2Q5DescriptorBase<T> L;
    T c_u = L::c[iPop][0] * u[0] + L::c[iPop][1] * u[1];

    return rho * L::t[iPop] * ( ( T )1 + c_u * L::invCs2 ) - L::t[iPop];
  }

  /// RLB advection diffusion collision step
  static T rlbCollision( T* cell,
                         T rho, const T u[2], T omega )
  {
    typedef descriptors::D2Q5DescriptorBase<T> L;
    const T uSqr = u[0] * u[0] + u[1] * u[1];

    const T Cs2 = ( T )1 / L::invCs2;

    T rho_1 = rho - ( T )1;
    cell[0] = ( ( T )1 - ( T )2 * Cs2 ) * rho_1; //f[0]=(1-2c_s^2)(rho-1)

    const T omega_ = ( T )1 - omega;

    const T f1_3 = ( T )0.5 * omega_ * ( cell[1] - cell[3] );
    const T f2_4 = ( T )0.5 * omega_ * ( cell[2] - cell[4] );

    rho_1 *= ( T )0.5 * Cs2;

    const T ux_ = ( T )0.5 * omega * rho * u[0];

    cell[1] = rho_1 + f1_3 - ux_; //f[1]=1/2*(c_s^2(rho-1)+(1-omega)*(f[1]-f[3])-omega*rho*u[x])
    cell[3] = rho_1 - f1_3 + ux_; //f[3]=1/2*(c_s^2(rho-1)-(1-omega)*(f[1]-f[3])+omega*rho*u[x])

    const T uy_ = ( T )0.5 * omega * rho * u[1];

    cell[2] = rho_1 + f2_4 - uy_; //f[2]=1/2*(c_s^2(rho-1)+(1-omega)*(f[2]-f[4])-omega*rho*u[y])
    cell[4] = rho_1 - f2_4 + uy_; //f[4]=1/2*(c_s^2(rho-1)-(1-omega)*(f[2]-f[4])+omega*rho*u[y])

    return uSqr;
  }

  // BGK advection diffusion collision step
  static T bgkCollision( T *cell,
                         T rho, const T u[2], T omega )
  {
    typedef descriptors::D2Q5DescriptorBase<T> L;

    const T Cs2 = ( T )1 / L::invCs2;
    const T uSqr = u[0] * u[0] + u[1] * u[1];

    const T omega_ = ( T )1 - omega;
    const T omega_2 = ( T )0.5 * omega;
    T rho_ = ( rho - ( T )1 );

    cell[0] = omega_ * cell[0] + omega * ( ( T )1 - ( T )2 * Cs2 ) * rho_;

    const T jx = rho * u[0];
    const T jy = rho * u[1];

    rho_ *= Cs2;
    cell[1] = omega_ * cell[1] + omega_2 * ( rho_ - jx );
    cell[2] = omega_ * cell[2] + omega_2 * ( rho_ - jy );
    cell[3] = omega_ * cell[3] + omega_2 * ( rho_ + jx );
    cell[4] = omega_ * cell[4] + omega_2 * ( rho_ + jy );

    return uSqr;
  }

};

} // namespace olb

#endif
