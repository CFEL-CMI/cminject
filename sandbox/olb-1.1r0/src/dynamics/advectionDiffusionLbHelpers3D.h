/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008,2017 Orestis Malaspinas, Andrea Parmigiani, Albert Mink
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
#ifndef ADVECTION_DIFFUSION_LB_HELPERS_3D_H
#define ADVECTION_DIFFUSION_LB_HELPERS_3D_H


namespace olb {

/** partial template specialization for D3Q7DescriptorBaseRTLB.
 */
template<typename T>
struct adLbDynamicsHelpers< T,descriptors::D3Q7DescriptorBaseRTLB<T> > {
  static T equilibrium( int iPop, T rho, const T u[3] )
  {
    typedef descriptors::D3Q7DescriptorBaseRTLB<T> L;
    return rho * L::t[iPop] - L::t[iPop];
  }

  // BGK advection diffusion collision step
  static T bgkCollision( T *cell, T rho, const T u[3], T omega )
  {
    const T Cs2 = 0.25;
    const T uSqr = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];

    const T omega_ = (T)1 - omega;
    const T omega_2 = (T)0.5 * omega;
    T rho_ = ( rho - (T)1 );
    cell[0] = omega_ * cell[0] + omega * ( (T)1 - (T)3 * Cs2 ) * rho_;
    const T jx = rho * u[0];
    const T jy = rho * u[1];
    const T jz = rho * u[2];

    rho_ *= Cs2;
    cell[1] = omega_ * cell[1] + omega_2 * ( rho_ - jx );
    cell[2] = omega_ * cell[2] + omega_2 * ( rho_ - jy );
    cell[3] = omega_ * cell[3] + omega_2 * ( rho_ - jz );
    cell[4] = omega_ * cell[4] + omega_2 * ( rho_ + jx );
    cell[5] = omega_ * cell[5] + omega_2 * ( rho_ + jy );
    cell[6] = omega_ * cell[6] + omega_2 * ( rho_ + jz );

    return uSqr;
  }

  /// Paper: Mink et al. 2016 DOI: 10.1016/j.jocs.2016.03.014
  /// omega is expected to be one
  static T sinkCollision( T *cell, T intensity, T omega, T sink )
  {
// equivalent code
//    for ( int iPop = 1; iPop < 7; ++iPop ) {
//      cell[iPop] = intensity*0.125 - sink*( cell[iPop] + 0.125 ) - 0.125;
//    }
//    cell[0] = intensity*0.25 - sink*( cell[0] + 0.125 ) - 0.25;
// equivalent code

    const T omega_ = (T)1 - omega;
    T ti = 0.25;
    cell[0] =  omega_ * ( cell[0] + ti ) + omega * intensity * ti - sink * ( cell[0] + ti ) - ti;

    // fi = (1-w) fi + w fi^{eq} - absorption fi;
    // where fi^{eq} = rho ti and t0 = 1/4, t1=...t6=1/8
    // note: cell[i] = fi - ti
    ti = 0.125;
    cell[1] = omega_ * ( cell[1] + ti ) + omega * intensity * ti - sink * ( cell[1] + ti ) - ti;
    cell[2] = omega_ * ( cell[2] + ti ) + omega * intensity * ti - sink * ( cell[2] + ti ) - ti;
    cell[3] = omega_ * ( cell[3] + ti ) + omega * intensity * ti - sink * ( cell[3] + ti ) - ti;
    cell[4] = omega_ * ( cell[4] + ti ) + omega * intensity * ti - sink * ( cell[4] + ti ) - ti;
    cell[5] = omega_ * ( cell[5] + ti ) + omega * intensity * ti - sink * ( cell[5] + ti ) - ti;
    cell[6] = omega_ * ( cell[6] + ti ) + omega * intensity * ti - sink * ( cell[6] + ti ) - ti;
    return 0.0;
  }
};


template<typename T>
struct adLbDynamicsHelpers<T,descriptors::D3Q7DescriptorBase<T> > {
  /// equilibrium distribution
  static T equilibrium( int iPop, T rho, const T u[3] )
  {
    typedef descriptors::D3Q7DescriptorBase<T> L;
    //T c_u = L::c[iPop][0] * u[0] + L::c[iPop][1] * u[1] + L::c[iPop][2] * u[2];

    return rho * L::t[iPop] - L::t[iPop];
  }

  /// RLB advection diffusion collision step
  static T rlbCollision( T* cell, T rho, const T u[3], T omega )
  {
    typedef descriptors::D3Q7DescriptorBase<T> L;
    const T uSqr = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];

    const T Cs2 = (T)1 / L::invCs2;

    T rho_1 = rho - (T)1;
    cell[0] = ( (T)1 - (T)3 * Cs2 ) * rho_1; //f[0]=(1-3c_s^2)(rho-1)

    const T omega_ = (T)1 - omega;

    const T f1_4 = omega_ * ( cell[1] - cell[4] );
    const T f2_5 = omega_ * ( cell[2] - cell[5] );
    const T f3_6 = omega_ * ( cell[3] - cell[6] );

    rho_1 *= Cs2;

    const T ux_ = omega * rho * u[0];

    cell[1] = (T)0.5 * ( rho_1 + f1_4 - ux_ ); //f[1]=1/2*(c_s^2(rho-1)+(1-omega)*(f[1]-f[4])-omega*rho*u[x])
    cell[4] = (T)0.5 * ( rho_1 - f1_4 + ux_ ); //f[4]=1/2*(c_s^2(rho-1)-(1-omega)*(f[1]-f[4])+omega*rho*u[x])

    const T uy_ = omega * rho * u[1];

    cell[2] = (T)0.5 * ( rho_1 + f2_5 - uy_ ); //f[2]=1/2*(c_s^2(rho-1)+(1-omega)*(f[2]-f[5])-omega*rho*u[y])
    cell[5] = (T)0.5 * ( rho_1 - f2_5 + uy_ ); //f[5]=1/2*(c_s^2(rho-1)-(1-omega)*(f[2]-f[5])+omega*rho*u[y])

    const T uz_ = omega * rho * u[2];

    cell[3] = (T)0.5 * ( rho_1 + f3_6 - uz_ ); //f[3]=1/2*(c_s^2(rho-1)+(1-omega)*(f[3]-f[6])-omega*rho*u[y])
    cell[6] = (T)0.5 * ( rho_1 - f3_6 + uz_ ); //f[6]=1/2*(c_s^2(rho-1)-(1-omega)*(f[3]-f[6])+omega*rho*u[y])

    return uSqr;
  }

  /// BGK advection diffusion collision step
  static T bgkCollision( T *cell, T rho, const T u[3], T omega )
  {
    typedef descriptors::D3Q7DescriptorBase<T> L;

    const T Cs2 = (T)1 / L::invCs2;
    const T uSqr = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];

    const T omega_ = (T)1 - omega;
    const T omega_2 = (T)0.5 * omega;
    T rho_ = ( rho - (T)1 );
    cell[0] = omega_ * cell[0] + omega * ( (T)1 - (T)3 * Cs2 ) * rho_;
    const T jx = rho * u[0];
    const T jy = rho * u[1];
    const T jz = rho * u[2];

    rho_ *= Cs2;
    cell[1] = omega_ * cell[1] + omega_2 * ( rho_ - jx );
    cell[2] = omega_ * cell[2] + omega_2 * ( rho_ - jy );
    cell[3] = omega_ * cell[3] + omega_2 * ( rho_ - jz );
    cell[4] = omega_ * cell[4] + omega_2 * ( rho_ + jx );
    cell[5] = omega_ * cell[5] + omega_2 * ( rho_ + jy );
    cell[6] = omega_ * cell[6] + omega_2 * ( rho_ + jz );

    return uSqr;
  }

};

} // namespace olb

#endif
