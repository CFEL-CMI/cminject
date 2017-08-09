/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, Orestis Malaspinas and Jonas Latt
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

#ifndef INAMURO_NEWTON_RAPHSON_DYNAMICS_HH
#define INAMURO_NEWTON_RAPHSON_DYNAMICS_HH

#include "inamuroNewtonRaphsonDynamics.h"
#include "dynamics/latticeDescriptors.h"
#include "core/util.h"
#include "dynamics/lbHelpers.h"
#include <cmath>


namespace olb {

using namespace descriptors;

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
InamuroNewtonRaphsonDynamics<T,Lattice,Dynamics,direction,orientation>::InamuroNewtonRaphsonDynamics (
  T omega, Momenta<T,Lattice>& momenta )
  : BasicDynamics<T,Lattice>(momenta),
    _boundaryDynamics(omega, momenta),
    clout(std::cout,"InamuroNewtonRaphsonDynamics")
{
  _xi[0] = (T)1;
  for (int iDim = 1; iDim < Lattice<T>::d; ++iDim) {
    _xi[iDim] = T();
  }
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
InamuroNewtonRaphsonDynamics<T,Lattice,Dynamics,direction,orientation>*
InamuroNewtonRaphsonDynamics<T,Lattice, Dynamics, direction, orientation>::clone() const
{
  return new InamuroNewtonRaphsonDynamics<T,Lattice,Dynamics,direction,orientation>(*this);
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
T InamuroNewtonRaphsonDynamics<T,Lattice, Dynamics, direction, orientation>::
computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const
{
  return _boundaryDynamics.computeEquilibrium(iPop, rho, u, uSqr);
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
void InamuroNewtonRaphsonDynamics<T,Lattice,Dynamics,direction,orientation>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  typedef Lattice<T> L;

  T rho, u[L::d];
  this->_momenta.computeRhoU(cell,rho,u);

  std::vector<int> missingIndexes = util::subIndexOutgoing<L,direction,orientation>();
  std::vector<int> knownIndexes;
  bool test[L::q];
  for (int iPop = 0; iPop < L::q; ++iPop) {
    test[iPop] = true;
  }

  for (unsigned iPop = 0; iPop < missingIndexes.size(); ++iPop) {
    test[missingIndexes[iPop]] = false;
  }
  for (int iPop = 0; iPop < L::q; ++iPop) {
    if (test[iPop]) {
      knownIndexes.push_back(iPop);
    }
  }

  T approxMomentum[L::d];

  computeApproxMomentum(approxMomentum,cell,rho,u,_xi,knownIndexes,missingIndexes);

  T error = computeError(rho, u,approxMomentum);
  int counter = 0;

  while (error > 1.0e-15) {
    ++counter;

    T gradError[L::d], gradGradError[L::d][L::d];
    computeGradGradError(gradGradError,gradError,rho,u,_xi,approxMomentum,missingIndexes);

    bool everythingWentFine = newtonRaphson(_xi,gradError,gradGradError);
    if ((counter > 100000) || everythingWentFine == false) {
      // if we need more that 100000 iterations or
      // if we have a problem with the inversion of the
      // jacobian matrix, we stop the program and
      // print this error message on the screen.
      clout << "Failed to converge...." << std::endl;
      clout << "Error = " << error << std::endl;
      clout << "u = (" << rho*u[0];
      for (int iD=1; iD<Lattice<T>::d; ++iD) {
        clout << ", " << rho*u[iD];
      }
      clout << ")" << std::endl;
      clout << "uApprox = (" << approxMomentum[0];
      for (int iD=1; iD<Lattice<T>::d; ++iD) {
        clout << ", " << approxMomentum[iD];
      }
      clout << ")" << std::endl;
      clout << "xi = (" << _xi[0];
      for (int iD=1; iD<Lattice<T>::d; ++iD) {
        clout << ", " << _xi[iD];
      }
      clout << ")" << std::endl;

      exit(1);
    }

    computeApproxMomentum(approxMomentum,cell,rho,u,_xi,knownIndexes,missingIndexes);
    error = computeError(rho, u,approxMomentum);

  }

  T uCs[L::d];
  int counterDim = 0;
  for (int iDim = 0; iDim < L::d; ++iDim) {
    if (direction == iDim) {
      ++counterDim;
      uCs[iDim] = u[iDim];
    } else {
      uCs[iDim] = u[iDim] + _xi[iDim+1-counterDim];
    }
  }

  T uCsSqr = util::normSqr<T,L::d>(uCs);

  for (unsigned iPop = 0; iPop < missingIndexes.size(); ++iPop) {
    cell[missingIndexes[iPop]] = computeEquilibrium(missingIndexes[iPop],_xi[0],uCs,uCsSqr);
  }

  _boundaryDynamics.collide(cell, statistics);
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
void InamuroNewtonRaphsonDynamics<T,Lattice,Dynamics,direction,orientation>::staticCollide (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  typedef Lattice<T> L;

  T rho = this->_momenta.computeRho(cell);

  std::vector<int> missingIndexes = util::subIndexOutgoing<L,direction,orientation>();
  std::vector<int> knownIndexes;
  bool test[L::q];
  for (int iPop = 0; iPop < L::q; ++iPop) {
    test[iPop] = true;
  }

  for (unsigned iPop = 0; iPop < missingIndexes.size(); ++iPop) {
    test[missingIndexes[iPop]] = false;
  }
  for (int iPop = 0; iPop < L::q; ++iPop) {
    if (test[iPop]) {
      knownIndexes.push_back(iPop);
    }
  }

  T approxMomentum[L::d];

  computeApproxMomentum(approxMomentum,cell,rho,u,_xi,knownIndexes,missingIndexes);

  T error = computeError(rho, u,approxMomentum);
  int counter = 0;

  while (error > 1.0e-15) {
    ++counter;

    T gradError[L::d], gradGradError[L::d][L::d];
    computeGradGradError(gradGradError,gradError,rho,u,_xi,approxMomentum,missingIndexes);

    bool everythingWentFine = newtonRaphson(_xi,gradError,gradGradError);
    if ((counter > 100000) || everythingWentFine == false) {
      // if we need more that 100000 iterations or
      // if we have a problem with the inversion of the
      // jacobian matrix, we stop the program and
      // print this error message on the screen.
      clout << "Failed to converge...." << std::endl;
      clout << "Error = " << error << std::endl;
      clout << "u = (" << rho*u[0];
      for (int iD=1; iD<Lattice<T>::d; ++iD) {
        clout << ", " << rho*u[iD];
      }
      clout << ")" << std::endl;
      clout << "uApprox = (" << approxMomentum[0];
      for (int iD=1; iD<Lattice<T>::d; ++iD) {
        clout << ", " << approxMomentum[iD];
      }
      clout << ")" << std::endl;
      clout << "xi = (" << _xi[0];
      for (int iD=1; iD<Lattice<T>::d; ++iD) {
        clout << ", " << _xi[iD];
      }
      clout << ")" << std::endl;

      exit(1);
    }

    computeApproxMomentum(approxMomentum,cell,rho,u,_xi,knownIndexes,missingIndexes);
    error = computeError(rho, u,approxMomentum);

  }

  T uCs[L::d];
  int counterDim = 0;
  for (int iDim = 0; iDim < L::d; ++iDim) {
    if (direction == iDim) {
      ++counterDim;
      uCs[iDim] = u[iDim];
    } else {
      uCs[iDim] = u[iDim] + _xi[iDim+1-counterDim];
    }
  }

  T uCsSqr = util::normSqr<T,L::d>(uCs);

  for (unsigned iPop = 0; iPop < missingIndexes.size(); ++iPop) {
    cell[missingIndexes[iPop]] = computeEquilibrium(missingIndexes[iPop],_xi[0],uCs,uCsSqr);
  }

  _boundaryDynamics.staticCollide(cell, u, statistics);

}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
T InamuroNewtonRaphsonDynamics<T,Lattice,Dynamics,direction,orientation>::getOmega() const
{
  return _boundaryDynamics.getOmega();
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
void InamuroNewtonRaphsonDynamics<T,Lattice,Dynamics,direction,orientation>::setOmega(T omega)
{
  _boundaryDynamics.setOmega(omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
T InamuroNewtonRaphsonDynamics<T,Lattice,Dynamics,direction,orientation>::getParameter(int whichParameter) const
{
  return _boundaryDynamics.getParameter(whichParameter);
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
void InamuroNewtonRaphsonDynamics<T,Lattice,Dynamics,direction,orientation>::setParameter(int whichParameter, T value)
{
  _boundaryDynamics.setParameter(whichParameter, value);
}


template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
void InamuroNewtonRaphsonDynamics<T,Lattice,Dynamics,direction,orientation>::
computeApproxMomentum(T approxMomentum[Lattice<T>::d],const Cell<T,Lattice> &cell,
                      const T &rho, const T u[Lattice<T>::d], const T xi[Lattice<T>::d],
                      const std::vector<int> knownIndexes,const std::vector<int> missingIndexes)
{
  typedef Lattice<T> L;

  T csVel[L::d];
  int counter = 0;
  for (int iDim = 0; iDim < L::d; ++iDim) {
    if (direction == iDim) {
      ++counter;
      csVel[iDim] = u[iDim];
    } else {
      csVel[iDim] = u[iDim] + xi[iDim+1-counter];
    }
  }

  T csVelSqr = util::normSqr<T,L::d>(csVel);

  for (int iDim = 0; iDim < L::d; ++iDim) {
    approxMomentum[iDim] = T();
    for (unsigned iPop = 0; iPop < knownIndexes.size(); ++iPop) {
      approxMomentum[iDim] += L::c[knownIndexes[iPop]][iDim] *
                              cell[knownIndexes[iPop]];
    }
    for (unsigned iPop = 0; iPop < missingIndexes.size(); ++iPop) {
      approxMomentum[iDim] += L::c[missingIndexes[iPop]][iDim] *
                              computeEquilibrium(missingIndexes[iPop],xi[0],csVel,csVelSqr);
    }
  }
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
T InamuroNewtonRaphsonDynamics<T,Lattice,Dynamics,direction,orientation>::
computeError(const T &rho, const T u[Lattice<T>::d], const T approxMomentum[Lattice<T>::d])
{
  typedef Lattice<T> L;

  T err = T();
  for (int iDim = 0; iDim < L::d; ++iDim) {
    err += (rho * u[iDim]-approxMomentum[iDim]) * (rho * u[iDim]-approxMomentum[iDim]);
  }
  return sqrt(err);
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
void InamuroNewtonRaphsonDynamics<T,Lattice,Dynamics,direction,orientation>::computeGradGradError(
  T gradGradError[Lattice<T>::d][Lattice<T>::d],T gradError[Lattice<T>::d],
  const T &rho, const T u[Lattice<T>::d],const T xi[Lattice<T>::d],
  const T approxMomentum[Lattice<T>::d],
  const std::vector<int> missingIndexes)
{
  typedef Lattice<T> L;

  T csVel[L::d];
  int counter = 0;
  for (int iDim = 0; iDim < L::d; ++iDim) {
    if (direction == iDim) {
      ++counter;
      csVel[iDim] = u[iDim];
    } else {
      csVel[iDim] = u[iDim] + xi[iDim+1-counter];
    }
  }
  T csVelSqr = util::normSqr<T,L::d>(csVel);

  std::vector<T> df[L::d];

  for (int iXi = 0; iXi < L::d; ++iXi) {
    for (unsigned iPop = 0; iPop < missingIndexes.size(); ++iPop) {
      df[iXi].push_back(T());
    }
  }

  // all the null terms are allready contained in df (no need for else in the ifs)

  for (unsigned iPop = 0; iPop < missingIndexes.size(); ++iPop) {
    T cu = T();
    for (int iDim = 0; iDim < L::d; ++iDim) {
      cu += L::c[missingIndexes[iPop]][iDim] * csVel[iDim];
    }
    df[0][iPop] = L::t[missingIndexes[iPop]]
                  * ((T)1+L::invCs2*cu
                     + 0.5 * L::invCs2 * L::invCs2 * cu * cu
                     - 0.5 * L::invCs2 * csVelSqr);
  }

  counter = 0;
  for (int iDim = 0; iDim < L::d; ++iDim) {
    if (direction == iDim) {
      ++counter;
    } else {
      for (unsigned iPop = 0; iPop < missingIndexes.size(); ++iPop) {
        T temp = T();
        for (int jDim = 0; jDim < L::d; ++jDim) {
          temp += L::c[missingIndexes[iPop]][jDim] * csVel[jDim];
        }

        df[iDim+1-counter][iPop] = xi[0]*L::t[missingIndexes[iPop]] *
                                   (L::invCs2*L::c[missingIndexes[iPop]][iDim]
                                    + L::invCs2*L::invCs2*L::c[missingIndexes[iPop]][iDim]*temp
                                    - L::invCs2*csVel[iDim]);
      }
    }
  }

  std::vector<T> ddf[L::d][L::d];

  for (int iAlpha = 0; iAlpha < L::d; ++iAlpha) {
    for (int iBeta = 0; iBeta < L::d; ++iBeta) {
      for (unsigned iPop = 0; iPop < missingIndexes.size(); ++iPop) {
        ddf[iAlpha][iBeta].push_back(T());
      }
    }
  }

  // ddf contains allready all the zero terms

  counter = 0;
  for (int iDim = 0; iDim < L::d; ++iDim) {
    if (direction == iDim) {
      ++counter;
    } else {
      for (unsigned iPop = 0; iPop < missingIndexes.size(); ++iPop) {
        T temp = T();
        for (int jDim = 0; jDim < L::d; ++jDim) {
          temp += L::c[missingIndexes[iPop]][jDim] * csVel[jDim];
        }

        T d_rho_sa = L::t[missingIndexes[iPop]] *
                     (L::invCs2*L::c[missingIndexes[iPop]][iDim]
                      + L::invCs2*L::invCs2*L::c[missingIndexes[iPop]][iDim]*temp
                      - L::invCs2*csVel[iDim]);

        ddf[iDim+1-counter][0][iPop] = d_rho_sa;
        ddf[0][iDim+1-counter][iPop] = d_rho_sa;
      }
    }
  }

  for (int iAlpha = 1; iAlpha < L::d; ++iAlpha) {
    for (int iBeta = 1; iBeta < L::d; ++iBeta) {
      for (unsigned iPop = 0; iPop < missingIndexes.size(); ++iPop) {
        ddf[iAlpha][iBeta][iPop] = L::t[missingIndexes[iPop]]*xi[0] *
                                   L::invCs2*L::invCs2*L::c[missingIndexes[iPop]][iAlpha]*L::c[missingIndexes[iPop]][iBeta];
        if (iAlpha == iBeta) {
          ddf[iAlpha][iBeta][iPop] -= L::t[missingIndexes[iPop]]*xi[0] * L::invCs2;
        }
      }
    }
  }

  // computation of the vector gradError
  counter = 0;
  for (int iDim = 0; iDim < L::d; ++iDim) {
    T du[L::d];
    gradError[iDim] = T();
    for (int jDim = 0; jDim < L::d; ++jDim) {
      du[jDim] = T();
      for (unsigned iPop = 0; iPop < missingIndexes.size(); ++iPop) {
        du[jDim] += L::c[missingIndexes[iPop]][jDim] * df[iDim][iPop];
      }
      gradError[iDim] += (approxMomentum[jDim]- rho * u[jDim]) * du[jDim];
    }
    gradError[iDim] *= (T)2;
  }

  // computation of the matrix gradGradError

  for (int iAlpha = 0; iAlpha < L::d; ++iAlpha) {
    for (int iBeta = 0; iBeta < L::d; ++iBeta) {
      gradGradError[iAlpha][iBeta] = T();

      T duAlpha[L::d], duBeta[L::d], ddu[L::d];
      for (int iDim = 0; iDim < L::d; ++iDim) {
        duAlpha[iDim] = T();
        duBeta[iDim] = T();
        ddu[iDim] = T();
        for (unsigned iPop = 0; iPop < missingIndexes.size(); ++iPop) {
          duAlpha[iDim] += L::c[missingIndexes[iPop]][iDim]
                           * df[iAlpha][iPop];

          duBeta[iDim] += L::c[missingIndexes[iPop]][iDim]
                          * df[iBeta][iPop];

          ddu[iDim] += L::c[missingIndexes[iPop]][iDim]
                       * ddf[iAlpha][iBeta][iPop];
        }
        gradGradError[iAlpha][iBeta] += duAlpha[iDim] * duBeta[iDim]
                                        + (approxMomentum[iDim]-rho * u[iDim]) * ddu[iDim];
      }
      gradGradError[iAlpha][iBeta] *= (T)2;
      // gradGradError = 2*sum_alpha(du_alpha/dxi_beta*du_alpha/dxi_gamma +
      //                  approxMomentum_alpha-u_alpha*d^2u_alpha/(dxi_gamma*dxi_beta))
    }
  }
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
bool InamuroNewtonRaphsonDynamics<T,Lattice,Dynamics,direction,orientation>::
newtonRaphson(T xi[Lattice<T>::d],
              const T gradError[Lattice<T>::d],
              const T gradGradError[Lattice<T>::d][Lattice<T>::d])
{
  typedef Lattice<T> L;

  T invGradGradError[L::d][L::d];
  bool inversion = invert(gradGradError,invGradGradError);

  for (int iXi = 0; iXi < L::d; ++iXi) {
    T correction = T();
    for (int iAlpha = 0; iAlpha < L::d; ++iAlpha) {
      correction += invGradGradError[iXi][iAlpha] * gradError[iAlpha];
    }

    xi[iXi] -= correction;
  }

  return inversion;
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
bool InamuroNewtonRaphsonDynamics<T,Lattice,Dynamics,direction,orientation>::
invert(const T a[2][2],T invA[2][2])
{
  T detA = a[0][0]*a[1][1] - a[0][1]*a[1][0];

  if (fabs(detA) < 1.0e-13) {
    clout << "error detA too small! = " << detA << std::endl;
    for (int iAlpha = 0; iAlpha < 2; ++iAlpha) {
      for (int iBeta = 0; iBeta < 2; ++iBeta) {
        clout << a[iAlpha][iBeta] << "\t";
      }
      clout << std::endl;
    }
    return false;
  } else {
    invA[0][0] = a[1][1];
    invA[1][1] = a[0][0];

    invA[0][1] = -a[1][0];
    invA[1][0] = -a[0][1];

    for (int iPop = 0; iPop < 2; ++iPop) {
      for (int jPop = 0; jPop < 2; ++jPop) {
        invA[iPop][jPop] /= detA;
      }
    }
    return true;
  }
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
bool InamuroNewtonRaphsonDynamics<T,Lattice,Dynamics,direction,orientation>::invert(const T a[3][3],T invA[3][3])
{
  T detA = a[0][0]*a[1][1]*a[2][2] + a[1][0]*a[2][1]*a[0][2] + a[2][0]*a[0][1]*a[1][2]
           - a[0][0]*a[2][1]*a[1][2] - a[2][0]*a[1][1]*a[0][2] - a[1][0]*a[0][1]*a[2][2];


  if (fabs(detA) < 1.0e-13) {
    clout << "Error: detA too small! = " << detA << std::endl;
    for (int iAlpha = 0; iAlpha < 3; ++iAlpha) {
      for (int iBeta = 0; iBeta < 3; ++iBeta) {
        clout << a[iAlpha][iBeta] << "\t";
      }
      clout << std::endl;
    }
    return false;
  } else {
    invA[0][0] = a[1][1]*a[2][2]-a[1][2]*a[2][1];
    invA[0][1] = a[0][2]*a[2][1]-a[0][1]*a[2][2];
    invA[0][2] = a[0][1]*a[1][2]-a[0][2]*a[1][1];

    invA[1][0] = a[1][2]*a[2][0]-a[1][0]*a[2][2];
    invA[1][1] = a[0][0]*a[2][2]-a[0][2]*a[2][0];
    invA[1][2] = a[0][2]*a[1][0]-a[0][0]*a[1][2];

    invA[2][0] = a[1][0]*a[2][1]-a[1][1]*a[2][0];
    invA[2][1] = a[0][1]*a[2][0]-a[0][0]*a[2][1];
    invA[2][2] = a[0][0]*a[1][1]-a[0][1]*a[1][0];

    for (int iPop = 0; iPop < 3; ++iPop) {
      for (int jPop = 0; jPop < 3; ++jPop) {
        invA[iPop][jPop] /= detA;
      }
    }
    return true;
  }
}

}  // namespace olb


#endif
