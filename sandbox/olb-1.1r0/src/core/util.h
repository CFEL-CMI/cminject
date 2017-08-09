/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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
 * Set of functions commonly used in LB computations
 *  -- header file
 */
#ifndef UTIL_H
#define UTIL_H

#include<sstream>
#include<algorithm>
#include<utilities/vectorHelpers.h>

// patch due to ploblems with older compilers
namespace std {
template<typename T>
std::string to_string(const T &n)
{
  std::ostringstream s;
  s << n;
  return s.str();
}
}

namespace olb {

namespace util {

template<typename T> T norm(const std::vector<T>& a);

template <typename T>
inline int sign(T val)
{
  return (0 < val) - (val < 0);
}

template <typename T>
inline T clamp_low(T val, T a)
{
  return val < a ? a : val;
}

template <typename T>
inline T clamp_high(T val, T b)
{
  return val > b ? b : val;
}

template <typename T>
inline T clamp(T val, T a, T b)
{
  return val < a ? a : (val > b ? b : val);
}

template <typename T>
inline bool aligned_to_x(std::vector<T> vec)
{
  return (vec[0]!=0 and vec[1]==0 and vec[2]==0);
}

template <typename T>
inline bool aligned_to_y(std::vector<T> vec)
{
  return (vec[0]==0 and vec[1]!=0 and vec[2]==0);
}

template <typename T>
inline bool aligned_to_z(std::vector<T> vec)
{
  return (vec[0]==0 and vec[1]==0 and vec[2]!=0);
}

template <typename T>
inline bool aligned_to_grid(std::vector<T> vec)
{
  return (aligned_to_x<T>(vec) or
          aligned_to_y<T>(vec) or
          aligned_to_z<T>(vec));
}





inline bool intersect (
  int x0, int x1, int y0, int y1,
  int x0_, int x1_, int y0_, int y1_,
  int& newX0, int& newX1, int& newY0, int& newY1 )
{
  newX0 = std::max(x0,x0_);
  newY0 = std::max(y0,y0_);

  newX1 = std::min(x1,x1_);
  newY1 = std::min(y1,y1_);

  return newX1>=newX0 && newY1>=newY0;
}

inline bool intersect (
  int x0, int x1, int y0, int y1, int z0, int z1,
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
  int& newX0, int& newX1, int& newY0, int& newY1, int& newZ0, int& newZ1 )
{
  newX0 = std::max(x0,x0_);
  newY0 = std::max(y0,y0_);
  newZ0 = std::max(z0,z0_);

  newX1 = std::min(x1,x1_);
  newY1 = std::min(y1,y1_);
  newZ1 = std::min(z1,z1_);

  return newX1>=newX0 && newY1>=newY0 && newZ1>=newZ0;
}

inline bool contained(int x, int y,
                      int x0, int x1, int y0, int y1)
{
  return x>=x0 && x<=x1 &&
         y>=y0 && y<=y1;
}

inline bool contained(int x, int y, int z,
                      int x0, int x1, int y0, int y1, int z0, int z1)
{
  return x>=x0 && x<=x1 &&
         y>=y0 && y<=y1 &&
         z>=z0 && z<=z1;
}


template<typename T>
T sqr(T arg)
{
  return arg*arg;
}

/// Compute norm square of a d-dimensional vector
template<typename T, int d>
T normSqr(const T u[d])
{
  T uSqr = T();
  for (int iD=0; iD<d; ++iD) {
    uSqr += u[iD]*u[iD];
  }
  return uSqr;
}

template<typename T, int d>
T scalarProduct(const T u1[d], const T u2[d])
{
  T prod = T();
  for (int iD=0; iD<d; ++iD) {
    prod += u1[iD]*u2[iD];
  }
  return prod;
}

template<typename T>
T scalarProduct(const std::vector<T>& u1, const std::vector<T>& u2)
{
  T prod = T();
  if (u1.size() == u2.size()) {
    for (int iD=0; iD<u1.size(); ++iD) {
      prod += u1[iD]*u2[iD];
    }
  }
  return prod;
}

/// Compute number of elements of a symmetric d-dimensional tensor
template <typename Descriptor> struct TensorVal {
  static const int n =
    (Descriptor::d*(Descriptor::d+1))/2; ///< result stored in n
};

/// Compute the opposite of a given direction
template <typename Descriptor> inline int opposite(int iPop)
{
  return Descriptor::opposite[iPop];
}

template <typename Descriptor, int index, int value>
class SubIndex {
private:
  SubIndex()
  {
    for (int iVel=0; iVel<Descriptor::q; ++iVel) {
      if (Descriptor::c[iVel][index]==value) {
        indices.push_back(iVel);
      }
    }
  }

  std::vector<int> indices;

  template <typename Descriptor_, int index_, int value_>
  friend std::vector<int> const& subIndex();
};

template <typename Descriptor, int index, int value>
std::vector<int> const& subIndex()
{
  static SubIndex<Descriptor, index, value> subIndexSingleton;
  return subIndexSingleton.indices;
}

template <typename Descriptor>
int findVelocity(const int v[Descriptor::d])
{
  for (int iPop=0; iPop<Descriptor::q; ++iPop) {
    bool fit = true;
    for (int iD=0; iD<Descriptor::d; ++iD) {
      if (Descriptor::c[iPop][iD] != v[iD]) {
        fit = false;
        break;
      }
    }
    if (fit) {
      return iPop;
    }
  }
  return Descriptor::q;
}

/**
* finds distributions incoming into the wall
* but we want the ones outgoing from the wall,
* therefore we have to take the opposite ones.
*/
template <typename Descriptor, int direction, int orientation>
class SubIndexOutgoing {
private:
  SubIndexOutgoing()   // finds the indexes outgoing from the walls
  {
    indices = util::subIndex<Descriptor,direction,orientation>();

    for (unsigned iPop = 0; iPop < indices.size(); ++iPop) {
      indices[iPop] = util::opposite<Descriptor>(indices[iPop]);
    }

  }

  std::vector<int> indices;

  template <typename Descriptor_, int direction_, int orientation_>
  friend std::vector<int> const& subIndexOutgoing();
};

template <typename Descriptor, int direction, int orientation>
std::vector<int> const& subIndexOutgoing()
{
  static SubIndexOutgoing<Descriptor, direction, orientation> subIndexOutgoingSingleton;
  return subIndexOutgoingSingleton.indices;
}

///finds all rthe remaining indexes of a lattice given some other indexes
template <typename Descriptor>
std::vector<int> remainingIndexes(const std::vector<int> &indices)
{
  std::vector<int> remaining;
  for (int iPop = 0; iPop < Descriptor::q; ++iPop) {
    bool found = false;
    for (unsigned jPop = 0; jPop < indices.size(); ++jPop) {
      if (indices[jPop] == iPop) {
        found = true;
      }
    }
    if (!found) {
      remaining.push_back(iPop);
    }
  }
  return remaining;
}

/// finds the indexes outgoing from a 2D corner
template <typename Descriptor, int xNormal, int yNormal>
class SubIndexOutgoingCorner2D {
private:
  SubIndexOutgoingCorner2D()
  {
    typedef Descriptor L;

    int vect[L::d] = {xNormal, yNormal};
    std::vector<int> knownIndexes;
    knownIndexes.push_back(util::findVelocity<L>(vect));
    vect[0] = xNormal;
    vect[1] = 0;
    knownIndexes.push_back(util::findVelocity<L>(vect));
    vect[0] = 0;
    vect[1] = yNormal;
    knownIndexes.push_back(util::findVelocity<L>(vect));
    vect[0] = 0;
    vect[1] = 0;
    knownIndexes.push_back(util::findVelocity<L>(vect));
    indices = util::remainingIndexes<L>(knownIndexes);
  }

  std::vector<int> indices;

  template <typename Descriptor_, int direction_, int orientation_>
  friend std::vector<int> const& subIndexOutgoingCorner2D();
};

template <typename Descriptor, int xNormal, int yNormal>
std::vector<int> const& subIndexOutgoingCorner2D()
{
  static SubIndexOutgoingCorner2D<Descriptor, xNormal, yNormal> subIndexOutgoingCorner2DSingleton;
  return subIndexOutgoingCorner2DSingleton.indices;
}

/// Util Function for Wall Model of Malaspinas
/// get link with smallest angle to a vector
template <typename T, template <typename U> class DESCRIPTOR>
int get_nearest_link(std::vector<T> vec)
{
  T max=-1;
  int max_index = 0;
  for (int iQ=1; iQ<DESCRIPTOR<T>::q; ++iQ) {
    std::vector<T> c_i(DESCRIPTOR<T>::c[iQ], DESCRIPTOR<T>::c[iQ]+3);
    T tmp = util::scalarProduct<T>(c_i, vec)/util::norm(c_i);
    if (tmp > max) {
      max = tmp;
      max_index = iQ;
    }
  }
  return max_index;
}

namespace tensorIndices2D {
enum { xx=0, xy=1, yy=2 };
}

namespace tensorIndices3D {
enum { xx=0, xy=1, xz=2, yy=3, yz=4, zz=5 };
}

}  // namespace util

}  // namespace olb

#endif
