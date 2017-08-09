/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Asher Zarth, Mathias J. Krause, Albert Mink
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
 * efficient implementation of a vector class
 */

#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <cstring>
#include "util.h"
#include "math.h"

typedef double S;

// Use BaseType double as the default BaseType.
#ifndef BaseType
///TODO
//check OLD_PRECONDITION
//change this typedef, double for testing only
//#define BaseType S
#define BaseType S
#endif



namespace olb {




/// descriptor base
template <unsigned Size = 3>
class DescriptorBase {

protected:
  /// number of dimensions
  enum { size = Size };
};

/** Efficient implementation of a vector class.
 * Static size.
 * Default value of vector element is T(0)
 */
template <typename T, unsigned Size>
class Vector : public DescriptorBase<Size> {
public:
  /// number of dimensions
  T data[Size];
public:
  /// Constructor from array
  Vector(const T ar[Size]);
  /// Constructor from std::vector<T>
  Vector(const std::vector<T>& v);
  /// Construct with all entries as explicit scalar
  Vector(const T& s = T());
  /// Construct 2D, 2 explicit arguments
  Vector(const T& a, const T& b);
  /// Construct 3D, 3 explicit arguments
  Vector(const T& a, const T& b, const T& c);
  /// copy constructor
  Vector(const Vector<T,Size>& v);
  /// copy assignment operator
  Vector& operator = (const Vector<T,Size>& v);
  /// destructor
  ~Vector() = default;
  /// element access
  T& operator[](unsigned n);
  /// element access (read only)
  const T& operator[](const unsigned n) const;
  /// cumulative add vector
  Vector& operator += (const Vector<T,Size>& v);
  /// cumulatively subtract vector
  Vector& operator -= (const Vector<T,Size>& v);
  /// add scalar to each entry
  Vector& operator += (const T& s);
  /// subtract scalar from each entry
  Vector& operator -= (const T& s);
  /// multiply by scalar
  Vector& operator *= (const T& s);
  /// equal operator returns true if all components match
  bool operator == (const Vector<T, Size>& v);
  /// unequal operator returns true if at least one components differs
  bool operator != (const Vector<T, Size>& v);
  /// get l2 norm of vector
  T norm();
  /// get normalized vector
  void normalize(const T& scale = T(1));
  /// check if vector has zero length
  bool closeToZero();
  //operator BaseType();
  /// convert to std::vector
  std::vector<T> toStdVector() const;
};



template <typename T, unsigned Size>
inline Vector<T, Size>::Vector(const T ar[Size]) : DescriptorBase<Size>()
{
  for (unsigned i = 0; i < Size; ++i) {
    data[i] = ar[i];
  }
}

template <typename T, unsigned Size>
inline Vector<T, Size>::Vector(const std::vector<T>& v) : DescriptorBase<Size>()
{
  //OLB_PRECONDITION(v.size() == Size);
  for (unsigned i = 0; i < Size; ++i) {
    data[i] = (T)v[i];
  }
}

template <typename T, unsigned Size>
inline Vector<T, Size>::Vector(const T& s) : DescriptorBase<Size>()
{
  for (unsigned i = 0; i < Size; ++i) {
    data[i] = s;
  }
}

template <typename T, unsigned Size>
inline Vector<T, Size>::Vector(const T& a, const T& b) : DescriptorBase<Size>()
{
  //OLB_PRECONDITION(this->size == 2S);
  data[0] = a;
  data[1] = b;
}

template <typename T, unsigned Size>
inline Vector<T, Size>::Vector(const T& a, const T& b, const T& c) : DescriptorBase<Size>()
{
  //OLB_PRECONDITION(this->size == 3);
  data[0] = a;
  data[1] = b;
  data[2] = c;
}

template <typename T, unsigned Size>
inline Vector<T, Size>::Vector(const Vector<T,Size>& v) : DescriptorBase<Size>()
{
  std::memcpy( data, v.data, Size*sizeof(T));
}

template <typename T, unsigned Size>
inline Vector<T, Size>& Vector<T, Size>::operator= (const Vector<T,Size>& v)
{
  for (unsigned i = 0; i < Size; ++i) {
    data[i] = v[i];
  }
  return *this;
}

template <typename T, unsigned Size>
inline Vector<T, Size>& Vector<T, Size>::operator+= (const Vector<T,Size>& v)
{
  for (unsigned i = 0; i < Size; ++i) {
    data[i] += v[i];
  }
  return *this;
}

template <typename T, unsigned Size>
inline Vector<T, Size>& Vector<T, Size>::operator+= (const T& s)
{
  for (unsigned i = 0; i < Size; ++i) {
    data[i] += s;
  }
  return *this;
}

template <typename T, unsigned Size>
inline Vector<T, Size>& Vector<T, Size>::operator-= (const Vector<T,Size>& v)
{
  for (unsigned i = 0; i < Size; ++i) {
    data[i] -= v[i];
  }
  return *this;
}

template <typename T, unsigned Size>
inline Vector<T, Size>& Vector<T, Size>::operator*= (const T& s)
{
  for (unsigned i = 0; i < Size; ++i) {
    data[i] *= s;
  }
  return *this;
}

template <typename T, unsigned Size>
inline Vector<T, Size>& Vector<T, Size>::operator-= (const T& s)
{
  for (unsigned i = 0; i < Size; ++i) {
    data[i] -= s;
  }
  return *this;
}

template <typename T, unsigned Size>
inline bool Vector<T, Size>::operator== (const Vector<T, Size>& v)
{
  bool theSame = true;
  for (unsigned i = 0; i < Size; ++i) {
    theSame &= data[i] == v.data[i];
  }
  return theSame;
}

template <typename T, unsigned Size>
inline bool Vector<T, Size>::operator!= (const Vector<T, Size>& v)
{
  return !(*this == v);
}

template <typename T, unsigned Size>
inline T& Vector<T, Size>::operator[](unsigned n)

{
  return data[n];
}

template <typename T, unsigned Size>
inline const T& Vector<T, Size>::operator[](const unsigned n) const
{
  return data[n];
}

template <typename T, unsigned Size>
inline bool Vector<T, Size>::closeToZero()
{
  bool zero = true;
  T EPSILON = std::numeric_limits<T>::epsilon();

  for (unsigned i = 0; i < Size; ++i) {
    if (fabs(data[i]) > EPSILON) {
      zero = false;
    }
  }
  return zero;
}

template <typename T, unsigned Size>
inline T Vector<T, Size>::norm()
{
  T v(0);
  for (unsigned iD = 0; iD < Size; ++iD) {
    v += data[iD]*data[iD];
  }
  v = sqrt(v);
  return v;
}

template <typename T, unsigned Size>
inline void Vector<T, Size>::normalize(const T& scale)
{
  T invScale = scale/this->norm();
  //assert(scale > 0);
  for (unsigned int iDim = 0; iDim < Size; ++iDim) {
    data[iDim] *= invScale;
  }
}

template <typename T, unsigned Size>
inline std::vector<T> Vector<T, Size>::toStdVector() const
{
  return std::vector<T> (data, data+Size);
}

template <class T, unsigned Size>
inline std::ostream& operator << (std::ostream& os, const Vector<T,Size>& o)
{
  if (Size != 0) {
    os << "[";
    for (unsigned i = 0; i < Size-1; ++i) {
      os << o[i] << " ";
    }
    os << o[Size-1]<<"]";
  } else {
    os << "[empty]";
  }
  return os;
}

template <typename T>
inline Vector<T, 3> crossProduct3D(const Vector<T, 3>& a, const Vector<T,3>& b)
{
  Vector<T, 3> v(a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]);
  return v;
}

template <typename T, unsigned Size>
inline T norm(const Vector<T, Size>& a)
{
  T v(0);
  for (unsigned iD = 0; iD < Size; ++iD) {
    v += a[iD]*a[iD];
  }
  v = sqrt(v);
  return v;
}

template <typename T, unsigned Size>
inline void normalize(Vector<T, Size>& a)
{
  T scale = norm(a);
  //assert(scale > 0);
  for (unsigned iDim = 0; iDim < Size; ++iDim) {
    a[iDim] /= scale;
  }
}




// indirect helper functions, providing general arithmetic
// overloading of operators for interaction with BaseType and itself
// addition, subtraction, multiplication and division

// PLUS
template <class T, unsigned Size>
inline Vector<T, Size> operator+ (const BaseType& a, const Vector<T, Size>& b)
{
  return Vector<T, Size>(b)+=a;
}

template <class T, unsigned Size>
inline Vector<T, Size> operator+ (const Vector<T, Size>& a, const BaseType& b)
{
  return Vector<T, Size>(a)+=b;
}

template <class T, unsigned Size>
inline Vector<T, Size> operator+ (const Vector<T, Size>& a, const Vector<T, Size>& b)
{
  return Vector<T, Size>(a)+=b;
}

// MINUS
template <class T, unsigned Size>
inline Vector<T, Size> operator- (const BaseType& a, const Vector<T, Size>& b)
{
  return Vector<T, Size>(a)-=b;
}

template <class T, unsigned Size>
inline Vector<T, Size> operator- (const Vector<T, Size>& a, const BaseType& b)
{
  return Vector<T, Size>(a)-=b;
}

template <class T, unsigned Size>
inline Vector<T, Size> operator- (const Vector<T, Size>& a, const Vector<T, Size>& b)
{
  return Vector<T, Size>(a)-=b;
}

// MULTIPLY
template <class T, unsigned Size>
inline Vector<T, Size> operator* (const BaseType& a, const Vector<T, Size>& b)
{
  return Vector<T, Size>(b)*=a;
}

template <class T, unsigned Size>
inline Vector<T, Size> operator* (const Vector<T, Size>& a, const BaseType& b)
{
  return Vector<T, Size>(a)*=b;
}

template <class T, unsigned Size>
inline T operator* (const Vector<T, Size>& a, const Vector<T, Size>& b)
{
  T scalarProduct = T();
  for (unsigned iD = 0; iD < Size; ++iD) {
    scalarProduct += a[iD]*b[iD];
  }
  return scalarProduct;
}

/// < Operator returns true if all components are smaller
template <class T, unsigned Size>
inline bool operator< (const Vector<T, Size>& lhs, const Vector<T, Size>& rhs)
{
  bool smaller = true;
  for (unsigned i = 0; i < Size; ++i) {
    smaller &= (lhs[i] < rhs[i]);
  }
  return smaller;
}

/// > Operator returns true if all components are greater
template <class T, unsigned Size>
inline bool operator> (const Vector<T, Size>& lhs, const Vector<T, Size>& rhs)
{
  return rhs < lhs;
}

/// <= Operator returns true if all components are smaller or equal
template <class T, unsigned Size>
inline bool operator<=(const Vector<T, Size>& lhs, const Vector<T, Size>& rhs)
{
  bool smallerEq = true;
  for (unsigned i = 0; i < Size; ++i) {
    smallerEq &= (lhs[i] <= rhs[i]);
  }
  return smallerEq;
}

/// >= Operator returns true if all components are smaller or equal
template <class T, unsigned Size>
inline bool operator>=(const Vector<T, Size>& lhs, const Vector<T, Size>& rhs)
{
  return rhs <= lhs;
}


template <typename T>
using Scalar   = Vector<T,1>;

template <typename T>
using Vector2D = Vector<T,2>;

template <typename T>
using Vector3D = Vector<T,3>;

template <typename T>
using Vector4D = Vector<T,4>;

} // end namespace olb

#endif
