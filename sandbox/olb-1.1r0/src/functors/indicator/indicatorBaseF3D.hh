/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2016 Cyril Masquelier, Albert Mink, Mathias J. Krause, Benjamin Förster
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

#ifndef INDICATOR_BASE_F_3D_HH
#define INDICATOR_BASE_F_3D_HH


#include<cmath>
#include "indicatorBaseF3D.h"
#include "utilities/vectorHelpers.h"
#include "math.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define M_PI2 1.57079632679489661923

namespace olb {

template <typename S>
IndicatorF3D<S>::IndicatorF3D() : GenericF<bool,S>(1, 3)
{}

template <typename S>
Vector<S,3>& IndicatorF3D<S>::getMin()
{
  return _myMin;
}

template <typename S>
Vector<S,3>& IndicatorF3D<S>::getMax()
{
  return _myMax;
}

template <typename S>
bool IndicatorF3D<S>::distance(S& distance, const Vector<S,3>& origin, const Vector<S,3>& direction, int iC)
{
  bool originValue;
  bool currentValue;
  S precision = .0001;
  S pitch = 0.5;

  // start at origin and move into given direction
  Vector<S,3> currentPoint(origin);

  (*this)(&originValue, origin.data);
  (*this)(&currentValue, currentPoint.data);

  while (currentValue == originValue && isInsideBox(currentPoint)) {
    currentPoint += direction;
    // update currentValue until the first point on the other side (inside/outside) is found
    (*this)(&currentValue, currentPoint.data);
  }

  // return false if no point was found in given direction
  if (!isInsideBox(currentPoint) && !originValue) {
    return false;
  }


  while (pitch >= precision) {
    if (!isInsideBox(currentPoint) && originValue) {
      currentPoint -= pitch * direction;
      pitch /= 2.;
    } else {
      (*this)(&currentValue, currentPoint.data);
      if (currentValue == originValue) {
        currentPoint += pitch * direction;
        pitch /= 2.;
      } else {
        currentPoint -= pitch * direction;
        pitch /= 2.;
      }
    }
  }

  distance = (currentPoint - origin).norm();
  return true;
}

/// dot product, only valid in 3d
template <typename T>
S dotProduct3D(const Vector<S,3>& a, const Vector<S,3>& b) // Was commented out
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

template <typename S>
bool IndicatorF3D<S>::rotOnAxis(Vector<S,3>& vec_rot, const Vector<S,3>& vec, const Vector<S,3>& axis, S& theta)
{
  /// http://mathworld.wolfram.com/RodriguesRotationFormula.html https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula

  //Vector<S,3> axisN(axis*(1/(axis).norm())); // normalize rotation axis
  Vector<S,3> axisN(axis*(1/const_cast<Vector<S,3>&> (axis).norm())); // normalize rotation axis

  vec_rot = vec ;

  Vector<S,3> crossProd;

  crossProd[0] = axisN[1]*vec[2] - axisN[2]*vec[1];
  crossProd[1] = axisN[2]*vec[0] - axisN[0]*vec[2];
  crossProd[2] = axisN[0]*vec[1] - axisN[1]*vec[0];

  S dotProd = axisN[0]*vec[0] + axisN[1]*vec[1] + axisN[2]*vec[2];

  //v_rot = std::cos(theta)*vec + (crossProd)*std::sin(theta) + axisN*(dotProduct3D(axisN,vec))*(1 - std::cos(theta));
  vec_rot = std::cos(theta)*vec + (crossProd)*std::sin(theta) + axisN*(dotProd)*(1 - std::cos(theta));

  return true;

}

template <typename S>
bool IndicatorF3D<S>::normal(Vector<S,3>& normal, const Vector<S,3>& origin, const Vector<S,3>& direction, int iC)
{
  //OstreamManager clout(std::cout,"IndicatorF3D");
  std::cout << "Calculating IndicatorF3D Normal" << std::endl;

  bool originValue;
  (*this)(&originValue, origin.data);
  Vector<S,3> currentPoint(origin);

  S precision = .0001;

  S dist;
  distance(dist, origin, direction, iC);


  //Vector<S,3> POS(origin + dist*direction*(1/const_cast<Vector<S,3>&> (direction).norm())); //Point on Surface
  //direction = direction*(1/const_cast<Vector<S,3>&> (direction).norm());
  Vector<S,3> directionN(direction*(1/const_cast<Vector<S,3>&> (direction).norm()));
  Vector<S,3> POS(origin + dist*directionN); //Point on Surface

  /// find perpendicular vector to direction
  Vector<S,3> directionPerp;
  if(    (util::nearZero(directionN[0]) && util::nearZero(directionN[1]) && !util::nearZero(directionN[2]))
      || (util::nearZero(directionN[0]) && !util::nearZero(directionN[1]) && util::nearZero(directionN[2]))    ) {
    directionPerp[0] = 1;
    directionPerp[1] = 0;
    directionPerp[2] = 0;
  } else if ( !util::nearZero(directionN[0]) && util::nearZero(directionN[1]) && util::nearZero(directionN[2]) ) {
    directionPerp[0] = 0;
    directionPerp[1] = 0;
    directionPerp[2] = 1;
  } else if ( ( !util::nearZero(directionN[0]) || !util::nearZero(directionN[1]) ) && !util::nearZero(directionN[2]) ) {
    directionPerp[0] = directionN[0];
    directionPerp[1] = directionN[1];
    directionPerp[2] = -(directionN[0] + directionN[1])/directionN[2];
  } else {
    std::cout << "Error: unknown case for perpendicular check" << std::endl;
    return false;
  }

  Vector<S,3> directionPerpN(directionPerp*(1/(directionPerp).norm()));


  Vector<S,3> point1;
  Vector<S,3> point2;
  Vector<S,3> point3;

  bool currentValue;

  /// Loop 3 times to find three points on the surface to use for normal calc.
  /// orthogonal to direction vector 120� to each other


  for (int n: {
         0,120,240
       }) {
    S thetaMain = n*M_PI/180.;

    /// rotate directionPerpN through 3 angles {0,120,240}
    Vector<S,3> perp;
    rotOnAxis(perp, directionPerpN, directionN, thetaMain);
    Vector<S,3> perpPoint(POS + perp);

    //std::cout << "perp = [" << perp[0] << "," << perp[1]  << "," << perp[2] << "]" << std::endl;

    S rotate = 90.;
    S pitch = rotate/2.;

    Vector<S,3> vec(perp);

    Vector<S,3> rotAxis;

    //rotAxis = cross(perp,directionN);
    rotAxis[0] = perp[1]*directionN[2] - perp[2]*directionN[1];
    rotAxis[1] = perp[2]*directionN[0] - perp[0]*directionN[2];
    rotAxis[2] = perp[0]*directionN[1] - perp[1]*directionN[0];

    //normal = rotAxis;
    //return true;
    //std::cout << "rotAxis = [" << rotAxis[0] << "," << rotAxis[1]  << "," << rotAxis[2] << "]" << std::endl;

    /// Find 'positive' angle
    Vector<S,3> testPOS;
    S testAngle(45.*M_PI/180.);
    rotOnAxis(testPOS, perp, rotAxis, testAngle);
    Vector<S,3> testPoint( POS + testPOS);
    //std::cout << "testPOS = [" << testPOS[0] << "," << testPOS[1]  << "," << testPOS[2] << "]" << std::endl;
    //std::cout << "testPoint = [" << testPoint[0] << "," << testPoint[1]  << "," << testPoint[2] << "]" << std::endl;

    S distTestPoint = (testPoint - origin).norm();
    S distPerpPoint = (perpPoint - origin).norm();

    S mod = 0;
    if (distTestPoint < distPerpPoint) { // pos. angle rotates towards
      mod = -1;
    } else {
      mod = 1;
    }

    while (std::abs(pitch) >= precision) {

      S theta(pitch*M_PI/180);

      currentPoint = POS + vec;
      (*this)(&currentValue, currentPoint.data);

      S temp;
      if (currentValue == originValue) {
        temp = mod*theta;
        rotOnAxis(vec, vec, rotAxis, temp);
      }  else {
        temp = -mod*theta;
        rotOnAxis(vec, vec, rotAxis, temp);
      }
      pitch /= 2.;



    }

    if (n == 0) {
      point1 = currentPoint;
    } else if (n == 120) {
      point2 = currentPoint;
    } else if (n == 240) {
      point3 = currentPoint;
    } else {
      std::cout << "Something broke" << std::endl;
      return false;
    }

  }

  /// Calculate Normal
  Vector<S,3> vec1 (point1 - point2);
  Vector<S,3> vec2 (point1 - point3);

  normal[0] = -(vec1[1]*vec2[2] - vec1[2]*vec2[1]);
  normal[1] = -(vec1[2]*vec2[0] - vec1[0]*vec2[2]);
  normal[2] = -(vec1[0]*vec2[1] - vec1[1]*vec2[0]);


  //S dist;
  //Vector<S,3> dist;
  //distance(dist, origin, direction, iC);
  //normal = Vector<S,3>(dist);
  //normal = POS;
  //normal = directionPerpN;
  //normal = directionN;
  return true;

}

template <typename S>
bool IndicatorF3D<S>::isInsideBox(Vector<S,3> point)
{
  return point >= _myMin && point <= _myMax;
}


// identity to "store results"
template <typename S>
IndicatorIdentity3D<S>::IndicatorIdentity3D(IndicatorF3D<S>& f) : _f(f)
{
  this->_myMin = _f.getMin();
  this->_myMax = _f.getMax();
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename S>
bool IndicatorIdentity3D<S>::operator() (bool output[], const S input[])
{
  return _f(output, input);
}

template <typename T, typename S>
SmoothIndicatorF3D<T,S>::SmoothIndicatorF3D()
  : AnalyticalF3D<T,S>(1),
    _myMin(S()), _myMax(S()), _center(S()),
    _vel(S()), _acc(S()), _theta(S()), _omega(S()), _alpha(S()), _mass(S()), _mofi(S())
{ }

template <typename T, typename S>
Vector<S,3>& SmoothIndicatorF3D<T,S>::getMin()
{
  return _myMin;
}

template <typename T, typename S>
Vector<S,3>& SmoothIndicatorF3D<T,S>::getMax()
{
  return _myMax;
}

template <typename T, typename S>
Vector<S,3>& SmoothIndicatorF3D<T,S>::getCenter()
{
  return _center;
};

template <typename T, typename S>
Vector<S,3>& SmoothIndicatorF3D<T,S>::getVel()
{
  return _vel;
}

template <typename T, typename S>
Vector<S,3>& SmoothIndicatorF3D<T,S>::getAcc()
{
  return _acc;
}

template <typename T, typename S>
Vector<S,6>& SmoothIndicatorF3D<T,S>::getAcc2()
{
  return _acc2;
}

template <typename T, typename S>
Vector<S,3>& SmoothIndicatorF3D<T,S>::getTheta()
{
  return _theta;
}

template <typename T, typename S>
Vector<S,3>& SmoothIndicatorF3D<T,S>::getOmega()
{
  return _omega;
}

template <typename T, typename S>
Vector<S,3>& SmoothIndicatorF3D<T,S>::getAlpha()
{
  return _alpha;
}

template <typename T, typename S>
S& SmoothIndicatorF3D<T,S>::getMass()
{
  return _mass;
}

template <typename T, typename S>
S& SmoothIndicatorF3D<T,S>::getMofi()
{
  return _mofi;
}

template <typename T, typename S>
S SmoothIndicatorF3D<T,S>::getDiam()
{
  return 2*_radius;
};

template <typename T, typename S>
S SmoothIndicatorF3D<T,S>::getRadius()
{
  return _radius+_epsilon;
};

template <typename T, typename S>
void SmoothIndicatorF3D<T,S>::setCenter(S centerX, S centerY, S centerZ)
{
  _center[0] = centerX;
  _center[1] = centerY;
  _center[2] = centerZ;
};

template <typename T, typename S>
void SmoothIndicatorF3D<T,S>::setTheta(S thetaX, S thetaY, S thetaZ)
{
  _theta[0] = thetaX;
  _theta[1] = thetaY;
  _theta[2] = thetaZ;
};

// identity to "store results"
template <typename T, typename S>
SmoothIndicatorIdentity3D<T,S>::SmoothIndicatorIdentity3D(SmoothIndicatorF3D<T,S>& f)
  : _f(f)
{
  for ( int i=0; i<3; i++) {
    this->_myMin[i] = _f.getMin()[i];
    this->_myMax[i] = _f.getMax()[i];
  }
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename T, typename S>
bool SmoothIndicatorIdentity3D<T,S>::operator() (T output[], const S input[])
{
  _f(output, input);
  return true;
}





///////////////////////////////////////
/////     ParticleIndicator3D     /////
///////////////////////////////////////

template <typename T, typename S>
ParticleIndicatorF3D<T,S>::ParticleIndicatorF3D()
  : AnalyticalF3D<T,S>(1),
    _myMin(S()), _myMax(S()), _pos(S()),
    _vel(S()), _acc(S()), _acc2(S()), _theta(S()), _omega(S()), _alpha(S()), _alpha2(S()), _mofi(S()), _rotMat(S()), _mass(S()), _radius(S())
{ }

template <typename T, typename S>
Vector<S,3>& ParticleIndicatorF3D<T,S>::getMin()
{
  return _myMin;
};


template <typename T, typename S>
Vector<S,3>& ParticleIndicatorF3D<T,S>::getMax()
{
  return _myMax;
};

template <typename T, typename S>
Vector<S,3>& ParticleIndicatorF3D<T,S>::getPos()
{
  return _pos;
};

template <typename T, typename S>
Vector<S,3>& ParticleIndicatorF3D<T,S>::getVel()
{
  return _vel;
};

template <typename T, typename S>
Vector<S,3>& ParticleIndicatorF3D<T,S>::getAcc()
{
  return _acc;
};

template <typename T, typename S>
Vector<S,3>& ParticleIndicatorF3D<T,S>::getAcc2()
{
  return _acc2;
};

template <typename T, typename S>
Vector<S,3>& ParticleIndicatorF3D<T,S>::getTheta()
{
  return _theta;
};

template <typename T, typename S>
Vector<S,3>& ParticleIndicatorF3D<T,S>::getOmega()
{
  return _omega;
};

template <typename T, typename S>
Vector<S,3>& ParticleIndicatorF3D<T,S>::getAlpha()
{
  return _alpha;
};

template <typename T, typename S>
Vector<S,3>& ParticleIndicatorF3D<T,S>::getAlpha2()
{
  return _alpha2;
};

template <typename T, typename S>
Vector<S,3>& ParticleIndicatorF3D<T,S>::getMofi()
{
  return _mofi;
};

template <typename T, typename S>
Vector<S,9>& ParticleIndicatorF3D<T,S>::getRotationMat()
{
  return _rotMat;
};

template <typename T, typename S>
S& ParticleIndicatorF3D<T,S>::getMass()
{
  return _mass;
};


template <typename T, typename S>
S& ParticleIndicatorF3D<T,S>::getRadius()
{
  return _radius;
};


template <typename T, typename S>
S ParticleIndicatorF3D<T,S>::getDiam()
{
  return 2.*_radius;
};


template <typename T, typename S>
ParticleIndicatorIdentity3D<T,S>::ParticleIndicatorIdentity3D(ParticleIndicatorF3D<T,S>& f)
  : _f(f)
{
  for ( int i=0; i<3; i++) {
    this->_myMin[i] = _f.getMin()[i];
    this->_myMax[i] = _f.getMax()[i];
  }
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename T, typename S>
bool ParticleIndicatorIdentity3D<T,S>::operator() (T output[], const S input[])
{
  _f(output, input);
  return true;
}

} // namespace olb

#endif
