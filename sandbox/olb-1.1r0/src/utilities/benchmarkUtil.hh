#ifndef BENCHMARK_UTIL_HH
#define BENCHMARK_UTIL_HH

#include "benchmarkUtil.h"
#include <iostream>
#include <string>
#include <sstream>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <cassert>


namespace olb {

namespace util {


/////////// Class ValueTracer ////////////////////////

template<typename T>
ValueTracer<T>::ValueTracer(T u, T L, T _epsilon)
  : deltaT((int)(L/u/2.)),
    epsilon(_epsilon),
    t(0),
    converged(false),
    clout(std::cout,"ValueTracer")
{ }

template<typename T>
ValueTracer<T>::ValueTracer(int _deltaT, T _epsilon)
  : deltaT(_deltaT),
    epsilon(_epsilon),
    t(0),
    converged(false),
    clout(std::cout,"ValueTracer")
{ }

template<typename T>
int ValueTracer<T>::getDeltaT() const
{
  return deltaT;
}

template<typename T>
void ValueTracer<T>::takeValue(T val, bool doPrint)
{
  values.push_back(val);
  if ((int)values.size() > abs(deltaT)) {
    values.erase(values.begin());
    if (doPrint && t%deltaT==0) {
      T average = computeAverage();
      T stdDev = computeStdDev(average);
      clout << "average=" << average << "; stdDev/average=" << stdDev/average << std::endl;
    }
  }
  ++t;
}

template<typename T>
void ValueTracer<T>::resetScale(T u, T L)
{
  t = t%deltaT;
  deltaT = (int) (L/u/2.);
  if ( (int)values.size() > abs(deltaT) ) {
    values.erase(values.begin(), values.begin() + (values.size()-deltaT) );
  }
}

template<typename T>
void ValueTracer<T>::resetValues()
{
  t = 0;
  if ((int)values.size() > 0) {
    values.erase(values.begin(), values.begin() + values.size() );
  }
}

template<typename T>
bool ValueTracer<T>::hasConverged() const
{
  if ((int)values.size() < abs(deltaT)) {
    return false;
  } else {
    T average = computeAverage();
    T stdDev = computeStdDev(average);
    if (!std::isnan(stdDev/average)) {
      return stdDev/average < epsilon;
    } else {
      clout << "simulation diverged." << std::endl;
      return true;
    }
  }
}

template<typename T>
bool ValueTracer<T>::hasConvergedMinMax() const
{
  if ((int)values.size() < abs(deltaT)) {
    return false;
  } else {
    T minEl = *min_element(values.begin(), values.end());
    T maxEl = *max_element(values.begin(), values.end());
    T average = computeAverage();
    return (maxEl - minEl)/average < epsilon;
  }
}

template<typename T>
T ValueTracer<T>::computeAverage() const
{
  return accumulate(values.begin(), values.end(), 0.) / values.size();
}

template<typename T>
T ValueTracer<T>::computeStdDev(T average) const
{
  int n = values.size();
  T sqrDev = 0.;
  for (int i=0; i<n; ++i) {
    sqrDev += (values[i]-average)*(values[i]-average);
  }
  return sqrt(sqrDev/(n-1));
}

template<typename T>
void ValueTracer<T>::setEpsilon(T epsilon_)
{
  epsilon = epsilon_;
}



/////////// Class BisectStepper ////////////////////////

template<typename T>
BisectStepper<T>::BisectStepper(T _iniVal, T _step)
  : iniVal(_iniVal), step(_step), state(first), clout(std::cout,"BisectStepper")
{
  if (util::nearZero(step)) {
    step = iniVal/5.;
  }
  assert(step>0.);
}

template<typename T>
T BisectStepper<T>::getVal(bool stable, bool doPrint)
{
  std::stringstream message;
  switch (state) {
  case first:
    if (stable) {
      currentVal = iniVal+step;
      state = up;
      message << "[" << iniVal << ",infty]";
    } else {
      currentVal = iniVal-step;
      state = down;
      message << "[-infty," << iniVal << "]";
    }
    break;
  case up:
    if (stable) {
      message << "[" << currentVal << ",infty]";
      currentVal += step;
    } else {
      lowerVal = currentVal-step;
      upperVal = currentVal;
      currentVal = 0.5*(lowerVal+upperVal);
      state = bisect;
      message << "[" << lowerVal << "," << upperVal << "]";
    }
    break;
  case down:
    if (!stable) {
      message << "[-infty," << currentVal << "]";
      currentVal -= step;
    } else {
      lowerVal = currentVal;
      upperVal = currentVal+step;
      currentVal = 0.5*(lowerVal+upperVal);
      state = bisect;
      message << "[" << lowerVal << "," << upperVal << "]";
    }
    break;
  case bisect:
    if (stable) {
      lowerVal = currentVal;
    } else {
      upperVal = currentVal;
    }
    currentVal = 0.5*(lowerVal+upperVal);
    message << "[" << lowerVal << "," << upperVal << "]";
    break;
  }
  if (doPrint) {
    clout << "Value in range " << message.str() << std::endl;
  }
  return currentVal;
}

template<typename T>
bool BisectStepper<T>::hasConverged(T epsilon) const
{
  return (state==bisect) && ((upperVal-lowerVal)/lowerVal < epsilon);
}

/////////// Class Newton 1D ////////////////////////

template<typename T>
Newton1D<T>::Newton1D(AnalyticalF1D<T,T>& f, T yValue, T eps, int maxIterations) : _f(f), _df(f,eps), _yValue(yValue), _eps(eps), _maxIterations(maxIterations)
{
}

template<typename T>
T Newton1D<T>::solve(T startValue, bool print)
{
  T fValue[1], dfValue[1], x[1];
  x[0] = startValue;
  _f(fValue,x);
  T eps = fabs(fValue[0] - _yValue);
  for (int i=0; i<_maxIterations && eps>_eps; i++) {
    _f(fValue,x);
    _df(dfValue,x);
    eps = fabs(fValue[0] - _yValue);
    if (print) {
      std::cout << "newtonStep=" << i << "; eps=" << eps << "; x=" << x[0]  << "; f(x)="<< fValue[0] << "; df(x)="<< dfValue[0] << std::endl;
    }
    x[0] -= (fValue[0] - _yValue)/dfValue[0];
  }
  return x[0];
}

template<typename T>
TrapezRuleInt1D<T>::TrapezRuleInt1D(AnalyticalF1D<T,T>& f) : _f(f)
{
}

template<typename T>
T TrapezRuleInt1D<T>::integrate(T min, T max, int nSteps)
{
  T integral[1], tmp[1], min_[1], max_[1], x[1], deltax[0];
  min_[0] = min;
  max_[0] = max;
  deltax[0] = (max - min) / nSteps;
  _f(tmp, min);
  integral[0] = 0.5 * tmp[0];
  x[0] = min;
  for (int i = 1; i < nSteps; i++) {
    x[0] += deltax[0];
    _f(tmp, x);
    integral[0] += tmp[0];
  }

  _f(tmp, max);
  integral[0] += 0.5*tmp[0];
  integral[0] *= deltax[0];

  std::cout << "Itnegral=" <<  integral[0] << std::endl;

  return integral[0];
}

} // namespace util

} // namespace olb

#endif
