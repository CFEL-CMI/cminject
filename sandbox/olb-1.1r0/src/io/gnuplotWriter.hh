/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Fabian Klemens
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

#ifndef GNUPLOT_WRITER_HH
#define GNUPLOT_WRITER_HH

#include <fstream>
#include <iostream>
#include <unistd.h>
#include "core/singleton.h"
#include "io/fileName.h"
#include "gnuplotWriter.h"
#include "utilities/vectorHelpers.h"

namespace olb {
/// Constructor with name of outputFiles
/// boolean true for real-time plotting //WARNING: experimental!
template< typename T >
Gnuplot<T>::Gnuplot(std::string name, bool liveplot)
  : _name(name),
    _dataFile(singleton::directories().getGnuplotOutDir()+"data/"+_name+".dat"),
    _dir(singleton::directories().getGnuplotOutDir())
{
  _liveplot = liveplot;
  if (singleton::mpi().getRank() == _rank) {
    std::ofstream fout;

    ///add (new) data file
    fout.open(_dataFile.c_str(), std::ios::trunc);
    fout.close();
  }
}

/// writes the data and plot file for two doubles (x and y)
template< typename T >
void Gnuplot<T>::setData(T xValue, T yValue, std::string name, std::string key)
{
  if (_init) {
    _dataSize = 1;
    _key = key;
    _nameList = {name};

    if (_liveplot) {
      writePlotFile("plot");
    }
  }
  writeDataFile(xValue, yValue);

  if (_liveplot && _init) {
    startGnuplot("plot");
  }

  _init = false;
  return;
}

/// writes the data and plot file for two doubles (x and y), where x is increasing integer
template< typename T >
void Gnuplot<T>::setData(bool noXvalue, T yValue, std::string name, std::string key)
{
  T xValue = _time;
  setData(xValue, yValue, name, key);
  _time++;
}

/// writes the data and plot file for a double and a list of doubles (x and y1,y2,...)
template< typename T >
void Gnuplot<T>::setData(T xValue, std::list<T> yValues, std::list<std::string> names, std::string key)
{
  if (_init) {
    _dataSize = yValues.size();
    _key = key;
    _nameList = names;
    if (_nameList.size() != _dataSize) {
      for (unsigned int i = 0; i < _dataSize - 1; i++) {
        _nameList.push_back("");
      }
    }
    if (_liveplot) {
      writePlotFile("plot");
    }
  }
  writeDataFile(xValue,yValues);

  if (_liveplot && _init) {
    startGnuplot("plot");
  }

  _init = false;
  return;
}

/// writes the data and plot file for a double and a list of doubles (x and y1,y2,...), where x is increasing integer
template< typename T >
void Gnuplot<T>::setData(bool noXvalue, std::list<T> yValues, std::list<std::string> names, std::string key)
{
  T xValue = _time;
  setData(xValue, yValues, names, key);
  _time++;
}


/// writes an PDF
template< typename T >
void Gnuplot<T>::writePDF()
{
  if (!_init) {
    writePlotFile("pdf");
    startGnuplot("plotPDF");
  }
  return;
}


/// writes PNGs
/// usage: first argument: numbering of png file, second argument: range for the x axis
/// no arguments: writes consecutive numbers with adaptive xrange
template< typename T >
void Gnuplot<T>::writePNG(int iT, double xRange)
{
  if (!_init) {
    _iT = iT;
    _xRange = xRange;

    writePlotFile("png");
    startGnuplot("plotPNG");
  }
  return;
}

template< typename T >
void Gnuplot<T>::writePlotFile(std::string type)
{
  if (singleton::mpi().getRank() == _rank ) {
    std::ofstream fout;

    std::string plotFile;
    if (_liveplot && type == "plot") {
      plotFile = singleton::directories().getGnuplotOutDir()+"data/plot.p";
    } else if (type == "pdf") {
      plotFile = singleton::directories().getGnuplotOutDir()+"data/plotPDF.p";
    } else if (type == "png") {
      plotFile = singleton::directories().getGnuplotOutDir()+"data/plotPNG.p";
    } else {
      std::cout << "WARNING: invalid Gnuplot type={'', 'plot'; 'pdf', 'png'}" << std::endl;
      exit(-1);
    }

    fout.open(plotFile.c_str(), std::ios::trunc);
    fout << "set key " << _key << "\n";

    if (type=="pdf") {
      fout << "set terminal pdf enhanced" << "\n"
           << "set output '"<<_dir<<_name<<".pdf'" << "\n";
    }
    if (type=="png") {
      if ( !util::nearZero(_xRange+1) ) {
        fout << "set xr[0:"<< _xRange <<"]" << "\n";
      }
      fout << "set terminal png" << "\n"
           << "set output '"<<_dir<<_name;
      if (_iT != -1) {
        fout <<"_"<<_iT;
      }
      fout <<".png'" << "\n";
    }
    std::list<std::string>::iterator string = _nameList.begin();
    fout << "plot '"<<_dataFile<<"' u 1:2 w l t '"<< *string << "'";
    for (unsigned int i = 0; i < _dataSize-1; ++i) {
      ++string;
      fout << ", '"<<_dataFile<<"' u 1:" << i+3 << " w l t '" << *string << "'";
    }
    fout << "\n";
    if (_liveplot && type=="plot") {
      fout << "pause -1" << "\n"
           << "reread" << "\n";
    }
    fout.close();
  }
  return;
}


/// writes the data file for two doubles (x and y)
template< typename T >
void Gnuplot<T>::writeDataFile(T xValue, T yValue)
{
  if (singleton::mpi().getRank() == _rank) {
    std::ofstream fout;
    fout.precision(6);
    fout.open(_dataFile.c_str(), std::ios::app);
    fout << xValue
         << " "
         << yValue
         << std::endl;
    fout.close();
  }
  return;
}


/// writes the data file for one double and a list of doubles (x and y1,y2,...)
template< typename T >
void Gnuplot<T>::writeDataFile(T xValue, std::list<T> list)
{
  if (singleton::mpi().getRank() == _rank) {
    std::ofstream fout;
    fout.precision(6);
    fout.open(_dataFile.c_str(), std::ios::app);
    fout << xValue;
    typename std::list<T>::iterator yValue;
    for ( yValue = list.begin(); yValue!=list.end(); ++yValue) {
      fout << " " << *yValue;
    }
    fout << "\n";
    fout.close();
  }
  return;
}


/// system command to start gnuplot (LINUX ONLY!)
template< typename T >
void Gnuplot<T>::startGnuplot(std::string plotFile)
{
#ifdef WIN32
  std::cout << "GNUPLOT WORKS ONLT WITH LINUX" << std::endl;
//  exit (EXIT_FAILURE);
  return;
#endif
#ifndef WIN32
  if (singleton::mpi().getRank() == _rank) {
    if (!system(NULL)) {
      exit (EXIT_FAILURE);
    }

    const std::string command = "gnuplot -persistent "+_dir+"data/"+plotFile+".p > /dev/null &";
    if ( system(command.c_str()) ) {
      std::cout << "Error at GnuplotWriter" << std::endl;
    }
  }
  return;
#endif
}
}  // namespace olb

#endif

