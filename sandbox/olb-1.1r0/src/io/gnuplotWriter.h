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

#ifndef GNUPLOT_WRITER_H
#define GNUPLOT_WRITER_H

#include <iomanip>
#include <iostream>
#include <list>

namespace olb {

template< typename T >
class Gnuplot {
public:
  /// Constructor with name for output files
  /// boolean true for real-time plotting //WARNING: experimental!
  Gnuplot(std::string name, bool liveplot = false);

  /// initialises the data file
  void init();

  /// sets the data and plot file for two doubles (x and y)
  void setData(T xValue, T yValue, std::string name = "", std::string key = "");
  /// if no x value is given, it is just an increasing integer
  void setData(bool noXvalue, T yValue, std::string name = "", std::string key = "");

  /// sets the data and plot file for a double and a list of doubles (x and {y1,y2,...})
  void setData(T xValue, std::list<T> yValues, std::list<std::string> names = {""}, std::string key = "right");
  /// if no x value is given, it is just an increasing integer
  void setData(bool noXvalue, std::list<T> yValues, std::list<std::string> names = {""}, std::string key = "right");

  /// writes an PDF
  void writePDF();

  /// writes PNGs
  /// usage: first argument: numbering of png file (optional), second argument: range for the x axis (optional)
  /// no arguments: writes in one file with adaptive xrange
  void writePNG(int iT = -1, double xRange = -1);

private:
  std::string _name;
  bool _liveplot;
  std::string _type;
  std::string _dataFile;
  std::string _dir;
  std::list<std::string> _nameList;
  std::string _key;
  bool _init = true;
  unsigned int _dataSize = 0;
  int _iT = -1;
  double _xRange = -1;
  T _time = 0.;

  int _rank = 0;

  /// writes a plot file for type {"plot", "png", "pdf")
  void writePlotFile(std::string type);

  /// writes the data file for two doubles (x and y)
  void writeDataFile(T xValue, T yValue);

  /// writes the data file for one double and a list of doubles (x and y1,y2,...)
  void writeDataFile(T xValue, std::list<T> list);

  /// system command to start gnuplot (LINUX ONLY!)
  void startGnuplot(std::string plotFile);

};

}  // namespace olb

#endif

