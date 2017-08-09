/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Mathias J. Krause
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
 * Dynamics for a generic 2D block -- header file.
 */
#ifndef BLOCK_STRUCTURE_2D_H
#define BLOCK_STRUCTURE_2D_H


namespace olb {

/** An empty hull with left bottom corner at (0,0,0).
 *
 * \param _nx extension in x direction
 * \param _ny extension in y direction
 * \param _nz extension in z direction
 *
 */
class BlockStructure2D {
protected:
  /// Block width
  int _nx;
  /// Block height
  int _ny;
public:
  BlockStructure2D(int nx, int ny) : _nx(nx), _ny(ny) {};
  /// Read only access to block width
  virtual int getNx() const
  {
    return _nx;
  };
  /// Read only access to block height
  virtual int getNy() const
  {
    return _ny;
  };
  virtual ~BlockStructure2D() { }
};

}  // namespace olb

#endif
