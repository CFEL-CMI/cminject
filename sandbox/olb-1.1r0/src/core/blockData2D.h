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
 * Dynamics for a generic 2D block data -- header file.
 */
#ifndef BLOCK_DATA_2D_H
#define BLOCK_DATA_2D_H


#include <stdio.h>
#include "blockStructure2D.h"
#include "serializer.h"


namespace olb {

// forward declaration is sufficient at this state, include is moved to .hh
template<typename BaseType> class BlockF2D;
template<typename T> class Cuboid2D;


/** A highly generic data structure class.
 *
 * Stored data is of type BaseType
 *
 * \param _size     data element size
 * \param _rawData  classic array representation of all elements
 * \param _filed    a matrix like representation of _rawData
 */
template<typename T, typename BaseType>
class BlockData2D : public BlockStructure2D, public Serializable {
protected:
  /// dimension of data element, vector, scalar, ...
  int _size;
  /// holds data as a 1D vector
  BaseType *_rawData;
  /** Pointer structure to a 3D data field.
   *  It can be interpreted as a 2D matrix[iX,iY] with elements of dimension _size.
   *  Those elements may be
   *    1. vector valued like velocity, (f0,f1,...,f8)
   *    2. scalar valued like density, pressure, ...
   *
   */
  BaseType ***_field;
public:
  virtual ~BlockData2D();
  /// Construct empty cuboid
  BlockData2D();
  /// Construct from cuboid
  BlockData2D(Cuboid2D<T>& cuboid, int size=1);
  /// Construct from X-Y node count
  BlockData2D(int nx, int ny, int size=1);
  /// Construct from Block Functor, attention!! operator() accesses functor data
  BlockData2D(BlockF2D<BaseType>& rhs);
  /// Copy Constructor
  BlockData2D(BlockData2D<T,BaseType> const& rhs);
  /// Assignment Operator
  BlockData2D<T,BaseType>& operator=(BlockData2D<T,BaseType> const& rhs);
  /// Move Operator
  BlockData2D<T,BaseType>& operator=(BlockData2D<T,BaseType>&& rhs);
  /// Move Constructor
  BlockData2D<T,BaseType>(BlockData2D<T,BaseType>&& rhs);
  /// Swap rhs Data into local fields
  void swap(BlockData2D<T,BaseType>& rhs);
  /// Memory Management
  bool isConstructed() const;
  void construct();
  void deConstruct();
  void reset();
  /// read and write access to data element [iX][iY][iSize]
  virtual BaseType& get(int iX, int iY, int iSize=0);
  /// read only access to data element [iX][iY][iSize]
  virtual BaseType const& get(int iX, int iY, int iSize=0) const;
  /// \return dataElement _rawData[ind], read and write
  BaseType& operator[] (int ind);
  /// \return dataElement _rawData[ind], read only
  BaseType const& operator[] (int ind) const;
  /// Write access to the memory of the data of the block data where (iX, iY) is the point providing the data iData
  bool* operator() (int iX, int iY, int iData);
  /// \return max of data, for vector valued data it determines the max component
  BaseType getMax();
  /// \return min of data, for vector valued data it determines the max component
  BaseType getMin();
  /// \return _rawData array
  BaseType* getRawData() const;
  /// \return _field
  BaseType*** getField() const;
  /// \return length of array _rawData or equivalent nX*nY*size
  virtual size_t getDataSize() const;
  /// Read only access to the dim of the data of the super structure
  int getSize() const;
  /// Number of data blocks for the serializable interface
  virtual std::size_t getNblock() const
  {
    return 4;
  };
  /// Binary size for the serializer
  virtual std::size_t getSerializableSize() const;
  /// Returns a pointer to the memory of the current block and its size for the serializable interface
  virtual bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode);
private:
  void allocateMemory(); // TODO
  void releaseMemory();
};


}  // namespace olb

#endif
