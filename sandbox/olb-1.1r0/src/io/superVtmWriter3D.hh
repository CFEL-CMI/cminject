/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016-2017 Albert Mink, Maximilian Gaedtke, Markus Morhard Mathias J. Krause
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
 * A method to write vtk data for cuboid geometries
 * (only for uniform grids) -- generic implementation.
 */

#ifndef SUPER_VTM_WRITER_3D_HH
#define SUPER_VTM_WRITER_3D_HH

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include "core/singleton.h"
#include "communication/loadBalancer.h"
#include "geometry/cuboidGeometry3D.h"
#include "communication/mpiManager.h"
#include "io/fileName.h"
#include "io/superVtmWriter3D.h"
#include "io/base64.h"

#include <stdio.h>
#include <assert.h>
#include "zlib.h"



namespace olb {


template<typename T, typename W>
SuperVTMwriter3D<T,W>::SuperVTMwriter3D( std::string name)
  : clout( std::cout,"SuperVTMwriter3D" ), _createFile(false), _name(name)
{}

template<typename T, typename W>
void SuperVTMwriter3D<T,W>::write(int iT)
{
  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif

  //  !!!!!!!!!!! check whether _pointerVec is empty
  if ( _pointerVec.empty() ) {
    clout << "Error: Did you add a Functor ?";
  } else {
    // to get first element _pointerVec
    // problem if functors with different SuperStructure are stored
    // since till now, there is only one origin
    auto it_begin = _pointerVec.cbegin();
    CuboidGeometry3D<T> const& cGeometry = (**it_begin).getSuperStructure().getCuboidGeometry();
    // no gaps between vti files (cuboids)
    (**it_begin).getSuperStructure().communicate();
    LoadBalancer<T>& load = (**it_begin).getSuperStructure().getLoadBalancer();

    // PVD, owns all
    if ( rank == 0 ) {
      std::string pathPVD = singleton::directories().getVtkOutDir()
                            + createFileName( _name ) + ".pvd";
      dataPVDmaster( iT, pathPVD,  "data/" + createFileName( _name, iT ) + ".vtm" );

      std::string pathVTM = singleton::directories().getVtkOutDir()
                            + "data/" + createFileName( _name, iT ) + ".vtm";
      preambleVTM(pathVTM);
      for (int iC = 0; iC < cGeometry.getNc(); iC++) {
        dataVTM( iC, pathVTM, createFileName( _name, iT, iC) + ".vti" );
      }
      closeVTM(pathVTM);
    }
    // VTI, each process writes his cuboids
    int originLatticeR[4] = {int()};
    for (int iCloc = 0; iCloc < load.size(); iCloc++) {
      int nx = cGeometry.get(load.glob(iCloc)).getNx();
      int ny = cGeometry.get(load.glob(iCloc)).getNy();
      int nz = cGeometry.get(load.glob(iCloc)).getNz();
      // to be changed into the following line once local refinement has been implemented
      // double deltaX = cGeometry.get(load.glob(iCloc)).getDeltaR();
      T delta = cGeometry.getMotherCuboid().getDeltaR();

      std::string fullNameVTI = singleton::directories().getVtkOutDir() + "data/"
                                + createFileName( _name, iT, load.glob(iCloc) ) + ".vti";

      // get dimension/extent for each cuboid
      originLatticeR[0] = load.glob(iCloc);
      T originPhysR[3] = {T()};
      cGeometry.getPhysR(originPhysR,originLatticeR);

      preambleVTI(fullNameVTI, -1,-1,-1,nx,ny,nz,
                  originPhysR[0],originPhysR[1],originPhysR[2], delta);
      for (auto it : _pointerVec) {
        dataArray(fullNameVTI, (*it), load.glob(iCloc), nx,ny,nz);
      }
      closePiece(fullNameVTI);
      closeVTI(fullNameVTI);
    }
  }
}

template<typename T, typename W>
void SuperVTMwriter3D<T,W>::write(SuperF3D<T,W>& f, int iT)
{
  CuboidGeometry3D<T> const& cGeometry = f.getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>& load = f.getSuperStructure().getLoadBalancer();
  // no gaps between vti files (cuboids)
  f.getSuperStructure().communicate();
  T delta = cGeometry.getMotherCuboid().getDeltaR();

  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif

  // write a pvd file, which links all vti files
  // each vti file is written by one thread, which may own severals cuboids
  if ( rank == 0 ) {
    // master only
    std::string pathVTM = singleton::directories().getVtkOutDir()
                          + createFileName( f.getName(), iT )  + ".vtm";

    preambleVTM(pathVTM);
    for (int iC = 0; iC < cGeometry.getNc(); iC++) {
      std::string nameVTI = "data/" + createFileName( f.getName(), iT, iC) + ".vti";
      // puts name of .vti piece to a .pvd file [fullNamePVD]
      dataVTM( iC, pathVTM, nameVTI );
    }
    closeVTM(pathVTM);
  } // master only

  int originLatticeR[4] = {int()};
  for (int iCloc = 0; iCloc < load.size(); iCloc++) {
    // cuboid
    int nx = cGeometry.get(load.glob(iCloc)).getNx();
    int ny = cGeometry.get(load.glob(iCloc)).getNy();
    int nz = cGeometry.get(load.glob(iCloc)).getNz();
    // to be changed into the following line once local refinement has been implemented
    // double deltaX = cGeometry.get(load.glob(iCloc)).getDeltaR();

    std::string fullNameVTI = singleton::directories().getVtkOutDir() + "data/"
                              + createFileName( f.getName(), iT, load.glob(iCloc) ) + ".vti";

    // get dimension/extent for each cuboid
    originLatticeR[0] = load.glob(iCloc);
    T originPhysR[3] = {T()};
    cGeometry.getPhysR(originPhysR,originLatticeR);

    preambleVTI(fullNameVTI, -1,-1,-1, nx,ny,nz,
                originPhysR[0],originPhysR[1],originPhysR[2], delta);

    dataArray(fullNameVTI, f, load.glob(iCloc), nx,ny,nz);
    closePiece(fullNameVTI);
    closeVTI(fullNameVTI);
  } // cuboid
}

template<typename T, typename W>
void SuperVTMwriter3D<T,W>::createMasterFile()
{
  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif
  if ( rank == 0 ) {
    std::string fullNamePVDmaster = singleton::directories().getVtkOutDir()
                                    + createFileName( _name ) + ".pvd";
    preamblePVD(fullNamePVDmaster);
    closePVD(fullNamePVDmaster);
    _createFile = true;
  }
}

template<typename T, typename W>
void SuperVTMwriter3D<T,W>::addFunctor(SuperF3D<T,W>& f)
{
  _pointerVec.push_back(&f);
}

template<typename T, typename W>
void SuperVTMwriter3D<T,W>::clearAddedFunctors()
{
  _pointerVec.clear();
}

template<typename T, typename W>
std::string SuperVTMwriter3D<T,W>::getName() const
{
  return _name;
}




////////////////////private member functions///////////////////////////////////
template<typename T, typename W>
void SuperVTMwriter3D<T,W>::preambleVTI (const std::string& fullName,
    int x0,int y0,int z0,int x1,int y1,int z1, T originX,T originY,T originZ, T delta)
{
  std::ofstream fout(fullName.c_str(), std::ios::trunc);
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }
  fout << "<?xml version=\"1.0\"?>\n";
  fout << "<VTKFile type=\"ImageData\" version=\"0.1\" "
       << "byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
  fout << "<ImageData WholeExtent=\""
       << x0 <<" "<< x1 <<" "
       << y0 <<" "<< y1 <<" "
       << z0 <<" "<< z1
       << "\" Origin=\"" << originX << " " << originY << " " << originZ
       << "\" Spacing=\"" << delta << " " << delta << " " << delta << "\">\n";
  fout << "<Piece Extent=\""
       << x0 <<" "<< x1 <<" "
       << y0 <<" "<< y1 <<" "
       << z0 <<" "<< z1 <<"\">\n";
  fout << "<PointData>\n";
  fout.close();
}

template<typename T, typename W>
void SuperVTMwriter3D<T,W>::closeVTI(const std::string& fullNamePiece)
{
  std::ofstream fout(fullNamePiece.c_str(), std::ios::app );
  if (!fout) {
    clout << "Error: could not open " << fullNamePiece << std::endl;
  }
  fout << "</ImageData>\n";
  fout << "</VTKFile>\n";
  fout.close();
}

template<typename T, typename W>
void SuperVTMwriter3D<T,W>::preamblePVD(const std::string& fullNamePVD)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::trunc);
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
  fout << "<?xml version=\"1.0\"?>\n";
  fout << "<VTKFile type=\"Collection\" version=\"0.1\" "
       << "byte_order=\"LittleEndian\">\n"
       << "<Collection>\n";
  fout.close();
}

template<typename T, typename W>
void SuperVTMwriter3D<T,W>::closePVD(const std::string& fullNamePVD)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::app );
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
  fout << "</Collection>\n";
  fout << "</VTKFile>\n";
  fout.close();
}

template<typename T, typename W>
void SuperVTMwriter3D<T,W>::preambleVTM(const std::string& fullNamePVD)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::trunc);
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
  fout << "<?xml version=\"1.0\"?>\n";
  fout << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" "
       << "byte_order=\"LittleEndian\">\n"
       << "<vtkMultiBlockDataSet>\n" ;
  fout.close();
}

template<typename T, typename W>
void SuperVTMwriter3D<T,W>::closeVTM(const std::string& fullNamePVD)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::app );
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
  fout << "</vtkMultiBlockDataSet>\n";
  fout << "</VTKFile>\n";
  fout.close();
}

template<typename T, typename W>
void SuperVTMwriter3D<T,W>::dataVTM(int iC, const std::string& fullNamePVD,
                                    const std::string& namePiece)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
  fout << "<Block index=\"" << iC << "\" >\n";
  fout << "<DataSet index= \"0\" " << "file=\"" << namePiece << "\">\n"
       << "</DataSet>\n";
  fout << "</Block>\n";
  fout.close();
}

template<typename T, typename W>
void SuperVTMwriter3D<T,W>::dataPVDmaster(int iT,
    const std::string& fullNamePVDMaster,
    const std::string& namePiece)
{
  std::ofstream fout(fullNamePVDMaster.c_str(), std::ios::in | std::ios::out | std::ios::ate);
  if (fout) {
    fout.seekp(-25,std::ios::end);    // jump -25 form the end of file to overwrite closePVD

    fout << "<DataSet timestep=\"" << iT << "\" "
         << "group=\"\" part=\"\" "
         << "file=\"" << namePiece << "\"/>\n";
    fout.close();
    closePVD(fullNamePVDMaster);
  } else {
    clout << "Error: could not open " << fullNamePVDMaster << std::endl;
  }
}

template<typename T, typename W>
void SuperVTMwriter3D<T,W>::dataArray(const std::string& fullName,
                                      SuperF3D<T,W>& f, int iC, int nx, int ny, int nz)
{
  const char* fileName = fullName.c_str();

  std::ofstream fout( fileName, std::ios::out | std::ios::app );
  if (!fout) {
    clout << "Error: could not open " << fileName << std::endl;
  }

  fout << "<DataArray ";
  if (f.getTargetDim() == 1) {
    fout << "type=\"Float32\" Name=\"" << f.getName() << "\" "
         << "format=\"binary\" encoding=\"base64\">\n";
  } else {
    fout << "type=\"Float32\" Name=\"" << f.getName() << "\" "
         << "format=\"binary\" encoding=\"base64\" "
         << "NumberOfComponents=\"" << f.getTargetDim() <<"\">\n";
  }

  int i[4] = {iC, 0, 0, 0};
  W evaluated[f.getTargetDim()];
  for (int iDim = 0; iDim < f.getTargetDim(); ++iDim) {
    evaluated[iDim] = W();
  }

  size_t numberOfFloats = f.getTargetDim() * (nx+2) * (ny+2) * (nz+2);
  uint32_t binarySize = static_cast<uint32_t>( numberOfFloats*sizeof(float) );

  std::unique_ptr<float[]> streamFloat(new float[numberOfFloats]);    // stack may be too small
  int itter = 0;
  // fill buffer with functor data
  for (i[3] = -1; i[3] < nz+1; ++i[3]) {
    for (i[2] = -1; i[2] < ny+1; ++i[2]) {
      for (i[1] = -1; i[1] < nx+1; ++i[1]) {
        f(evaluated,i);
        for (int iDim = 0; iDim < f.getTargetDim(); ++iDim) {
          streamFloat[itter] = float( evaluated[iDim] );
          itter ++;
        }
      }
    }
  }

  // char buffer for functor data
  const unsigned char* charData = reinterpret_cast<unsigned char*>(streamFloat.get());
  // buffer for compression
  std::unique_ptr<unsigned char[]> comprData(new unsigned char[ binarySize ]);    // stack may be too small

  // compress data (not yet decoded as base64) by zlib
  uLongf sizeCompr = compressBound(binarySize);
  compress2( comprData.get(), &sizeCompr, charData, binarySize, -1);

  // encode prefix to base64 documented in  http://www.earthmodels.org/software/vtk-and-paraview/vtk-file-formats
  Base64Encoder<uint32_t> prefixEncoder(fout, 4);
  uint32_t prefix[4] = {1,binarySize,binarySize,static_cast<uint32_t>(sizeCompr)};
  prefixEncoder.encode(prefix, 4);

  // encode compressed data to base64
  Base64Encoder<unsigned char> dataEncoder( fout, sizeCompr );
  dataEncoder.encode(comprData.get(), sizeCompr);

  fout.close();

  std::ofstream ffout( fileName,  std::ios::out | std::ios::app );
  if (!ffout) {
    clout << "Error: could not open " << fileName << std::endl;
  }
  ffout << "\n</DataArray>\n";
  ffout.close();
}

template<typename T, typename W>
void SuperVTMwriter3D<T,W>::closePiece(const std::string& fullNamePiece)
{
  std::ofstream fout(fullNamePiece.c_str(), std::ios::app );
  if (!fout) {
    clout << "Error: could not open " << fullNamePiece << std::endl;
  }
  fout << "</PointData>\n";
  fout << "</Piece>\n";
  fout.close();
}


}  // namespace olb

#endif
