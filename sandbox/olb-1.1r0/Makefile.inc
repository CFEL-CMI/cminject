#  This file is part of the OpenLB library
#
#  Copyright (C) 2007 Mathias Krause
#  E-mail contact: info@openlb.net
#  The most recent release of OpenLB can be downloaded at
#  <http://www.openlb.net/>
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public 
#  License along with this program; if not, write to the Free 
#  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
#  Boston, MA  02110-1301, USA.

###########################################################################
###########################################################################
## DEFINITIONS TO BE CHANGED

CXX             := g++
#CXX             := icpc -D__aligned__=ignored
#CXX             := mpiCC
#CXX             := mpic++

CC              := gcc                                          # necessary for zlib, for Intel use icc

OPTIM           := -O3 -Wall -march=native -mtune=native        # for gcc
#OPTIM           := -O3 -Wall -xHost                            # for Intel compiler
DEBUG           := -g -DOLB_DEBUG

CXXFLAGS        := $(OPTIM)
#CXXFLAGS        := $(DEBUG)

# to enable std::shared_ptr in functor arithmetik
# works in gcc 4.3 and later, source https://gcc.gnu.org/projects/cxx0x.html
CXXFLAGS        += -std=c++0x
# works in gcc 4.7 and later (recommended)
#CXXFLAGS        += -std=c++11

#CXXFLAGS        += -fdiagnostics-color=auto
#CXXFLAGS        += -std=gnu++14

ARPRG           := ar
#ARPRG           := xiar                  # mandatory for intel compiler

LDFLAGS         :=

PARALLEL_MODE   := OFF
#PARALLEL_MODE   := MPI
#PARALLEL_MODE   := OMP
#PARALLEL_MODE   := HYBRID

MPIFLAGS        :=
OMPFLAGS        := -fopenmp

#BUILDTYPE       := precompiled
BUILDTYPE       := generic


###########################################################################
## conditional settings

ifeq ($(BUILDTYPE), precompiled)
   CXXFLAGS := -DOLB_PRECOMPILED $(CXXFLAGS)
endif

ifeq ($(PARALLEL_MODE), MPI)
   CXXFLAGS := -DPARALLEL_MODE_MPI $(MPIFLAGS) $(CXXFLAGS)
endif

ifeq ($(PARALLEL_MODE), OMP)
   CXXFLAGS := -DPARALLEL_MODE_OMP $(OMPFLAGS) $(CXXFLAGS)
   LDFLAGS  := $(OMPFLAGS) $(LDFLAGS)
endif

ifeq ($(PARALLEL_MODE), HYBRID)
   CXXFLAGS := -DPARALLEL_MODE_OMP -DPARALLEL_MODE_MPI $(OMPFLAGS) $(MPIFLAGS) $(CXXFLAGS)
   LDFLAGS  := $(OMPFLAGS) $(LDFLAGS)
endif

###########################################################################
## defines shell

SHELL           := /bin/sh

###########################################################################
## dependencies, object, library directory and library name

DEPENDDIR       := build/$(BUILDTYPE)/dep
OBJDIR          := build/$(BUILDTYPE)/obj
LIBDIR          := build/$(BUILDTYPE)/lib
LIB             := olb
LIBS            := -l$(LIB) -lz

###########################################################################
## search directories

SUBDIRS         := src/boundary \
                   src/communication \
                   src/dynamics \
                   src/core \
                   src/geometry \
                   src/external/tinyxml \
                   src/external/zlib \
                   src/functors \
                   src/functors/indicator \
                   src/io \
                   src/particles \
                   src/particles/forces \
                   src/particles/boundaries \
                   src/utilities

EXAMPLEDIRS     := examples/aorta3d \
                   examples/bifurcation3d/eulerEuler \
                   examples/bifurcation3d/eulerLagrange \
                   examples/bstep2d \
                   examples/bstep3d \
                   examples/cavity2d/sequential \
                   examples/cavity2d/parallel \
                   examples/cavity3d/sequential \
                   examples/cavity3d/parallel \
                   examples/cylinder2d \
                   examples/cylinder3d \
                   examples/multiComponent2d \
                   examples/multiComponent3d \
                   examples/nozzle3d \
                   examples/phaseSeparation2d \
                   examples/phaseSeparation3d \
                   examples/poiseuille2d/bgk \
                   examples/poiseuille2d/forced \
                   examples/poiseuille2d/mrt \
                   examples/powerLaw2d \
                   examples/tgv3d \
                   examples/thermal2d \
                   examples/thermal3d \
                   examples/venturi3d

INCLUDEDIRS     := src \
                   src/ \
                   src/external \
                   src/external/zlib
			
BUILDTYPEDIRS   := build/precompiled \
                   build/generic

SUBDIRSH        := $(foreach d,$(SUBDIRS),$(ROOT)/$(d))
IDIR            := $(foreach d,$(INCLUDEDIRS),-I$(ROOT)/$(d))

###########################################################################
###########################################################################
