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
#  of the License, or (at your option) any later version.
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

include Makefile.inc

###########################################################################
## all

all: lib

###########################################################################
## compile

compile:
	@for i in $(SUBDIRS); \
	do \
	  (cd $$i; \
           echo "-------------------------------------------------------------"; \
           echo "-- Entering $$i (compile)"; \
           $(MAKE) compile; \
           echo "-- Leaving $$i (compile)"; \
           echo "-------------------------------------------------------------") \
	done

###########################################################################
## depend

depend:
	@for i in $(SUBDIRS); \
	do \
	  (cd $$i; \
           echo "-------------------------------------------------------------"; \
           echo "-- Entering $$i (depend)"; \
           $(MAKE) depend; \
           echo "-- Leaving $$i (depend)"; \
           echo "-------------------------------------------------------------") \
	done

###########################################################################
## clean

clean:
	@rm -f *~ core
	@for i in $(BUILDTYPEDIRS); \
	do \
	  (cd $$i; \
           echo "-------------------------------------------------------------"; \
           echo "-- Clean object, dependencies and library files in $$i"; \
           rm -f lib/*.a; \
           rm -f dep/*.d; \
           rm -f obj/*.o; \
           echo "-------------------------------------------------------------") \
	done

cleanbuild:
	@echo "-------------------------------------------------------------";
	@echo "-- Clean object, dependencies and library files in ´"$(BUILDTYPE)"'";
	@echo "-------------------------------------------------------------";
	@rm -f $(LIBDIR)/*.a
	@rm -f $(DEPENDDIR)/*.d
	@rm -f $(OBJDIR)/*.o
	@for i in $(SUBDIRS); \
	do \
	  (cd $$i; \
           echo "-------------------------------------------------------------"; \
           echo "-- Entering $$i (clean)"; \
           $(MAKE) clean; \
           echo "-- Leaving $$i (clean)"; \
           echo "-------------------------------------------------------------") \
	done

cleansamples:
	@for i in $(EXAMPLEDIRS); \
	do \
	  (cd $$i; \
           echo "-------------------------------------------------------------"; \
           echo "-- Entering $$i (clean)"; \
           $(MAKE) clean; \
           echo "-- Leaving $$i (clean)"; \
           echo "-------------------------------------------------------------") \
	done

###########################################################################
## lib

lib: depend compile $(LIBDIR)/lib$(LIB).a $(LIBDIR)/libz.a
$(LIBDIR)/lib$(LIB).a: $(wildcard $(OBJDIR)/*.o)
	@echo Build lib$(LIB).a:
	@$(ARPRG) -rusv $(LIBDIR)/lib$(LIB).a $(OBJDIR)/*.o

$(LIBDIR)/libz.a: src/external/zlib/libz.a
	@cp src/external/zlib/libz.a $(LIBDIR)/libz.a

###########################################################################
## examples

samples:
	@for i in $(EXAMPLEDIRS); \
	do \
	  (cd $$i; \
           echo "-------------------------------------------------------------"; \
           echo "-- Entering $$i (make examples)"; \
           $(MAKE) all; \
           echo "-- Leaving $$i (make examples)"; \
           echo "-------------------------------------------------------------") \
	done

###########################################################################
## user guide documentation

userguide:
	@cd doc/olb-ug-latest/; \
	latexmk -pdf -silent -f olb-ug.tex

###########################################################################
## doxygen documentation

doxygen:
	doxygen doc/DoxygenConfig

###########################################################################
## makefile options file

Makefile.inc:
	cp Makefile.git.inc Makefile.inc

###########################################################################
###########################################################################
