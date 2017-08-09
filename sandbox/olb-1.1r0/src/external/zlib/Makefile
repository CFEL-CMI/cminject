# Makefile for zlib
# Copyright (C) 1995-2017 Jean-loup Gailly, Mark Adler
# For conditions of distribution and use, see copyright notice in zlib.h

ROOT    := ../../..

include $(ROOT)/Makefile.inc

CFLAGS=$(filter-out -std=%,$(CXXFLAGS))

AR=$(ARPRG)
ARFLAGS=rc

SRCDIR=
ZINC=

OBJZ = adler32.o crc32.o deflate.o infback.o inffast.o inflate.o inftrees.o trees.o zutil.o
OBJG = compress.o uncompr.o gzclose.o gzlib.o gzread.o gzwrite.o
OBJC = $(OBJZ) $(OBJG)

all: static

static: libz.a

check: test

compile: static

depend:

libz.a: $(OBJC)
	$(AR) $(ARFLAGS) $@ $(OBJC)

adler32.o: $(SRCDIR)adler32.c
	$(CC) $(CFLAGS) $(ZINC) -c -o $@ $(SRCDIR)adler32.c

crc32.o: $(SRCDIR)crc32.c
	$(CC) $(CFLAGS) $(ZINC) -c -o $@ $(SRCDIR)crc32.c

deflate.o: $(SRCDIR)deflate.c
	$(CC) $(CFLAGS) $(ZINC) -c -o $@ $(SRCDIR)deflate.c

infback.o: $(SRCDIR)infback.c
	$(CC) $(CFLAGS) $(ZINC) -c -o $@ $(SRCDIR)infback.c

inffast.o: $(SRCDIR)inffast.c
	$(CC) $(CFLAGS) $(ZINC) -c -o $@ $(SRCDIR)inffast.c

inflate.o: $(SRCDIR)inflate.c
	$(CC) $(CFLAGS) $(ZINC) -c -o $@ $(SRCDIR)inflate.c

inftrees.o: $(SRCDIR)inftrees.c
	$(CC) $(CFLAGS) $(ZINC) -c -o $@ $(SRCDIR)inftrees.c

trees.o: $(SRCDIR)trees.c
	$(CC) $(CFLAGS) $(ZINC) -c -o $@ $(SRCDIR)trees.c

zutil.o: $(SRCDIR)zutil.c
	$(CC) $(CFLAGS) $(ZINC) -c -o $@ $(SRCDIR)zutil.c

compress.o: $(SRCDIR)compress.c
	$(CC) $(CFLAGS) $(ZINC) -c -o $@ $(SRCDIR)compress.c

uncompr.o: $(SRCDIR)uncompr.c
	$(CC) $(CFLAGS) $(ZINC) -c -o $@ $(SRCDIR)uncompr.c

gzclose.o: $(SRCDIR)gzclose.c
	$(CC) $(CFLAGS) $(ZINC) -c -o $@ $(SRCDIR)gzclose.c

gzlib.o: $(SRCDIR)gzlib.c
	$(CC) $(CFLAGS) $(ZINC) -c -o $@ $(SRCDIR)gzlib.c

gzread.o: $(SRCDIR)gzread.c
	$(CC) $(CFLAGS) $(ZINC) -c -o $@ $(SRCDIR)gzread.c

gzwrite.o: $(SRCDIR)gzwrite.c
	$(CC) $(CFLAGS) $(ZINC) -c -o $@ $(SRCDIR)gzwrite.c


# The ranlib in install is needed on NeXTSTEP which checks file times
# ldconfig is for Linux
clean:
	rm -f *.o *.lo *~ \
	   libz.* foo.gz so_locations \
	   _match.s maketree contrib/infback9/*.o
	rm -rf objs
	rm -f *.gcda *.gcno *.gcov
	rm -f contrib/infback9/*.gcda contrib/infback9/*.gcno contrib/infback9/*.gcov

maintainer-clean: distclean
distclean: clean zconf zconf.h.cmakein docs
	rm -f Makefile zlib.pc configure.log
	-@rm -f .DS_Store
	@if [ -f Makefile.in ]; then \
	printf 'all:\n\t-@echo "Please use ./configure first.  Thank you."\n' > Makefile ; \
	printf '\ndistclean:\n\tmake -f Makefile.in distclean\n' >> Makefile ; \
	touch -r $(SRCDIR)Makefile.in Makefile ; fi
	@if [ ! -f zconf.h.in ]; then rm -f zconf.h zconf.h.cmakein ; fi
	@if [ ! -f zlib.3 ]; then rm -f zlib.3.pdf ; fi

tags:
	etags $(SRCDIR)*.[ch]

adler32.o zutil.o: $(SRCDIR)zutil.h $(SRCDIR)zlib.h zconf.h
gzclose.o gzlib.o gzread.o gzwrite.o: $(SRCDIR)zlib.h zconf.h $(SRCDIR)gzguts.h
compress.o example.o minigzip.o uncompr.o: $(SRCDIR)zlib.h zconf.h
crc32.o: $(SRCDIR)zutil.h $(SRCDIR)zlib.h zconf.h $(SRCDIR)crc32.h
deflate.o: $(SRCDIR)deflate.h $(SRCDIR)zutil.h $(SRCDIR)zlib.h zconf.h
infback.o inflate.o: $(SRCDIR)zutil.h $(SRCDIR)zlib.h zconf.h $(SRCDIR)inftrees.h $(SRCDIR)inflate.h $(SRCDIR)inffast.h $(SRCDIR)inffixed.h
inffast.o: $(SRCDIR)zutil.h $(SRCDIR)zlib.h zconf.h $(SRCDIR)inftrees.h $(SRCDIR)inflate.h $(SRCDIR)inffast.h
inftrees.o: $(SRCDIR)zutil.h $(SRCDIR)zlib.h zconf.h $(SRCDIR)inftrees.h
trees.o: $(SRCDIR)deflate.h $(SRCDIR)zutil.h $(SRCDIR)zlib.h zconf.h $(SRCDIR)trees.h
