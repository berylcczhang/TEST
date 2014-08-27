# Makefile for the slug code, v2
.PHONY: all debug exec clean

MACHINE	=
FITS ?= ENABLE_FITS

all:
	cd src && $(MAKE) all MACHINE=$(MACHINE) FITS=$(FITS)
	@(if [ ! -e bin ]; \
	then \
		mkdir bin; \
	fi)
	@(if [ ! -e output ]; \
        then \
                mkdir output; \
        fi)
	@(cp src/slug bin)

debug:
	cd src && $(MAKE) debug MACHINE=$(MACHINE) FITS=$(FITS)
	@(if [ ! -e bin ]; \
	then \
		mkdir bin; \
	fi)
	@(if [ ! -e output ]; \
	then \
		mkdir output; \
	fi)
	@(cp src/slug bin)

clean:
	cd src && $(MAKE) clean
	@(if [ ! -e bin ]; \
	then \
		rm -f bin/slug; \
	fi)