# Makefile for the slug code, v2
.PHONY: all debug clean

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
	cd bayesphot && $(MAKE) all MACHINE=$(MACHINE)
	@(cp bayesphot/bayesphot.* slugpy/bayesphot)

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
	cd bayesphot && $(MAKE) debug MACHINE=$(MACHINE)
	@(cp bayesphot/bayesphot.* slugpy/bayesphot)

clean:
	cd src && $(MAKE) clean
	@(if [ ! -e bin ]; \
	then \
		rm -f bin/slug; \
	fi)
	cd bayesphot && $(MAKE) clean
	@(rm -f slugpy/bayesphot/bayesphot.so)
	@(rm -f slugpy/bayesphot/bayesphot.dylib)
	@(rm -f slugpy/bayesphot/bayesphot.dll)
