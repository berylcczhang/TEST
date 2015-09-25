# Makefile for the slug code, v2
.PHONY: all debug clean bayesphot slug bayesphot-debug slug-debug

MACHINE	=
FITS ?= ENABLE_FITS

all: slug bayesphot

debug: slug-debug bayesphot-debug

bayesphot:
	cd slugpy/bayesphot/bayesphot_c && $(MAKE) MACHINE=$(MACHINE)
	@(cp slugpy/bayesphot/bayesphot_c/bayesphot.* slugpy/bayesphot)

slug:
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

bayesphot-debug:
	cd slugpy/bayesphot/bayesphot_c && $(MAKE) debug MACHINE=$(MACHINE)
	@(cp slugpy/bayesphot/bayesphot_c bayesphot.* slugpy/bayesphot)

slug-debug:
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
	cd src/bayesphot && $(MAKE) clean
	@(rm -f slugpy/bayesphot/bayesphot.so)
	@(rm -f slugpy/bayesphot/bayesphot.dylib)
	@(rm -f slugpy/bayesphot/bayesphot.dll)
