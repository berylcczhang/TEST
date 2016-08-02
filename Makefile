# Makefile for the slug code, v2
.PHONY: all debug clean bayesphot slug bayesphot-debug slug-debug exe lib lib-debug

MACHINE	=
FITS ?= ENABLE_FITS
GSLVERSION ?= 2

all: slug bayesphot lib

exe: slug

debug: slug-debug bayesphot-debug

bayesphot:
	cd slugpy/bayesphot/bayesphot_c && $(MAKE) MACHINE=$(MACHINE)
	@(cp slugpy/bayesphot/bayesphot_c/bayesphot.* slugpy/bayesphot)

slug:
	cd src && $(MAKE) MACHINE=$(MACHINE) FITS=$(FITS) GSLVERSION=$(GSLVERSION)
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
	cd src && $(MAKE) debug MACHINE=$(MACHINE) FITS=$(FITS) GSLVERSION=$(GSLVERSION)
	@(if [ ! -e bin ]; \
	then \
		mkdir bin; \
	fi)
	@(if [ ! -e output ]; \
	then \
		mkdir output; \
	fi)
	@(cp src/slug bin)

lib:
	cd src && $(MAKE) lib MACHINE=$(MACHINE) FITS=$(FITS) GSLVERSION=$(GSLVERSION)

lib-debug:
	cd src && $(MAKE) lib-debug MACHINE=$(MACHINE) FITS=$(FITS) GSLVERSION=$(GSLVERSION)

clean:
	cd src && $(MAKE) clean
	@(if [ ! -e bin ]; \
	then \
		rm -f bin/slug; \
	fi)
	cd slugpy/bayesphot/bayesphot_c && $(MAKE) clean
	@(rm -f slugpy/bayesphot/bayesphot.so)
	@(rm -f slugpy/bayesphot/bayesphot.dylib)
	@(rm -f slugpy/bayesphot/bayesphot.dll)
	@(rm -f src/libslug.so)
	@(rm -f src/libslug.dylib)
	@(rm -f src/libslug.dll)
