# Makefile for the slug code, v2
.PHONY: all debug clean bayesphot slug bayesphot-debug slug-debug exe lib lib-debug libstatic libstatic-debug tools tools-debug

MACHINE	=
FITS ?= ENABLE_FITS
GSLVERSION ?= 2
MPI ?= NO_MPI

exe: slug bayesphot

all: slug bayesphot lib libstatic

debug: slug-debug bayesphot-debug

tools: write_isochrone

tools-debug: write_isochrone-debug

slug:
	cd src && $(MAKE) MACHINE=$(MACHINE) FITS=$(FITS) GSLVERSION=$(GSLVERSION) MPI=$(MPI)
	@(if [ ! -e bin ]; \
	then \
		mkdir bin; \
	fi)
	@(if [ ! -e output ]; \
        then \
                mkdir output; \
        fi)
	@(cp src/slug bin)

bayesphot:
	cd slugpy/bayesphot/bayesphot_c && $(MAKE) MACHINE=$(MACHINE)
	@(cp slugpy/bayesphot/bayesphot_c/bayesphot.* slugpy/bayesphot)

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
	cd src && $(MAKE) lib MACHINE=$(MACHINE) FITS=$(FITS) GSLVERSION=$(GSLVERSION) MPI=$(MPI)

lib-debug:
	cd src && $(MAKE) lib-debug MACHINE=$(MACHINE) FITS=$(FITS) GSLVERSION=$(GSLVERSION) MPI=$(MPI)

libstatic:
	cd src && $(MAKE) libstatic MACHINE=$(MACHINE) FITS=$(FITS) GSLVERSION=$(GSLVERSION) MPI=$(MPI)

libstatic-debug:
	cd src && $(MAKE) libstatic-debug MACHINE=$(MACHINE) FITS=$(FITS) GSLVERSION=$(GSLVERSION) MPI=$(MPI)

write_isochrone:
	cd tools/c/write_isochrone && $(MAKE) MACHINE=$(MACHINE) FITS=$(FITS) GSLVERSION=$(GSLVERSION)
	@(if [ ! -e bin ]; \
	then \
		mkdir bin; \
	fi)
	@(cp tools/c/write_isochrone/write_isochrone bin)

write_isochrone-debug:
	cd tools/c/write_isochrone && $(MAKE) debug MACHINE=$(MACHINE) FITS=$(FITS) GSLVERSION=$(GSLVERSION)
	@(if [ ! -e bin ]; \
	then \
		mkdir bin; \
	fi)
	@(cp tools/c/write_isochrone/write_isochrone bin)

clean:
	cd src && $(MAKE) clean
	@(if [ ! -e bin ]; \
	then \
		rm -f bin/slug; \
		rm -f bin/write_isochrone; \
	fi)
	cd slugpy/bayesphot/bayesphot_c && $(MAKE) clean
	cd tools/c/write_isochrone && $(MAKE) clean
	@(rm -f slugpy/bayesphot/bayesphot.so)
	@(rm -f slugpy/bayesphot/bayesphot.dylib)
	@(rm -f slugpy/bayesphot/bayesphot.dll)
	@(rm -f src/libslug.so)
	@(rm -f src/libslug.dylib)
	@(rm -f src/libslug.dll)
	@(rm -f src/libslug.a)
	@(rm -f src/libslug.lib)
