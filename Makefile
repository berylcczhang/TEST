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
	cd cluster_slug && $(MAKE) all MACHINE=$(MACHINE)
	@(cp cluster_slug/cluster_slug.* slugpy/cluster_slug)

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
	cd cluster_slug && $(MAKE) debug MACHINE=$(MACHINE)
	@(cp cluster_slug/cluster_slug.* slugpy/cluster_slug)

clean:
	cd src && $(MAKE) clean
	@(if [ ! -e bin ]; \
	then \
		rm -f bin/slug; \
	fi)
	cd cluster_slug && $(MAKE) clean
	@(rm -f slugpy/cluster_slug/cluster_slug.so)
	@(rm -f slugpy/cluster_slug/cluster_slug.dylib)
