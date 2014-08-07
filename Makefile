# Makefile for the slug code, v2

.PHONY: all debug exec clean

all:
	cd src && $(MAKE) all
	@(if [ ! -e bin ]; \
	then \
		mkdir bin; \
	fi)
	@(cp src/slug bin)

debug:
	cd src && $(MAKE) debug
	@(if [ ! -e bin ]; \
	then \
		mkdir bin; \
	fi)
	@(cp src/slug bin)

clean:
	cd src && $(MAKE) clean
	@(if [ ! -e bin ]; \
	then \
		mkdir bin; \
	fi)
	@(rm bin/slug)
