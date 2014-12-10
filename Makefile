# Makefile for bayesphot
.PHONY: all debug clean

MACHINE =

all:
	cd csrc && $(MAKE) all MACHINE=$(MACHINE)
	@(cp csrc/bayesphot.* bayesphot)

debug:
	cd csrc && $(MAKE) debug MACHINE=$(MACHINE)
	@(cp csrc/bayesphot.* bayesphot)

clean:
	cd csrc && $(MAKE) clean
	@(rm -f bayesphot/bayesphot.so)
	@(rm -f bayesphot/bayesphot.dylib)
	@(rm -f bayesphot/bayesphot.dll)
