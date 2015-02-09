# Makefile for the c part of bayesphot

# Did the user tell us to use a particular machine? If so, use that.
ifdef MACHINE
     include Make.mach.$(MACHINE)
else
     # Machine not specified, so try to guess
     UNAME              = $(shell uname)
     UNAMEN             = $(shell uname -n)

     # Do we have a makefile that matches the machine name? If so, use
     # that. If not, use a generic makefile depending on the OS.
     ifeq ($(UNAMEN), hyades.ucsc.edu)
          include Make.mach.ucsc-hyades
     else ifeq ($(UNAMEN), cosma-e)
          include Make.mach.icc-cosma4
     else ifeq ($(UNAME), Linux)
          include Make.mach.linux-gnu
     else ifeq ($(UNAME), Darwin)
          include Make.mach.darwin
     else
          $(info Cannot detect system type. Suggest you specify MACHINE= manually.)
          include Make.mach.generic
     endif
endif

# Set compiler
CC		= $(MACH_CC)

# Set optimization mode flags
COPTFLAGS	= $(MACH_COPTFLAGS) $(MACH_OPENMPFLAGS) -DNDEBUG \
	-DHAVE_INLINE -MMD -MP -std=c99
LDOPTFLAGS	= $(MACH_LDOPTFLAGS) $(MACH_OPENMPFLAGS)

# Set debug mode flags
CDEBFLAGS	= $(MACH_CDEBFLAGS) $(MACH_OPENMPFLAGS) -MMD -MP -std=c99
LDDEBFLAGS	= $(MACH_LDDEBFLAGS) $(MACH_OPENMPFLAGS)

# Uncomment for diagnostic mode
#COPTFLAGS	+= -DDIAGNOSTIC
#CDEBFLAGS	+= -DDIAGNOSTIC

# Include flags
ifdef C_INCLUDE_PATH
     INCFLAGS += -I$(C_INCLUDE_PATH)
endif
ifdef GSL_HDR_PATH
     INCFLAGS += -I$(GSL_HDR_PATH)
endif

# Link flags
LDLIBFLAGS      = -lgsl
ifdef LD_LIBRARY_PATH
     LDLIBFLAGS += -L$(LD_LIBRARY_PATH)
endif
ifdef GSL_LIB_PATH
     LDLIBFLAGS += -L$(GSL_LIB_PATH)
endif

# Read any user overrides
-include Make.config.override

# Name for compiled library
LIBNAME		= bayesphot$(LIB_EXTENSION)

# Pointers to sources
SOURCES		= $(wildcard *.c)
OBJECTS		= $(SOURCES:%.c=%.o)
DEPS		= $(SOURCES:%.c=%.d)

# Default target
.PHONY: all debug
all: CFLAGS += $(COPTFLAGS) $(CLIBFLAGS) $(INCFLAGS) $(DEFINES)
all: LDFLAGS += $(LDOPTFLAGS) $(LDLIBFLAGS) $(DYNLIBFLAG)
all: $(LIBNAME)

# Include dependencies
-include $(DEPS)

debug: CFLAGS += $(CDEBFLAGS) $(CLIBFLAGS) $(INCFLAGS) $(DEFINES)
debug: LDFLAGS += $(LDDEBFLAGS) $(LDLIBFLAGS) $(DYNLIBFLAG)
debug: $(LIBNAME)

$(LIBNAME):	$(OBJECTS)
	$(CC) -o $(LIBNAME) $^ $(LDFLAGS)

compile:	$(OBJECTS)

clean:
	rm -f $(LIBNAME) $(OBJECTS) $(DEPS)
