# Main makefile for slug2

# Did the user tell us to use a particular machine? If so, use that.
ifdef MACHINE
     include Make.mach.$(MACHINE)
else
     # Machine not specified, so try to guess
     UNAME		= $(shell uname)
     UNAMEN		= $(shell uname -n)

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
CXX		= $(MACH_CXX)

# Set optimization mode flags
CXXOPTFLAGS	= $(MACH_CXXOPTFLAGS) $(MACH_C11FLAG) -DNDEBUG \
	-DBOOST_DISABLE_ASSERTS -DHAVE_INLINE \
	-MMD -MP
LDOPTFLAGS	= $(MACH_LDOPTFLAGS) $(MACH_CXXFLAG)

# Set debug mode flags
CXXDEBFLAGS     = $(MACH_CXXDEBFLAGS) $(MACH_C11FLAG) -MMD -MP
LDDEBFLAGS	= $(MACH_LDDEBFLAGS) $(MACH_CXXFLAG)

# Read any user overrides
-include Make.config.override

# Include flags
ifdef C_INCLUDE_PATH
     INCFLAGS += -I$(C_INCLUDE_PATH)
endif
ifdef CXX_INCLUDE_PATH
     INCFLAGS += -I$(CXX_INCLUDE_PATH)
endif
ifdef BOOST_HDR_PATH
     INCFLAGS += -I$(BOOST_HDR_PATH)
endif
ifdef GSL_HDR_PATH
     INCFLAGS += -I$(GSL_HDR_PATH)
endif

# Link flags
LDLIBFLAGS      = -lgsl -lgslcblas -lboost_system$(MACH_BOOST_TAG) \
	-lboost_filesystem$(MACH_BOOST_TAG) \
	-lboost_regex$(MACH_BOOST_TAG)
ifdef LD_LIBRARY_PATH
     LDLIBFLAGS += -L$(LD_LIBRARY_PATH)
endif
ifdef LIBRARY_PATH
     LDLIBFLAGS += -L$(LIBRARY_PATH)
endif
ifdef BOOST_LIB_PATH
     LDLIBFLAGS += -L$(BOOST_LIB_PATH)
endif
ifdef GSL_LIB_PATH
     LDLIBFLAGS += -L$(GSL_LIB_PATH)
endif

# Turn FITS on or off
ifeq ($(FITS), ENABLE_FITS)
     LDLIBFLAGS += -lcfitsio
     DEFINES    += -DENABLE_FITS
endif
ifdef FITS_LIB_PATH
     LDLIBFLAGS += -L$(FITS_LIB_PATH)
endif

# Specify which version of the GSL we're using
DEFINES += -DGSLVERSION=$(GSLVERSION)

CXXFLAGS +=  $(INCFLAGS) $(DEFINES)
LDFLAGS  +=  $(LDLIBFLAGS)

# Name for executable
EXENAME		= slug

# Name for compiled library
LIBNAME		= lib$(EXENAME)$(LIB_EXTENSION)

# Pointers to sources
SOURCES		= $(wildcard *.cpp)
OBJECTS		= $(SOURCES:%.cpp=%.o)
DEPS		= $(SOURCES:%.cpp=%.d)

# Default target
.PHONY: all exe debug lib lib-debug compile clean


exe: CXXFLAGS += $(CXXOPTFLAGS)
exe: LDFLAGS  += $(LDOPTFLAGS)
exe: $(EXENAME) 

debug: CXXFLAGS += $(CXXDEBFLAGS)
debug: LDFLAGS  += $(LDDEBFLAGS)
debug: $(EXENAME)

lib: CXXFLAGS += $(CXXOPTFLAGS) $(CLIBFLAGS)
lib: LDFLAGS  += $(LDOPTFLAGS)
lib: $(LIBNAME)

lib-debug: CXXFLAGS += $(CXXDEBFLAGS) $(CLIBFLAGS)
lib-debug: LDFLAGS  += $(LDDEBFLAGS)
lib-debug: $(LIBNAME)

all: exe
all: lib

# Include dependencies
-include $(DEPS)

$(EXENAME): $(OBJECTS)
	$(CXX) -o $(EXENAME) $^ $(LDFLAGS) $(LDLIBFLAGS)

$(LIBNAME): $(OBJECTS)
	$(CXX) -o $(LIBNAME) $^ $(LDFLAGS) $(LDLIBFLAGS) $(DYNLIBFLAG)

compile:	$(OBJECTS)

clean:
	rm -f $(EXENAME) $(OBJECTS) $(DEPS) $(LIBNAME)
