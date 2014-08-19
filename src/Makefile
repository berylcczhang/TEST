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
LDOPTFLAGS	= $(MACH_LDOPTFLAGS)

# Set debug mode flags
CXXDEBFLAGS     = $(MACH_CXXDEBFLAGS) $(MACH_C11FLAG) -MMD -MP
LDDEBFLAGS	= $(MACH_LDDEBFLAGS)

# Include flags
ifdef CXX_INCLUDE_PATH
     INCFLAGS += -I $(CXX_INCLUDE_PATH)
endif
ifdef BOOST_HDR_PATH
     INCFLAGS += -I $(BOOST_HDR_PATH)
endif
ifdef GSL_HDR_PATH
     INCFLAGS += -I $(GSL_HDR_PATH)
endif

# Link flags
LDLIBFLAGS      = -lgsl -lgslcblas -lboost_system$(MACH_BOOST_TAG) \
	-lboost_filesystem$(MACH_BOOST_TAG) \
	-lboost_regex$(MACH_BOOST_TAG)
ifdef LD_LIBRARY_PATH
     LDLIBFLAGS += -L $(LD_LIBRARY_PATH)
endif
ifdef BOOST_LIB_PATH
     LDLIBFLAGS += -L $(BOOST_LIB_PATH)
endif
ifdef GSL_LIB_PATH
     LDLIBFLAGS += -L $(GSL_LIB_PATH)
endif

# Turn FITS on or off
ifeq ($(FITS), ENABLE_FITS)
     LDLIBFLAGS += -lcfitsio
     DEFINES    += -DENABLE_FITS
endif

# Read any user overrides
-include Make.config.override

# Name for executable
EXENAME		= slug

# Pointers to sources
SOURCES		= $(wildcard *.cpp)
OBJECTS		= $(SOURCES:%.cpp=%.o)
DEPS		= $(SOURCES:%.cpp=%.d)

# Default target
.PHONY: all debug
all: CXXFLAGS += $(CXXOPTFLAGS) $(INCFLAGS) $(DEFINES)
all: LDFLAGS += $(LDOPTFLAGS) $(LDLIBFLAGS)
all: $(EXENAME)

# Include dependencies
-include $(DEPS)

debug: CXXFLAGS += $(CXXDEBFLAGS) $(INCFLAGS) $(DEFINES)
debug: LDFLAGS += $(LDDEBFLAGS) $(LDLIBFLAGS)
debug: $(EXENAME)

$(EXENAME):	$(OBJECTS)
	$(CXX) $(LDFLAGS) -o $(EXENAME) $^

compile:	$(OBJECTS)

clean:
	rm -f $(EXENAME) $(OBJECTS) $(DEPS)
