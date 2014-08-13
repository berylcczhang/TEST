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

# Set optimization mode flags
CXXOPTFLAGS	= -c $(MACH_CXXOPTFLAGS) $(MACH_C11FLAG) -DNDEBUG \
	-DBOOST_DISABLE_ASSERTS -DHAVE_INLINE \
	-MMD -MP -Wall
LDOPTFLAGS	= $(MACH_LDOPTFLAGS)

# Set debug mode flags
CXXDEBFLAGS     = -c $(MACH_CXXDEBFLAGS) $(MACH_C11FLAG) -MMD -MP -Wall
LDDEBFLAGS	= $(MACH_LDDEBFLAGS)

# Read any user overrides
-include Make.config.override

# Include and load flags
INCFLAGS	= -I $(CXX_INCLUDE_PATH) -I $(BOOST_HDR_PATH) -I $(GSL_HDR_PATH)
LDLIBFLAGS	= -lm -lgsl -lboost_system-mt -lboost_filesystem-mt \
	-lboost_regex-mt -L $(LD_LIBRARY_PATH) $(BOOST_LIB_PATH) \
	$(GSL_LIB_PATH)

# Name for executable
EXENAME		= slug

# Pointers to sources
SOURCES		= $(wildcard *.cpp)
OBJECTS		= $(SOURCES:%.cpp=%.o)
DEPS		= $(SOURCES:%.cpp=%.d)

# Default target
.PHONY: all debug
all: CXXFLAGS += $(CXXOPTFLAGS) $(INCFLAGS)
all: LDFLAGS += $(LDOPTFLAGS) $(LDLIBFLAGS)
all: $(EXENAME)

# Include dependencies
-include $(DEPS)

debug: CXXFLAGS += $(CXXDEBFLAGS) $(INCFLAGS)
debug: LDFLAGS += $(LDDEBFLAGS) $(LDLIBFLAGS)
debug: $(EXENAME)

$(EXENAME):	$(OBJECTS)
	$(CXX) $(LDFLAGS) -o $(EXENAME) $^

compile:	$(OBJECTS)

clean:
	rm -f $(EXENAME) $(OBJECTS) $(DEPS)
