# Set flags for debug and non-debug mode
CXXOPTFLAGS     = -O3 -MMD -MP -I $(CXX_INCLUDE_PATH) -I $(BOOST_HOME) -Wall -DNDEBUG
LDOPTFLAGS      = -O3 -lm -lboost_system-mt -lboost_filesystem-mt -L $(LD_LIBRARY_PATH)
CXXDEBFLAGS     = -g -MMD -MP -I $(CXX_INCLUDE_PATH) -I $(BOOST_HOME) -Wall
LDDEBFLAGS      = -g -lm -lboost_system-mt -lboost_filesystem-mt -L $(LD_LIBRARY_PATH)

# Set default compilation flags to non-debug ones; can be overridden below
CXXFLAGS	= $(CXXOPTFLAGS)
LDFLAGS		= $(LDOPTFLAGS)

# Name for standalone executable output
EXENAME		= slug

# Default target
.PHONY: debug exec
exec:	$(EXENAME)

# Include dependencies
-include $(DEPS)

# Pointers to sources
SOURCES		= $(wildcard *.cpp)
OBJECTS		= $(SOURCES:%.cpp=%.o)
DEPS		= $(SOURCES:%.cpp=%.d)

debug: CXXFLAGS = $(CXXDEBFLAGS)
debug: LDFLAGS = $(LDDEBFLAGS)
debug: exec

exec:	$(EXENAME)

$(EXENAME):	$(OBJECTS)
	$(CXX) $(LDFLAGS) -o $(EXENAME) $^

compile:	$(OBJECTS)

clean:
	rm -f $(EXENAME) $(OBJECTS) $(DEPS)
