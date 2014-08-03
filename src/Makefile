# Set flags for debug and non-debug mode
CXXOPTFLAGS     = -O3 -ffast-math -MMD -MP -I $(CXX_INCLUDE_PATH) -I $(BOOST_HOME) -Wall -std=c++11 -DNDEBUG -DBOOST_DISABLE_ASSERTS -DHAVE_INLINE
LDOPTFLAGS      = -O4 -ffast-math -lm -lgsl -lboost_system-mt -lboost_filesystem-mt -lboost_regex-mt -L $(LD_LIBRARY_PATH)
CXXDEBFLAGS     = -g -MMD -MP -I $(CXX_INCLUDE_PATH) -I $(BOOST_HOME) -Wall -std=c++11
LDDEBFLAGS      = -g -lm -lgsl -lboost_system-mt -lboost_filesystem-mt -lboost_regex-mt -L $(LD_LIBRARY_PATH)

# Name for executable
EXENAME		= slug

# Pointers to sources
SOURCES		= $(wildcard *.cpp)
OBJECTS		= $(SOURCES:%.cpp=%.o)
DEPS		= $(SOURCES:%.cpp=%.d)

# Default target
.PHONY: all debug exec
all: CXXFLAGS += $(CXXOPTFLAGS)
all: LDFLAGS += $(LDOPTFLAGS)
all:	exec

# Include dependencies
-include $(DEPS)

$(EXENAME):	all

debug: CXXFLAGS += $(CXXDEBFLAGS)
debug: LDFLAGS += $(LDDEBFLAGS)
debug: exec

exec:	$(OBJECTS)
	$(CXX) $(LDFLAGS) -o $(EXENAME) $^

compile:	$(OBJECTS)

clean:
	rm -f $(EXENAME) $(OBJECTS) $(DEPS)
