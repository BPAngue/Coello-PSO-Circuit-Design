CXX      := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -pthread
LDFLAGS  := -pthread

CORE_SRCS := circuits.cpp matrixpso.cpp random.cpp statistics.cpp
CORE_OBJS := $(CORE_SRCS:.cpp=.o)

.PHONY: all clean

all: psomatrixcircuit multithreaded

########################################################
# Original single-threaded executable
########################################################

psomatrixcircuit: $(CORE_OBJS) psomatrixcircuit.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

########################################################
# Compile psomatrixcircuit.cpp for the parallel program
#
# This keeps all Swarm implementations but renames
# the original main() so it does not conflict with
# multithreaded.cpp's main().
########################################################

psomatrixcircuit_lib.o: psomatrixcircuit.cpp
	$(CXX) $(CXXFLAGS) -Dmain=pso_original_main_unused -c -o $@ $<

########################################################
# Multithreaded executable
########################################################

multithreaded: $(CORE_OBJS) psomatrixcircuit_lib.o multithreaded.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

########################################################
# Normal .cpp -> .o rule
########################################################

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

########################################################
# Clean
########################################################

clean:
	del /Q *.o psomatrixcircuit.exe multithreaded.exe *.csv 2>NUL