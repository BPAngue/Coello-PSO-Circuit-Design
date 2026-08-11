########################################################
# CINVESTAV - IPN | PSO Combinational Circuit Design
########################################################

CXX      := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -pthread
LDFLAGS  := -pthread

BUILD_DIR := build

########################################################
# Shared source files
########################################################

CORE_SRCS := circuits.cpp matrixpso.cpp random.cpp statistics.cpp
CORE_OBJS := $(CORE_SRCS:%.cpp=$(BUILD_DIR)/%.o)

########################################################
# Targets
########################################################

.PHONY: all clean

all: psomatrixcircuit.exe multithreaded.exe

########################################################
# Create build directory
########################################################

$(BUILD_DIR):
	if not exist "$(BUILD_DIR)" mkdir "$(BUILD_DIR)"

########################################################
# Original single-threaded executable
########################################################

psomatrixcircuit.exe: $(CORE_OBJS) $(BUILD_DIR)/psomatrixcircuit.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

########################################################
# Compile psomatrixcircuit.cpp for parallel version
#
# Rename main() so it does not conflict with the
# main() in multithreaded.cpp.
########################################################

$(BUILD_DIR)/psomatrixcircuit_lib.o: psomatrixcircuit.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -Dmain=pso_original_main_unused -c -o $@ $<

########################################################
# Multithreaded executable
########################################################

multithreaded.exe: \
	$(CORE_OBJS) \
	$(BUILD_DIR)/psomatrixcircuit_lib.o \
	$(BUILD_DIR)/multithreaded.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

########################################################
# Compile .cpp -> build/.o
########################################################

$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

########################################################
# Clean build artifacts
########################################################

clean:
	if exist "$(BUILD_DIR)" rmdir /S /Q "$(BUILD_DIR)"
	del /Q psomatrixcircuit.exe multithreaded.exe 2>NUL