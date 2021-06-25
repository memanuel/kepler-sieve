# *************************************************************************************************
# Tools and libraries used
# *************************************************************************************************
# gcc version 10.3.0
# Boost 1.71
# (both gcc and Boost current version installed through apt)
# MariaDB Connector/C version 3.1.11 and Connector/C++ 1.0
# MariaDB server is on version 10.5.10

# Documentation for MariaDB Connector/C++:
# https://mariadb.com/docs/clients/connector-cpp/

# *************************************************************************************************
# Set environment
# *************************************************************************************************
# Include the configuration file
include build/config.mk

# Include the file dependencies
include build/Makefile.dep

# Directory for source and executables; see
# https://stackoverflow.com/questions/2908057/can-i-compile-all-cpp-files-in-src-to-os-in-obj-then-link-to-binary-in
SRC_DIR := src
OBJ_DIR := obj

# Set verbosity flag; uncomment this line to show lists of source files, object files and executables.
# VERBOSE := true

# *************************************************************************************************
# Source, object and executable files
# *************************************************************************************************
# List of C++ source files; filename.cpp.  Includes both executables and intermediates.
SRC_CPP := $(wildcard $(SRC_DIR)/*.cpp)
# Debug listing of SRC_CPP
ifdef VERBOSE
$(info SRC_CPP is $(SRC_CPP) $(NEWLINE))
endif

# List of object files; filename.o
OBJ := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC_CPP))
# Debug listing of object files
ifdef VERBOSE
$(info OBJ is $(OBJ) $(NEWLINE))
endif

# List of executable files is at the top of Makefile.dep!
# Debug listing of executable prefix names
ifdef VERBOSE
$(info EXEC is $(EXEC) $(NEWLINE))
endif

# List of executables files on this platform; e.g. myprogram.x on Linux, myprogram.exe on Windows
EXEC_FILES := $(patsubst %, %.$(EXEC_SUFFIX), $(EXEC))

# *************************************************************************************************
# Convenience targets: all, clean, default
# *************************************************************************************************
# Define target "all" to build all executables
all: $(EXEC)

# The PHONY command tells GNU Make than "clean" is a phony target
.PHONY: clean

# Define target "clean" to remove all the built files
clean:
	rm -f $(OBJ) $(EXEC_FILES)

# Set default goal to "all"
.DEFAULT_GOAL := all

# *************************************************************************************************
# Rules to build object and executable filess
# *************************************************************************************************
# Build object files from source (including executables)
# %.o: %.cpp
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXX_OUT_OBJ) $(CXX_FLAGS) $(INCLUDE) $(NEWLINE)

# Link executables from object files
# $(EXEC): $(OBJ)
$(EXEC):
	$(CXX) $(CXX_OUT_EXE) $(CXX_FLAGS) $(LD_FLAGS) $(LD_LIBS) $(NEWLINE)
