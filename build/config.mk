# *************************************************************************************************
# Michael S. Emanuel
# 2021-06-24
# *************************************************************************************************

# *************************************************************************************************
# Newline and tab characters for legibile printing to screen
# *************************************************************************************************
# newline character
define NEWLINE


endef

# tab character
TAB := $(shell printf '\011')

# *************************************************************************************************
# Check the operating system.  Supported OS Linux, Windows
# *************************************************************************************************
# https://stackoverflow.com/questions/714100/os-detecting-makefile

# Is the OS Linux?
ifeq ($(shell uname), Linux)
	detected_OS := Linux
	OBJECT_SUFFIX := o
	EXEC_SUFFIX := x
endif

# Is the OS Windows?
ifeq ($(OS), Windows_NT)
    detected_OS := Windows
	OBJECT_SUFFIX := obj
	EXEC_SUFFIX := exe
endif

# Show the operating system
$(info OS is $(detected_OS))

# If the OS is not Windows or Linux, quit early: unsupported.
ifndef detected_OS
	$(error Unspported operating system!)
endif

# *************************************************************************************************
# Make settings, environment variables, initial status message
# *************************************************************************************************
# Get the number of jobs for make from the environment; default to all of them
ifndef MAKE_JOBS
	NUMPROC := $(shell grep -c ^processor /proc/cpuinfo)
	MAKE_JOBS := $(shell echo $(NUMPROC)/1 + 0 | bc)
endif

# Make settings: warn on unset variables, use parallel processing
MAKEFLAGS=--warn-undefined-variables --jobs=$(MAKE_JOBS)

# Show the number of jobs
$(info Running make with up to $(MAKE_JOBS) parallel jobs.$(NEWLINE))

# *************************************************************************************************
# Load required modules (e.g. boost)
# *************************************************************************************************

# TODO figure this out

# *************************************************************************************************
# Compiler settings
# *************************************************************************************************
# C++ compiler
CXX=g++

# Compilation flags
CXX_FLAGS= \$(NEWLINE) $(TAB) -Wall -std=c++20 -O3

# Output command for object files
# Note $< is a shorthand for the first dependency
CXX_OUT_OBJ = -c $< -o $@

# Output command for executables
# Note $@ is a shorthand for the file to be built
# Note $^ is a shorthand for all the dependencies
CXX_OUT_EXE = -o $@.$(EXEC_SUFFIX) \$(NEWLINE) $(TAB) $^

# *************************************************************************************************
# INCLUDE: Include paths for additional software libraries 
# *************************************************************************************************
# Root directory for manually installed software libraries; 
# found in environment variable SOFTWARE_LIBRARY_DIR
ifdef SOFTWARE_LIBRARY_DIR
	INCLUDE_USR := $(addprefix -I,$(SOFTWARE_LIBRARY_DIR))
endif	

# Boost flags are read using environment variable BOOST_DIR
ifdef BOOST_DIR
	INCLUDE_BOOST=$(addprefix -I,$(BOOST_DIR))
endif

# *************************************************************************************************
# LD: Linker flags for additional libraries
# *************************************************************************************************
# Directory with library files (.a) for additional software libraries
LD_FLAGS_USR := $(addprefix -L, $(addsuffix /lib, $(SOFTWARE_LIBRARY_DIR)))

# fmt (string formattng and printing)
LD_FMT := -lfmt

# Boost
LD_LIB_BOOST := -lboost_program_options

# Math library
LD_LIB_MATH := -lm

# BLAS/LAPACK (linear algebra)
LD_LIB_LAPACK := -llapack -lblas

# GNU Scientific Library (numerical computing)
LD_LIB_GSL := -lgsl -lgslcblas

# MariaDB (database connectivity)
LD_LIB_MARIADB := -lmariadbcpp

# *************************************************************************************************
# Generate command line arguments 
# -I: additional include directories
# -L: additional library directories
# -l: additional libraries in library search path named lib<library_name>.a
# *************************************************************************************************
# Additional include directories -I
# Initialize INCLUDE to an empty string
INCLUDE := 

# Include user software libraries if provided
ifdef SOFTWARE_LIBRARY_DIR
	INCLUDE := $(INCLUDE) \$(NEWLINE) $(TAB) $(INCLUDE_USR)
endif

# Only include boost if it was supplied manually
# Note: on Ubuntu, if Boost was installed using b2, this is not necessary
ifdef BOOST_DIR
	INCLUDE := $(INCLUDE) \$(NEWLINE) $(TAB) $(INCLUDE_BOOST)
endif

# Additional library directories -L
# LD_FLAGS := \
# \$(NEWLINE) $(TAB) $(LD_FLAGS_USR) 
LD_FLAGS := $(LD_FLAGS_USR) 

# Additional libraries -l
# LD_LIBS := \
# \$(NEWLINE) $(TAB) $(LD_LIB_MARIADB)
LD_LIBS := $(LD_FMT) $(LD_LIB_BOOST) $(LD_LIB_GSL) $(LD_LIB_MARIADB)
