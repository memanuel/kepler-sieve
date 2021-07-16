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

# Set verbosity flag; uncomment this line to show lists of source files, object files and executables.
# VERBOSE := 1

# *************************************************************************************************
# Compiler settings
# *************************************************************************************************
# C++ compiler
CXX=g++

# C++ Compilation flags
CXX_FLAGS= \$(NEWLINE) $(TAB) -std=c++20 -Wall -O3

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
	INCLUDE_BOOST := $(addprefix -I,$(BOOST_DIR))
endif

# Include directory for the rebound library
INCLUDE_REBOUND := -I/usr/include/rebound

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

# rebound (integration of gravitational problem)
LD_LIB_REBOUND := -lrebound

# *************************************************************************************************
# Gather additional include directories with -I flag
# *************************************************************************************************

# Initialize INCLUDE to default usr/include directory
INCLUDE := -I/usr/include

# Include user software libraries if SOFTWARE_LIBRARY_DIR environment variable was set
ifdef SOFTWARE_LIBRARY_DIR
	INCLUDE := $(INCLUDE) $(INCLUDE_USR)
endif

# Only include boost if BOOST_DIR environment variable was set
# Note: on Ubuntu, if Boost was installed using b2, this is not necessary
# Loading a boost module will set this, e.g. $ module load boost/1.75
ifdef BOOST_DIR
	INCLUDE := $(INCLUDE) $(INCLUDE_BOOST)
endif

# Add rebound to the include directory
INCLUDE := $(INCLUDE) $(INCLUDE_REBOUND)

# *************************************************************************************************
# Generate linker flags
# -L: additional library directories
# -l: additional libraries in library search path named lib<library_name>.a
# *************************************************************************************************

# Additional library directories -L and other miscellanoues flags (e.g. address saniter)
LD_FLAGS := $(LD_FLAGS_USR) $(LD_ASAN)

# Additional libraries -l
LD_LIBS := $(LD_FMT) $(LD_LIB_BOOST) $(LD_LIB_GSL) $(LD_LIB_MARIADB) $(LD_LIB_REBOUND)

# *************************************************************************************************
# Address sanitizer configuration
# *************************************************************************************************

# Set flag whether to use address sanitizer; uncomment to turn this on
# USE_ASAN := 1
# The option passed to both compiler and linker when using address sanitizer
ASAN_OPT := -fsanitize=address

# Add this option if address sanitizer is being used
ifdef ASAN
	CXX_FLAGS := $(CXX_FLAGS) $(ASAN_OPT)
	LD_LIBS := $(LD_LIBS) $(ASAN_OPT)
endif
