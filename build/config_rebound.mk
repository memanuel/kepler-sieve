# MSE configuration
CC := gcc
OPT+= -Wpointer-arith -D_GNU_SOURCE -O3 -march=native
FFTW := 1
OPENGL := 1
OPENMP := 0
# Check operating system
ifndef OS
	OS=$(shell uname)
endif

# Additional flags for Linux
ifeq ($(OS), Linux)
	OPT+= -Wall -g -Wno-unused-result
	LIB+= -lm -lrt
endif

# Check for FFTW
ifeq ($(FFTW), 1)
	PREDEF+= -DFFTW
	LIB+= -lfftw3
endif

# Check for OPENGL
ifeq ($(OPENGL), 1)
	PREDEF+= -DOPENGL
	LIB+= -lglfw
endif

# Check for OPENMp
ifeq ($(OPENMP), 1)
	PREDEF+= -DOPENMP
	OPT+= -fopenmp
	LIB+= -fopenmp
endif

# Check for GITHASH
ifndef GITHASH
	GITHASH = $(shell git rev-parse HEAD || echo '0000000000gitnotfound0000000000000000000')
	PREDEF+= -DGITHASH=$(GITHASH)
endif
