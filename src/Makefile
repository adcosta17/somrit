## Makefile based on 10xtrim by Joanna Pineda, found here: https://github.com/jopineda/10xtrim/
# Sub directories containing source code, except for the main programs
SUBDIRS :=. 

LIBS = -lz -lm -llzma -lbz2  -lpthread -lcurl -lssl -lcrypto -lgcov -ldeflate -lparasail
CXXFLAGS ?= -g -O3 -Wall -Wextra
CXXFLAGS += -std=c++14 
CFLAGS ?= -O3 -std=c99
CXX ?= g++
CC ?= gcc

HTS_LIB=./htslib/libhts.a
HTS_INCLUDE=-I./htslib -I./htslib/htslib

PARASAIL_LIB=./parasail/libparasail.a
PARASAIL_INCLUDE= -I ./parasail/include -L ./parasail/lib

# Include the src subdirectories
NP_INCLUDE=$(addprefix -I./, $(SUBDIRS))

# Add include flags
CPPFLAGS += $(NP_INCLUDE) $(HTS_INCLUDE)

#
# Build libhts
#
htslib/libhts.a:
	cd htslib && make || exit 255


#
# Build libparasail
#

parasail/libparasail.a:
	cd parasail && autoreconf -fi && make || exit 255

#
# Source files
#


Align: ReAlign.cpp TsvRecord.cpp TsvRecord.hpp
	g++ -g -std=c++14 -c -o librealign.o ReAlign.cpp TsvRecord.cpp $(HTS_LIB) $(PARASAIL_LIB) $(LIBS) $(LDFLAGS)
	g++ -shared -o librealign.so librealign.o
	