# Sub directories containing source code, except for the main programs
SUBDIRS :=. 

LIBS = -lz -lm -llzma -lbz2 -lpthread -lcurl -lssl -lcrypto -lgcov -ldeflate -labpoa -lminimap2 -lz -lhts

HTS?=install

ifeq ($(HTS), install)
	HTS_INCLUDE=-I./htslib -I./htslib/htslib -L./htslib/htslib -L./htslib
else
	HTS_INCLUDE=
endif	

PARASAIL_INCLUDE= -I./parasail/include -L./parasail/lib

MINIMAP_INCLUDE= -I./minimap2 -L./minimap2

ABPOA_INCLUDE= -I./abPOA/include -I./abPOA -I./abPOA/src -L./abPOA/lib

# Build libhts
#

all: librealign.so

clean:
	rm librealign.so ReAlign.o

# If htslib needs to be compiled
htslib/libhts.so:
	cd htslib && make || exit 255

abPOA/lib/libabpoa.a:
	cd abPOA && make

minimap2/libminimap2.a:
	cd minimap2 && make CC='gcc' CFLAGS=' -g -Wall -O2 -Wc++-compat -fPIC'

ReAlign.o: ReAlign.cpp ReAlign.hpp abPOA/lib/libabpoa.a
	g++ -c -g -Wall -Wextra -O2 -std=c++14 -fPIC ReAlign.cpp $(HTS_INCLUDE) $(MINIMAP_INCLUDE) $(ABPOA_INCLUDE) -lhts -labpoa -lminimap2 -o ReAlign.o

librealign.so: ReAlign.o
	g++ -shared -g -Wall -Wextra -std=c++14 -fPIC -o librealign.so ReAlign.o $(HTS_INCLUDE) $(MINIMAP_INCLUDE) $(ABPOA_INCLUDE) $(LIBS) $(LDFLAGS)
	
