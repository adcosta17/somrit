# Sub directories containing source code, except for the main programs
SUBDIRS :=. 

LIBS = -lz -lm -llzma -lbz2 -lpthread -lcurl -lssl -lcrypto -lgcov -ldeflate -labpoa -lminimap2 -lspoa -lz -lhts

HTS_STATIC_LIB=./htslib/libhts.so
HTS_LIB=./htslib/libhts.a
HTS_INCLUDE=-I./htslib -I./htslib/htslib -L./htslib/htslib -L./htslib

PARASAIL_INCLUDE= -I./parasail/include -L./parasail/lib

MINIMAP_INCLUDE= -I./minimap2 -L./minimap2

ABPOA_INCLUDE= -I./abPOA/include -I./abPOA -I./abPOA/src -L./abPOA/lib

WFA_LIB=./WFA/build/libwfa.so
WFA_INCLUDE=-I./WFA -I./WFA/bin -I./WFA/build/ -L./WFA/build/ -L./WFA/bin/ -L./WFA

EDLIB_INCLUDE=-I./edlib/edlib/include -L./edlib

GABA_INCLUDE=-I./libgaba -L./libgaba

SPOA_INCLUDE=-Ispoa/include/ -Ispoa -Lspoa/build/lib/  

# Build libhts
#
htslib/libhts.so:
	cd htslib && make || exit 255

abPOA/lib/libabpoa.a:
	cd abPOA && make

gaba_wrapper.o: gaba_wrapper.c
	gcc -O3 -Wall -Wno-unused-function -Wno-unused-label -std=c99 -pipe -c gaba_wrapper.c $(GABA_INCLUDE) libgaba/*.o -o gaba_wrapper.o

# Source files
# g++ -c -Wall -Wextra -O2 -std=c++14 -fPIC ReAlign.cpp $(HTS_INCLUDE) $(PARASAIL_INCLUDE) $(WFA_INCLUDE) -lhts -lparasail -lwfa -o ReAlign.o
ReAlign: ReAlign.cpp ReAlign.hpp
	#g++ -Wall -Wextra -O2 -std=c++14 -fPIC -c ./edlib/edlib/src/edlib.cpp $(EDLIB_INCLUDE) -o edlib.o
	g++ -c -g -Wall -Wextra -O2 -std=c++14 -fPIC ReAlign.cpp $(HTS_INCLUDE) $(MINIMAP_INCLUDE) $(SPOA_INCLUDE) $(ABPOA_INCLUDE) -lhts -labpoa -lminimap2 -lspoa -o ReAlign.o

# g++ -shared -Wall -Wextra -std=c++14 -fPIC ReAlign.o $(HTS_LIB) $(HTS_INCLUDE) $(PARASAIL_INCLUDE) $(WFA_INCLUDE) $(LIBS) $(LDFLAGS) -o librealign.so
librealign.so: ReAlign.o
	g++ -shared -g -Wall -Wextra -std=c++14 -fPIC -o librealign.so ReAlign.o $(HTS_STATIC_LIB) $(HTS_INCLUDE) $(MINIMAP_INCLUDE) $(SPOA_INCLUDE) $(ABPOA_INCLUDE) $(LIBS) $(LDFLAGS)
	
test: ReAlign.cpp ReAlign.hpp
	g++ -g -Wall -Wextra -O2 -std=c++14 -fPIC ReAlign.cpp $(HTS_LIB) $(HTS_INCLUDE) $(MINIMAP_INCLUDE) $(SPOA_INCLUDE) $(ABPOA_INCLUDE) $(LIBS) $(LDFLAGS) -o test
	#g++ -g -Wall -Wextra -O2 -std=c++14 -fPIC -o test ReAlign.o $(HTS_INCLUDE) $(MINIMAP_INCLUDE) $(LIBS) $(LDFLAGS)


	