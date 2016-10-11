######################################
# Makefile written by Zev Kronenberg #
#     zev.kronenberg@gmail.com       #
######################################

CXX=g++
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)
CFLAGS       = -Wall -DVERSION=\"$(GIT_VERSION)\" -std=c++0x -g
INCLUDE      = -I lib/
CPP_FILES    = $(wildcard lib/*cpp)
OBJ_FILES    = $(addprefix lib/,$(notdir $(CPP_FILES:.cpp=.o)))
BIN_FILES    = $(wildcard src/*cpp)
BIN_NAMES    = $(notdir $(BIN_FILES:.cpp=))

.PHONY: bins

bins: createBin bin/layout bin/qChain bin/synblocks

bin/synblocks: $(OBJ_FILES)
	$(CXX) $(CFLAGS) $(INCLUDE) $(OBJ_FILES) -o bin/synblocks src/synblocks.cpp

bin/layout: $(OBJ_FILES)
	$(CXX) $(CFLAGS) $(INCLUDE) $(OBJ_FILES) -o bin/layout src/layout.cpp

bin/qChain: $(OBJ_FILES)
	$(CXX) $(CFLAGS) $(INCLUDE) $(OBJ_FILES) -o bin/qChain src/chaining.cpp

src/%.o: src/%.cpp $(BAMTOOLS_LIB)
	$(CXX) $(CFLAGS) $(INCLUDE) -c -o $@ $<

createBin:
	-mkdir bin

clean:
	-cd lib && rm *.o && rm -rf ../bin