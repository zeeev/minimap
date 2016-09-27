######################################
# Makefile written by Zev Kronenberg #
#     zev.kronenberg@gmail.com       #
######################################

CXX=g++
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)


.PHONY: all clean

all: bin/chain bin/layout

bin/chain: createBin
	cd chain && $(CXX) *cpp -o ../bin/qChain
bin/layout: createBin
	cd layout && $(CXX) *cpp -o ../bin/layout

createBin:
	-mkdir bin

clean:
	 -cd chain && rm *.o && cd ../layout && rm *.o && cd ../bin && rm -rf *