######################################
# Makefile written by Zev Kronenberg #
#     zev.kronenberg@gmail.com       #
######################################

CXX=g++ -Wreorder
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)


.PHONY: all clean

all: bin/chain bin/layout

bin/chain:
	cd chain && $(CXX) *cpp -o ../bin/qChain
bin/layout:
	cd layout && $(CXX) *cpp -o ../bin/layout

clean:
	 -cd chain && rm *.o && cd ../layout && rm *.o && cd ../bin && rm -rf *