BAMTOOLS_DIR := ./bamtools-2.3.0/
ABS_BAMTOOLS_DIR := $(realpath $(BAMTOOLS_DIR))
CXX := g++
CXXFLAGS := -Wno-deprecated -Wall -g -Wl,-rpath,$(ABS_BAMTOOLS_DIR)/lib/
INCLUDES := -I$(ABS_BAMTOOLS_DIR)/include/ -L$(ABS_BAMTOOLS_DIR)/lib/

all: bamtools Microassembler

.PHONY : Microassembler
Microassembler:
	cd Microassembler; make; cd ../

.PHONY : bamtools
bamtools:
	mkdir $(ABS_BAMTOOLS_DIR)/build; cd $(ABS_BAMTOOLS_DIR)/build; cmake ..; make; cd ../../

.PHONY : cleanbamtools
cleanbamtools:
	cd $(ABS_BAMTOOLS_DIR)/build; make clean; cd ../../

#.PHONY : clean
clean:
	rm -rf $(ABS_BAMTOOLS_DIR)/build; rm -f Microassembler/Microassembler
