
#This is the makefile for the Discrete Potential Simulator
COMP=g++
FILE=ESPSim.cpp Simulator.cpp Scheduler.cpp Engine.cpp pugixml/pugixml.cpp
FLAGS=-O3
OUT=ESPSim
LINK=/usr/lib/boost-1_46/libboost_program_options-mt.so

build: $(FILE)
	$(COMP) $(FILE) -o $(OUT) $(FLAGS) $(LINK)

debug: $(FILE)
	$(COMP) $(FILE) -o $(OUT) $(LINK) -g

clean:
	rm -f *.o core
rebuild: clean build
