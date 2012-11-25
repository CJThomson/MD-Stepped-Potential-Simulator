
#This is the makefile for the Discrete Potential Simulator
COMP=g++
FILE=ESPSim.cpp Simulator.cpp Scheduler.cpp Engine.cpp Sampler.cpp NeighbourList.cpp
FILEWPUGI=ESPSim.cpp Simulator.cpp Scheduler.cpp Engine.cpp Sampler.cpp NeighbourList.cpp pugixml/pugixml.cpp
OBJECTS=$(FILE:.cpp=.o)
PUGI = pugixml/pugixml.o
FLAGS=-O3
OUT=ESPSim
LINK=-lboost_program_options

full: $(FILE)
	$(COMP) $(FILE) $(PUGI) -o $(OUT) $(LINK) $(FLAGS)
static: $(FILE)
	$(COMP) $(FLAGS) $(FILE) $(PUGI) -o $(OUT) /usr/lib64/libboost_program_options-mt-a
sim: $(FILE) $(OUT)
$(OUT): $(OBJECTS)
	$(COMP) $(OBJECTS) $(PUGI) $(FLAGS) -o $(OUT) $(LINK)

wpugi: $(FILEWPUGI)
	$(COMP) $(FILEWPUGI) -o $(OUT) $(LINK) $(FLAGS)
debug: $(FILE)
	$(COMP) $(FILE) $(PUGI) -o $(OUT) $(LINK) -g
clean:
	rm -f *.o core
rebuild: clean build

.cpp.o:
	$(COMP) -c $(FLAGS) $< -o $@