
#This is the makefile for the Discrete Potential Simulator
COMP=g++
FILE=ESPSim.cpp Simulator.cpp Scheduler.cpp Engine.cpp Sampler.cpp
FILEWPUGI=ESPSim.cpp Simulator.cpp Scheduler.cpp Engine.cpp Sampler.cpp pugixml/pugixml.cpp
OBJECTS=$(FILE:.cpp=.o)
PUGI = pugixml/pugixml.o
FLAGS=-O3
OUT=ESPSim
LINK=-l boost_program_options

sim: $(FILE) $(OUT)
$(OUT): $(OBJECTS)
	$(COMP) $(OBJECTS) $(PUGI) $(FLAGS) -o $(OUT) $(LINK)

full: $(FILE)
	$(COMP) $(FILE) $(PUGI) -o $(OUT) $(LINK) $(FLAGS)
wpugi: $(FILEWPUGI)
	$(COMP) $(FILEWPUGI) -o $(OUT) $(LINK) $(FLAGS)
debug: $(FILE)
	$(COMP) $(FILE) -o $(OUT) $(LINK) -g
clean:
	rm -f *.o core
rebuild: clean build

.cpp.o:
	$(COMP) -c $(FLAGS) $< -o $@