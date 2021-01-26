SOURCE = $(wildcard *.cpp)
OBJ = $(patsubst %.cpp, %.o, $(SOURCE))

CXXFLAGS = -Wall 

LDFLAGS = -lm 

ifeq (x$(debug), x)
	CXXFLAGS += -O3 -fopenmp
else
	CXXFLAGS += -g -DDEBUG
endif

ifneq (x$(cpp11), x)
	CXXFLAGS += -std=c++11
endif

ifneq (x$(cpp0x), x)
	CXXFLAGS += -std=c++0x
endif


all: main

main : $(OBJ) 
	$(CXX) -o $@ $(OBJ) $(CXXFLAGS) $(LDFLAGS)

%.o : %.cpp 
	$(CXX) -c -o $@ $(CPPFLAGS) $(CXXFLAGS) $< 

.PHONY: clean
clean:
	-rm $(OBJ)
	-rm main
