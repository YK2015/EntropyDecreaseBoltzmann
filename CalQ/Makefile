HEADER = $(wildcard *.h)
BIN = $(patsubst %.cpp, %, $(wildcard *.cpp))

CXX = g++-8

CXXFLAGS = -w

LDFLAGS = -lm -lfftw3 -llapack

ifeq (x$(debug), x)
	CXXFLAGS += -O3 -std=c++14 
else
	CXXFLAGS += -g -DDEBUG
endif

ifneq (x$(cpp11), x)
	CXXFLAGS += -std=c++11
endif

ifneq (x$(cpp0x), x)
	CXXFLAGS += -std=c++0x
endif


all: $(BIN)
$(BIN): % : %.cpp $(HEADER)
	$(CXX) -o $@ $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS)

.PHONY: clean
clean:
	-rm $(BIN)
