SRCS=$(wildcard *.cpp)
BINS=$(SRCS:.cpp=)

CXXFLAGS := $(shell pkg-config --cflags ibex-affine)
LIBS	 := $(shell pkg-config --libs  ibex-affine)

ifeq ($(DEBUG), yes)
CXXFLAGS := $(CXXFLAGS) -O0 -g -pg -Wall
else
CXXFLAGS := $(CXXFLAGS) -O3 -DNDEBUG
endif

all: $(BINS)

% :	%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(LIBS)

clean:
	rm -f $(BINS)
	
