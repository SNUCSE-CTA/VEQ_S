CXX = g++
CXXFLAGS = -std=c++14 -march=native -flto -O3 -w -DNDEBUG -shared -fPIC
CPPFLAGS := -I../include
EXTRAFLAGS := -DDYNAMIC_ORDERING -DFILTERING_BY_NEIGHBOR_SAFETY -DLEAF_ADAPTIVE_MATCHING -DN_OPTIMIZATION -DPRUNING_BY_EQUIVALENCE_SETS
LDFLAGS = -Wl,-undefined,dynamic_lookup

SRC = libVEQ.cpp
TARGET = libVEQ.so

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(EXTRAFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm -f $(TARGET)

.PHONY: all clean
