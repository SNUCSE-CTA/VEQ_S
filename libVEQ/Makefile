CXX = g++
CXXFLAGS = -std=c++14 -march=native -flto -O3 -w -DNDEBUG -shared -fPIC
CPPFLAGS := -I../include
LDFLAGS = -Wl,-undefined,dynamic_lookup

SRC = libVEQ.cpp
TARGET = libVEQ.so

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm -f $(TARGET)

.PHONY: all clean
