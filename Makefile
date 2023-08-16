CXX := g++
CXXFLAGS := -std=c++11 -O3 -w -DNDEBUG
CPPFLAGS := -Iinclude
EXTRA_FLAGS := -DDYNAMIC_ORDERING -DFILTERING_BY_NEIGHBOR_SAFETY -DLEAF_ADAPTIVE_MATCHING -DN_OPTIMIZATION -DPRUNING_BY_EQUIVALENCE_SETS
TEST_LIBS := -lgtest -lpthread

SRC := src
OBJ := obj

SRCS=$(wildcard $(SRC)/*.cpp)
OBJS=$(SRCS:$(SRC)/%.cpp=$(SRC)/%.o)

VEQ_S := VEQ_S

.PHONY: all clean run

all: $(VEQ_S)

$(VEQ_S): $(SRC)/main.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(EXTRA_FLAGS) -c $(SRC)/main.cpp -o $(SRC)/main.o
	$(CXX) $(CXXFLAGS) $(SRC)/main.o $(TEST_LIBS) -o $@

$(SRC)/%.o: $(SRC)/%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

clean:
	$(RM) -rv $(VEQ_S)
	$(RM) -rv $(OBJS)

run: $(VEQ_S)
	./$(VEQ_S) -dg graph/data/COLLAB.gfu -qg graph/query/COLLAB/randomwalk/8/q30.gfu

test: $(VEQ_S)
	./$(VEQ_S) -test
