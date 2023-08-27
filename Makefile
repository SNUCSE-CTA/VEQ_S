CXX := g++
CXXFLAGS := -std=c++14 -O3 -w -DNDEBUG
CPPFLAGS := -Iinclude
EXTRA_FLAGS := -DDYNAMIC_ORDERING -DFILTERING_BY_NEIGHBOR_SAFETY -DLEAF_ADAPTIVE_MATCHING -DN_OPTIMIZATION -DPRUNING_BY_EQUIVALENCE_SETS
TEST_LIBS := -lgtest -lpthread

SRC := src
OBJ := obj
TEST_SRC := tests

SRCS=$(wildcard $(SRC)/*.cpp)
OBJS=$(SRCS:$(SRC)/%.cpp=$(SRC)/%.o)

TEST_SRCS=$(wildcard $(TEST_SRC)/*.cpp)
TEST_OBJS=$(TEST_SRCS:$(TEST_SRC)/%.cpp=$(TEST_SRC)/%.o)

VEQ_S := VEQ_S
TESTP := test_run

.PHONY: all clean run

all: $(VEQ_S)

$(VEQ_S): $(SRC)/main.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(EXTRA_FLAGS) -c $(SRC)/main.cpp -o $(SRC)/main.o
	$(CXX) $(CXXFLAGS) $(SRC)/main.o -o $@

$(SRC)/%.o: $(SRC)/%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(TESTP): $(LIBGI) $(TEST_OBJS)
	$(CXX) -o $(TESTP) $(TEST_OBJS) $(LIBGI) $(TEST_LIBS)

$(TEST_SRC)/%.o: $(TEST_SRC)/%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(CPPTESTFLAGS) -c $< -o $@

clean:
	$(RM) -rv $(VEQ_S)
	$(RM) -rv $(OBJS)
	$(RM) -rv $(TEST_OBJS) $(TESTP)

run: $(VEQ_S)
	./$(VEQ_S) -dg graph/data/COLLAB.gfu -qg graph/query/COLLAB/randomwalk/8/q30.gfu

test: $(TESTP)
	./$(TESTP)
