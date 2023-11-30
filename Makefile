XX := g++
CXXFLAGS := -std=c++14 -flto -march=native -O3 -w -DNDEBUG
CPPFLAGS := -Iinclude
CPPTESTFLAGS := -TEST
EXTRAFLAGS := -DDYNAMIC_ORDERING -DFILTERING_BY_NEIGHBOR_SAFETY -DLEAF_ADAPTIVE_MATCHING -DN_OPTIMIZATION -DPRUNING_BY_EQUIVALENCE_SETS
TEST_LIBS := -lgtest -lpthread

SRC := src
OBJ := obj
TEST_SRC := tests

SRCS=$(wildcard $(SRC)/*.cpp)
OBJS=$(SRCS:$(SRC)/%.cpp=$(SRC)/%.o)

DEFAULT_TEST_SRCS=$(addprefix $(TEST_SRC)/, test_main.cpp test_compare.cpp)
CI_TEST_SRCS=$(addprefix $(TEST_SRC)/, test_main.cpp test_readGraph.cpp)
ALL_TEST_OBJS=$(wildcard $(TEST_SRC)/*.o)

ifeq ($(TARGET), CI)
TEST_SRCS=$(CI_TEST_SRCS)
else
TEST_SRCS=$(DEFAULT_TEST_SRCS)
endif
TEST_OBJS=$(TEST_SRCS:$(TEST_SRC)/%.cpp=$(TEST_SRC)/%.o)

VEQ_S := VEQ_S
VEQ_M := VEQ_M
TESTP := test_run

.PHONY: all clean run

all: $(VEQ_S) $(VEQ_M)

$(VEQ_S): $(SRC)/main.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(EXTRAFLAGS) -c $(SRC)/main.cpp -o $(SRC)/main.o
	$(CXX) $(CXXFLAGS) $(SRC)/main.o -o $@

$(VEQ_M): $(SRC)/main.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(EXTRAFLAGS) -DSUBGRAPH_MATCHING -c $(SRC)/main.cpp -o $(SRC)/main_matching.o
	$(CXX) $(CXXFLAGS) $(SRC)/main_matching.o -o $@

$(SRC)/%.o: $(SRC)/%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(TESTP): $(LIBGI) $(TEST_OBJS)
	$(CXX) -o $(TESTP) $(TEST_OBJS) $(LIBGI) $(TEST_LIBS)

$(TEST_SRC)/%.o: $(TEST_SRC)/%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(CPPTESTFLAGS) $(EXTRAFLAGS) -c $< -o $@

clean:
	$(RM) -rv $(VEQ_S)
	$(RM) -rv $(VEQ_M)
	$(RM) -rv $(OBJS)
	$(RM) -rv $(ALL_TEST_OBJS) $(TESTP)

runs: $(VEQ_S)
	./$(VEQ_S) -dg graph/search/data/COLLAB.gfu -qg graph/search/query/COLLAB/randomwalk/8/q30.gfu

runm: $(VEQ_M)
	./$(VEQ_M) -dg graph/matching/data/yeast.gfu -qg graph/matching/query/yeast/sparse/50/q30.gfu

test: $(TESTP)
	./$(TESTP)
