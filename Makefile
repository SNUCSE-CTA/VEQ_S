CXX := g++
CXXFLAGS := -std=c++14 -O3 -w -DNDEBUG
CPPFLAGS := -Iinclude
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
TESTP := test_run

.PHONY: all clean run

all: $(VEQ_S)

$(VEQ_S): $(SRC)/main.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $(SRC)/main.cpp -o $(SRC)/main.o
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
	$(RM) -rv $(ALL_TEST_OBJS) $(TESTP)

run: $(VEQ_S)
	./$(VEQ_S) -dg graph/data/COLLAB.gfu -qg graph/query/COLLAB/randomwalk/8/q30.gfu

test: $(TESTP)
	./$(TESTP)
