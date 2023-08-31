#include "compare.h"

TEST(compare_test, compare_test_1) {
  string queryGraphFile = "graph/query/IMDB-MULTI/bfs/8/q0.gfu";
  string resultFile = "testdata/IMDB-MULTI/bfs/8/q0.txt";

  compare("graph/data/IMDB-MULTI.gfu", queryGraphFile, resultFile);
}

TEST(compare_test, compare_test_2) {
  string queryGraphFile = "graph/query/IMDB-MULTI/bfs/16/q0.gfu";
  string resultFile = "testdata/IMDB-MULTI/bfs/16/q0.txt";

  compare("graph/data/IMDB-MULTI.gfu", queryGraphFile, resultFile);
}
