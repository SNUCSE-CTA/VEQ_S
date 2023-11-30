#include <gtest/gtest.h>

#include "run.h"

TEST(graph_test, graph_test_readGraph_1) {
  string queryGraphFiles[1] = {"graph/search/query/IMDB-MULTI/bfs/8/q0.gfu"};
  ReadGFUFormat("graph/search/data/IMDB-MULTI.gfu", dataGraph);
  ASSERT_EQ(false, dataGraph.back()->fail);
}
