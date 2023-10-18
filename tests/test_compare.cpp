#include "compare.h"

TEST(compare_test, compare_test_1) {
  string queryGraphFile = "graph/query/IMDB-MULTI/bfs/8/q0.gfu";
  string resultFile = "testdata/IMDB-MULTI/bfs/8/q0.txt";

  EXPECT_TRUE(compare("graph/data/IMDB-MULTI.gfu", queryGraphFile, resultFile));
}

TEST(compare_test, compare_test_2) {
  string queryGraphFile = "graph/query/IMDB-MULTI/bfs/16/q0.gfu";
  string resultFile = "testdata/IMDB-MULTI/bfs/16/q0.txt";

  EXPECT_TRUE(compare("graph/data/IMDB-MULTI.gfu", queryGraphFile, resultFile));
}

TEST(compare_test, compare_test_3) {
  string queryGraphFile = "graph/query/IMDB-MULTI/bfs/32/q0.gfu";
  string resultFile = "testdata/IMDB-MULTI/bfs/32/q0.txt";

  EXPECT_TRUE(compare("graph/data/IMDB-MULTI.gfu", queryGraphFile, resultFile));
}

TEST(compare_test, compare_test_4) {
  string queryGraphFile = "graph/query/IMDB-MULTI/bfs/64/q0.gfu";
  string resultFile = "testdata/IMDB-MULTI/bfs/64/q0.txt";

  EXPECT_TRUE(compare("graph/data/IMDB-MULTI.gfu", queryGraphFile, resultFile));
}

TEST(compare_test, compare_test_5) {
  string queryGraphFile = "graph/query/IMDB-MULTI/randomwalk/8/q0.gfu";
  string resultFile = "testdata/IMDB-MULTI/randomwalk/8/q0.txt";

  EXPECT_TRUE(compare("graph/data/IMDB-MULTI.gfu", queryGraphFile, resultFile));
}

TEST(compare_test, compare_test_6) {
  string queryGraphFile = "graph/query/IMDB-MULTI/randomwalk/16/q0.gfu";
  string resultFile = "testdata/IMDB-MULTI/randomwalk/16/q0.txt";

  EXPECT_TRUE(compare("graph/data/IMDB-MULTI.gfu", queryGraphFile, resultFile));
}

TEST(compare_test, compare_test_7) {
  string queryGraphFile = "graph/query/IMDB-MULTI/randomwalk/32/q0.gfu";
  string resultFile = "testdata/IMDB-MULTI/randomwalk/32/q0.txt";

  EXPECT_TRUE(compare("graph/data/IMDB-MULTI.gfu", queryGraphFile, resultFile));
}

TEST(compare_test, compare_test_8) {
  string queryGraphFile = "graph/query/IMDB-MULTI/randomwalk/64/q0.gfu";
  string resultFile = "testdata/IMDB-MULTI/randomwalk/64/q0.txt";

  EXPECT_TRUE(compare("graph/data/IMDB-MULTI.gfu", queryGraphFile, resultFile));
}

TEST(compare_test, compare_test_9) {
  string queryGraphFile = "graph/query/pcms/bfs/8/q0.gfu";
  string resultFile = "testdata/pcms/bfs/8/q0.txt";

  EXPECT_TRUE(compare("graph/data/pcms.gfu", queryGraphFile, resultFile));
}

TEST(compare_test, compare_test_10) {
  string queryGraphFile = "graph/query/pcms/bfs/16/q0.gfu";
  string resultFile = "testdata/pcms/bfs/16/q0.txt";

  EXPECT_TRUE(compare("graph/data/pcms.gfu", queryGraphFile, resultFile));
}

TEST(compare_test, compare_test_11) {
  string queryGraphFile = "graph/query/pcms/bfs/32/q0.gfu";
  string resultFile = "testdata/pcms/bfs/32/q0.txt";

  EXPECT_TRUE(compare("graph/data/pcms.gfu", queryGraphFile, resultFile));
}

TEST(compare_test, compare_test_12) {
  string queryGraphFile = "graph/query/pcms/bfs/64/q0.gfu";
  string resultFile = "testdata/pcms/bfs/64/q0.txt";

  EXPECT_TRUE(compare("graph/data/pcms.gfu", queryGraphFile, resultFile));
}

TEST(compare_test, compare_test_13) {
  string queryGraphFile = "graph/query/pcms/randomwalk/8/q0.gfu";
  string resultFile = "testdata/pcms/randomwalk/8/q0.txt";

  EXPECT_TRUE(compare("graph/data/pcms.gfu", queryGraphFile, resultFile));
}

TEST(compare_test, compare_test_14) {
  string queryGraphFile = "graph/query/pcms/randomwalk/16/q0.gfu";
  string resultFile = "testdata/pcms/randomwalk/16/q0.txt";

  EXPECT_TRUE(compare("graph/data/pcms.gfu", queryGraphFile, resultFile));
}

TEST(compare_test, compare_test_15) {
  string queryGraphFile = "graph/query/pcms/randomwalk/32/q0.gfu";
  string resultFile = "testdata/pcms/randomwalk/32/q0.txt";

  EXPECT_TRUE(compare("graph/data/pcms.gfu", queryGraphFile, resultFile));
}

TEST(compare_test, compare_test_16) {
  string queryGraphFile = "graph/query/pcms/randomwalk/64/q0.gfu";
  string resultFile = "testdata/pcms/randomwalk/64/q0.txt";

  EXPECT_TRUE(compare("graph/data/pcms.gfu", queryGraphFile, resultFile));
}
