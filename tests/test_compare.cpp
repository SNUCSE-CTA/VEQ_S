#include <gtest/gtest.h>

#include "run.h"

void compare(string dataGraphFile, const int queryCnt, string queryGraphFiles[],
             string resultFiles[]) {
  ReadGFUFormat(dataGraphFile, dataGraph);

  ProcessDataGraphs();
  AllocateForDataGraph();

  for (int i = 0; i < queryCnt; i++) {
    ReadGFUFormat(queryGraphFiles[i], queryGraph);
    ProcessQuery();

    ifstream inFile(resultFiles[i]);
    if (!inFile) {
      cout << "Failed to open a data result file " << resultFiles[i] << endl;
      exit(1);
    }

    vector<int> expectedResult;
    while (!inFile.eof()) {
      int num;
      inFile >> num;
      expectedResult.push_back(num);
    }

    EXPECT_EQ(recursiveCallCount, expectedResult[0]);
    EXPECT_EQ(nCandidate, expectedResult[1]);
    EXPECT_EQ(nMatch, expectedResult[2]);
    vector<int>::iterator expectedResultIter = expectedResult.begin() + 3;
    for (int g = 0; g < nGraph; ++g) {
      bool inExpectedResult = expectedResultIter != expectedResult.end() &&
                              *expectedResultIter == g;
      EXPECT_EQ(answer[g], inExpectedResult);
      if (inExpectedResult) expectedResultIter++;
    }
  }
}

TEST(compare_test, compare_test_1) {
  string queryGraphFiles[1] = {"graph/query/IMDB-MULTI/bfs/8/q0.gfu"};
  string resultFiles[1] = {"testdata/IMDB-MULTI/bfs/8/q0.txt"};

  compare("graph/data/IMDB-MULTI.gfu", 1, queryGraphFiles, resultFiles);
}
