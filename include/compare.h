#include <gtest/gtest.h>

#include "run.h"

void compare(string dataGraphFile, string queryGraphFile, string resultFile) {
  ReadGFUFormat(dataGraphFile, dataGraph);
  ProcessDataGraphs();
  AllocateForDataGraph();

  ReadGFUFormat(queryGraphFile, queryGraph);
  ProcessQuery();

  ifstream inFile(resultFile);
  if (!inFile) {
    cout << "Failed to open a data result file " << resultFile << endl;
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
    bool inExpectedResult =
        expectedResultIter != expectedResult.end() && *expectedResultIter == g;
    EXPECT_EQ(answer[g], inExpectedResult);
    if (inExpectedResult) expectedResultIter++;
  }

  initialize(dataGraph, queryGraph);
}
