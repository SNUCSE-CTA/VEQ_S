#include "run.h"

bool compare(string dataGraphFile, string queryGraphFile, string resultFile) {
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

  if (recursiveCallCount != expectedResult[0]) return false;
  if (nCandidate != expectedResult[1]) return false;
  if ((long long)nMatch != expectedResult[2]) return false;
  vector<int>::iterator expectedResultIter = expectedResult.begin() + 3;
  for (int g = 0; g < nGraph; ++g) {
    bool inExpectedResult =
        expectedResultIter != expectedResult.end() && *expectedResultIter == g;
    if (answer[g] != inExpectedResult) return false;
    if (inExpectedResult) expectedResultIter++;
  }

  initialize(dataGraph, queryGraph);

  return true;
}
