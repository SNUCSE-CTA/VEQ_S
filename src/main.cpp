// Adaptive matching order with new extendable candidate definition
#include <gtest/gtest.h>

#include "DAF.h"

void writeAnswer(string& answerFile) {
  ofstream outFile(answerFile);
  for (int g = 0; g < nGraph; ++g) {
    if (answer[g]) outFile << g << endl;
  }
  outFile.close();
}

inline bool BuildCS() {
  if (!FilterByCount(*currQ, *currG)) return false;
  root = SelectRootVertex(*currQ, *currG);
  ExtractResidualStructures(*currQ);
  BuildDAG(*currQ);
  // 1st refinement without neighbor-safety
  if (!TopDownInitial(*currQ, *currG)) return false;
  // Neighbor-safety in the 2nd refinement for subgraph search
  if (!filteringWithDAG(*currQ, *currG, bottomUp, childNgb, false))
    return false;
  // Neighbor-safety in the 3rd refinement
  if (!filteringWithDAG(*currQ, *currG, topDown, parentNgb, true)) return false;
  if (!ConstructAdjacencyList(*currQ, *currG, true)) return false;
  return true;
}

inline void ProcessQuery() {
  bool isCandidate;
  chrono::high_resolution_clock::time_point startClock, endClock;
  for (int x = 0; x < queryGraph.size(); ++x) {
    exceedTimeLimit = false;
    filteringTime = 0;
    verificationTime = 0;
    nMatch = 0;
    nCandidate = 0;
    recursiveCallCount = 0;
    currQ = queryGraph[x];
    AllocateForQueryGraph(*currQ);
    GetQueryDataStructure(*currQ);
    isTree = true;
    for (int i = 0; i < nQueryVertex; i++) {
      if (currQ->core[i] >= 2) {
        isTree = false;
        break;
      }
    }
    for (int y = 0; y < dataGraph.size(); ++y) {
      SetQueryGraphResource(*currQ);
      currG = dataGraph[y];
      // CS filtering
      startClock = GetClock();
      isCandidate = BuildCS();
      endClock = GetClock();
      filteringTime += TimeDiffInMilliseconds(startClock, endClock);
      if (!isCandidate) continue;
      ++nCandidate;
      // Backtracking
      startClock = GetClock();
      Backtrack(*currQ, *currG, y);
      endClock = GetClock();
      verificationTime += TimeDiffInMilliseconds(startClock, endClock);
    }
    totalTime = filteringTime + verificationTime;
    cout << "Total Recursive Call Count: " << recursiveCallCount << endl;
    cout << "Number of Candidate Graphs: " << nCandidate << endl;
    cout << "Number of Answer Graphs: " << nMatch << endl;
    writeAnswer(answerFile);
    cout << "Filtering Time (ms): " << filteringTime << endl;
    cout << "Verification Time (ms): " << verificationTime << endl;
    cout << "Processing Time (ms): " << totalTime << endl;
  }
}

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

// Subgraph query
// Usage: ./Program -d[dataFormat] [dataFile] -q[queryFormat] [queryFile] -o
// [answerFile]

// Subgraph matching
// Usage: ./Program -d[dataFormat] [dataFile] -q[queryFormat] [queryFile]

// dataFormat:
// -i: Igraph format
// -g: GFU format

// queryFormat:
// -c: CFL format
// -i: Igraph format
// -g: GFU format

int main(int argc, char* argv[]) {
  setvbuf(stdout, NULL, _IOLBF, 0);
  int queryFormat = GFU_FORMAT;
  int dataFormat = GFU_FORMAT;
  for (int i = 1; i < argc; ++i) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
        case 'd':
          switch (argv[i][2]) {
            case 'i':
              dataFormat = IGRAPH_FORMAT;
              break;
            case 'g':
              dataFormat = GFU_FORMAT;
              break;
          }
          dataGraphFile = argv[i + 1];
          break;
        case 'q':
          switch (argv[i][2]) {
            case 'c':
              queryFormat = CFL_FORMAT;
              break;
            case 'i':
              queryFormat = IGRAPH_FORMAT;
              break;
            case 'g':
              queryFormat = GFU_FORMAT;
              break;
          }
          queryGraphFile = argv[i + 1];
          break;
        case 'o':
          answerFile = argv[i + 1];
          break;
        case 't':
          ::testing::InitGoogleTest(&argc, argv);
          return RUN_ALL_TESTS();
      }
    }
  }
  cout << "Data file: " << dataGraphFile << endl;
  cout << "Query file: " << queryGraphFile << endl;
  cout << "Output file: " << answerFile << endl;
  if (dataFormat == IGRAPH_FORMAT)
    ReadIgraphFormat(dataGraphFile, dataGraph);
  else
    ReadGFUFormat(dataGraphFile, dataGraph);

  cout << "ReadGFUFormat" << endl;

  ProcessDataGraphs();

  cout << "ProcessDataGraphs" << endl;

  AllocateForDataGraph();

  cout << "AllocateForDataGraph" << endl;

  if (queryFormat == CFL_FORMAT)
    ReadCFLFormat(queryGraphFile, queryGraph);
  else if (queryFormat == IGRAPH_FORMAT)
    ReadIgraphFormat(queryGraphFile, queryGraph);
  else
    ReadGFUFormat(queryGraphFile, queryGraph);

  ProcessQuery();

  cout << "ProcessQuery" << endl;

  return 0;
}
