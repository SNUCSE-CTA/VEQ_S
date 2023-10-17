// Adaptive matching order with new extendable candidate definition
#include "VEQ.h"

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
  buildNbrSafetyStructure(*currQ, *currG);
  if (!filteringWithDAG(*currQ, *currG, topDown, parentNgb, true)) {
    FreeNbrSafetyStructure(*currG);
    return false;
  }
  FreeNbrSafetyStructure(*currG);
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
