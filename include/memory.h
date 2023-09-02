#include <vector>

#include "structure.h"

using namespace std;

unordered_map<Array, int, ArrayHash> cellToID;
int* cellPos;
int* cellList;
unordered_map<Array, int, ArrayHash>::iterator uIter;
unordered_map<int, Array>::iterator iIter;

Stack element[MAX_NUM_VERTEX];
Stack* currE;

Queue queueFllExt;

// Variables for maximal matching
int nodeIndex = 0;
int dist[MAX_NUM_VERTEX + 1];
int* uCand;
int* necMapping;
int* pairU;
int* pairV;
number_type leafCandSize = 0;
int leafCandIndex = 0;
vector<int> leafCand;
pair<int, int>* leafCandInfo;
int** sumNECCand;  // added for the optimization of performing all MM first
int* nSumNECCand;  // added for the optimization of performing all MM first
int* vCandIndex;   // length: maxNumDataVertex. its each entry stores the number
                   // of query nodes contain the corresponding candidate.
int posSumNECCand;
vector<pair<int, int>>* vCandInfo;  // a combination of v_cands and v_cands_pos
pair<int, int>* uCandInfo;          // two-d array:
// [i,j] means the j-th candidate pair of the i-th nec node(in the nec-node
// categorized by label) each candidate pair stores a candidate and the position
// of "the i-th nec node" in this candidate's u_cand set.
int* flagSumNECCand;  // indicate whether or not a candidate has been added into
                      // the candidate already.

int* arrayToClean;
int toCleanIndex = 0;
char* isVisited;  // data count
int* candToPos;   // data count

int uSequence[MAX_NUM_VERTEX];
int uSequenceSize = 0;
int orderDAG[MAX_NUM_VERTEX];

// Main result variables
int nCandidate = 0;
double nMatch, nCurrMatch, nRemainingMatch;
long long int recursiveCallCount = 0;
long long int recursiveCallCountPerGraph;
bool* answer;

inline void AllocateForDataGraph() {
#ifdef HUGE_GRAPH
  DATAV = maxNumDataVertex + 10;
  LABEL = nUniqueLabel + 10;
  color = new int*[LABEL];
  for (int i = 0; i < LABEL; i++) {
    color[i] = new int[2];
    for (int j = 0; j < 2; j++) {
      color[i][j] = 0;
    }
  }
  set_color = new int[LABEL];
  for (int i = 0; i < LABEL; i++) {
    set_color[i] = 0;
  }
  vis_adj = new int[DATAV];
  for (int i = 0; i < DATAV; i++) {
    vis_adj[i] = 0;
  }
  num_ngb_existence = new int*[DATAV];
  for (int i = 0; i < DATAV; i++) {
    num_ngb_existence[i] = new int[2];
    for (int j = 0; j < 2; j++) {
      num_ngb_existence[i][j] = 0;
    }
  }
  cnt_included_cs = new lint[DATAV];
  for (int i = 0; i < DATAV; i++) {
    cnt_included_cs[i] = (lint)(0);
  }
#endif
  pairV = new int[maxNumDataVertex];

  vCandIndex = new int[maxNumDataVertex];
  memset(vCandIndex, 0, sizeof(int) * maxNumDataVertex);

  flagSumNECCand = new int[maxNumDataVertex];
  memset(flagSumNECCand, -1, sizeof(int) * maxNumDataVertex);

  vCandInfo = new vector<pair<int, int>>[maxNumDataVertex];

  sumNECCand = new int*[nUniqueLabel + 1];
  for (int i = 0; i < (nUniqueLabel + 1); i++)
    sumNECCand[i] = new int[maxNumCandidate];
  nSumNECCand = new int[nUniqueLabel + 1]();

  mapTo = new int[maxNumDataVertex];
  memset(mapTo, -1, sizeof(int) * maxNumDataVertex);

  arrayToClean = new int[maxNumDataVertex];
  isVisited = new char[maxNumDataVertex];
  memset(isVisited, 0, sizeof(char) * maxNumDataVertex);
  candToPos = new int[maxNumDataVertex];
  memset(candToPos, -1, sizeof(int) * maxNumDataVertex);

  answer = new bool[nGraph];
  memset(answer, false, sizeof(bool) * nGraph);
  cellPos = new int[maxNumCandidate];
  cellList = new int[maxNumDataVertex];
}

inline void AllocateForQueryGraph(Graph& query) {
  nQueryVertex = query.nVertex;
  query.NECMapping = new int[nUniqueLabel * nQueryVertex];
  query.NECElement = new Element[nUniqueLabel * nQueryVertex];
  query.NECMap = new int[nQueryVertex];
  query.isProblemLeaf = new bool[nQueryVertex];
  FAILING_SET_SIZE = (nQueryVertex + 63) >> 6;

  leafCandSize = nQueryVertex * maxNumDataVertex;
  leafCand.resize(leafCandSize);  // tight upper bound
  leafCandInfo = new pair<int, int>[(nUniqueLabel + 1) * nQueryVertex];
  necMapping = new int[nQueryVertex + 1];
  uCand = new int[nQueryVertex * maxNumCandidate];  // definite upper bound

  pairU = new int[nQueryVertex + 1];

  uCandInfo = new pair<int, int>[nQueryVertex *
                                 maxNumCandidate];  // definite upper bound

  memset(localFlagQuery, 0, sizeof(char) * nQueryVertex);

  for (int i = 0; i < maxNumDataVertex; i++) vCandInfo[i].resize(nQueryVertex);

  nMaxCell = maxNumCandidate;
  CELL_SIZE = ((maxNumCandidate + 63) >> 6);
  for (int i = 0; i < nQueryVertex; i++) {
    element[i].failingSet = new uint64_t[FAILING_SET_SIZE];
    element[i].conflictCell = new uint64_t[CELL_SIZE]();
    element[i].isPruned = new bool[maxNumCandidate]();
    element[i].nPruned = 0;
    element[i].pruned = new int[maxNumCandidate];
    element[i].problemChild = new int[query.maxDegree];
    element[i].nProblemChild = 0;
    for (int j = 0; j < MAX_QUERY_DEGREE; j++) iec[i][j] = new int[maxDegree];
    if (useFailingSet) ancestors[i] = new uint64_t[FAILING_SET_SIZE];
  }

  // For BuildDAG
  for (int i = 0; i < nQueryVertex; i++) {
    DAG_parent_query[i] = new int[query.maxDegree];
  }
  for (int i = 0; i < nQueryVertex; i++) {
    DAG_ngb_query[i] = new int[query.maxDegree];
    DAG_ngb_query_ngb_index[i] = new int[query.maxDegree];
  }

  // For ConstructAdjacencyList
  for (int u = 0; u < nQueryVertex; ++u) {
    CandidateSpace& currSet = candSpace[u];
    currSet.candidates = new int[maxNumCandidate];
    currSet.weight = new long long[maxNumCandidate];
    currSet.nAdjacent = new int*[query.maxDegree];
    currSet.adjacent = new int**[query.maxDegree];
    currSet.capacity = new int*[query.maxDegree];
    for (int k = 0; k < query.maxDegree; ++k) {
      currSet.nAdjacent[k] = new int[maxNumCandidate];
      currSet.adjacent[k] = new int*[maxNumCandidate];
      currSet.capacity[k] = new int[maxNumCandidate];
      memset(currSet.capacity[k], 0, sizeof(int) * maxNumCandidate);
    }
    currSet.cell = new int[maxNumCandidate];
    currSet.cellVertex = new int[maxNumCandidate];
  }
  posArray = new int[maxNumCandidate];
  candArray = new int[maxNumCandidate];
  candOffset = new int[maxNumCandidate];
}

inline void Deallocate(Graph& query) {
  delete[] query.NECMapping;
  delete[] query.NECElement;
  delete[] query.NECMap;
  delete[] query.isProblemLeaf;
}

inline void SetQueryGraphResource(Graph& query) {
  globalCellID = 1;
  cellToID.clear();
  cellPos[0] = 0;
  memset(visitedForQuery, 0, sizeof(char) * nQueryVertex);
  memset(cntLabel, 0, sizeof(int) * nUniqueLabel);
  for (int u = 0; u < nQueryVertex; ++u) {
    candSpace[u].size = 0;
    candSpace[u].nCellVertex = 0;
    memset(candSpace[u].cell, -1, sizeof(int) * maxNumCandidate);
  }
}

// for test
inline void initialize(vector<Graph*>& dataGraph, vector<Graph*>& queryGraph) {
  // memory.h
  nodeIndex = 0;
  leafCandSize = 0;
  leafCandIndex = 0;
  toCleanIndex = 0;
  nCandidate = 0;
  recursiveCallCount = 0;

  // process.h
  globalCellID = 1;
  maxNumCandidate = 0;
  maxDegree = 0;
  maxNumDataVertex = 0;
  nQueryVertex = 0;

  // util.h
  labelID = -1;
  FAILING_SET_SIZE = ((MAX_NUM_VERTEX + 63) >> 6);
  // Variables for filtering by neighbor safety
  on = false;
  index_vis_color = 0;
  index_set_color = 0;
  index_vis_adj = 0;
  index_ngb_existence = 0;
  s1 = 0;
  s2 = 0;

  // data graph and query graph
  leafCand.clear();
  delete[] leafCandInfo;
  delete[] necMapping;
  delete[] uCand;
  delete[] pairU;
  delete[] uCandInfo;
  for (int i = 0; i < maxNumDataVertex; ++i) vCandInfo[i].clear();

  for (int i = 0; i < nQueryVertex; i++) {
    delete[] element[i].failingSet;
    delete[] element[i].conflictCell;
    delete[] element[i].isPruned;
    delete[] element[i].pruned;
    delete[] element[i].problemChild;
    for (int j = 0; j < MAX_QUERY_DEGREE; j++) delete[] iec[i][j];
    if (useFailingSet) delete[] ancestors[i];
  }

  // For BuildDAG
  for (int i = 0; i < nQueryVertex; i++) {
    delete[] DAG_parent_query[i];
  }
  for (int i = 0; i < nQueryVertex; i++) {
    delete[] DAG_ngb_query_ngb_index[i];
    delete[] DAG_ngb_query[i];
  }

  // For ConstructAdjacencyList
  for (int u = 0; u < nQueryVertex; ++u) {
    CandidateSpace& currSet = candSpace[u];
    delete[] currSet.candidates;
    delete[] currSet.weight;
    for (int k = 0; k < queryGraph[0]->maxDegree; ++k) {
      delete[] currSet.nAdjacent[k];
      delete[] currSet.adjacent[k];
      delete[] currSet.capacity[k];
    }
    delete[] currSet.capacity;
    delete[] currSet.nAdjacent;
    delete[] currSet.adjacent;

    delete[] currSet.cell;
    delete[] currSet.cellVertex;
  }
  delete[] posArray;
  delete[] candArray;
  delete[] candOffset;

  for (Graph* graphPtr : dataGraph) Deallocate(*graphPtr);
  for (Graph* graphPtr : queryGraph) Deallocate(*graphPtr);
  dataGraph.clear();
  queryGraph.clear();
}
