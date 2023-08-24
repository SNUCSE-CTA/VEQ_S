#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

#include "graph.h"
using namespace std;

CandidateSpace candSpace[MAX_NUM_VERTEX];
int uI;
int nCandOffset;
int globalCellID = 1;
int* posArray;
int* candArray;
int* candOffset;

// Basic variables for data graphs
Graph* currG;
int nGraph;
int maxNumCandidate = 0;
int maxDegree = 0;
int maxNumDataVertex = 0;
vector<Graph*> dataGraph;
unordered_map<string, int> labelMap;

// Basic variables for query graph
Graph* currQ;
int nQueryVertex = 0;
vector<Graph*> queryGraph;
// char * isMappedDataVertex;
int* mapTo;
bool isMapped[MAX_NUM_VERTEX];
int order[MAX_NUM_VERTEX];
int extendable[MAX_NUM_VERTEX];
// Core Decomposition of query
int arrForQuery[MAX_NUM_VERTEX];
int binQuery[MAX_NUM_VERTEX];
int posQuery[MAX_NUM_VERTEX];
int vertQuery[MAX_NUM_VERTEX];
char localFlagQuery[MAX_NUM_VERTEX];
char visitedForQuery[MAX_NUM_VERTEX];
int visitedForQueryInt[MAX_NUM_VERTEX];
vector<pair<int, int>>
    NECElementOffset;  // include a redunant element to set the end
unordered_map<pair<int, int>, int, pairHash> NECChild;
unordered_map<pair<int, int>, int, pairHash> NECSize;

// Basic variables for query DAG
int root;
#ifdef NEIGHBOR_LABEL_FREQUENCY
pair<int, int>* NLFArray[MAX_NUM_VERTEX];
int numNLF[MAX_NUM_VERTEX];
#endif
int* DAG_parent_query[MAX_NUM_VERTEX];  // DAG_parent_query[ui] contains parents
                                        // of ui in query DAG.
int* DAG_child_query[MAX_NUM_VERTEX];   // DAG_child_query[ui] contains children
                                        // of ui in query DAG.
int* DAG_child_query_parent_index
    [MAX_NUM_VERTEX];  // DAG_child_query_parent_index[ui][j] contains parent
                       // index of a query vertex ui in DAG_parent_query[cj]
                       // where cj is the j-th child node of ui.
int DAG_parent_query_size[MAX_NUM_VERTEX];  // DAG_parent_query_size[ui]
                                            // contains the number of parents of
                                            // ui in query DAG.
int DAG_child_query_size[MAX_NUM_VERTEX];   // DAG_child_query_size[ui] contains
                                            // the number of children of ui in
                                            // query DAG.

int* DAG_ngb_query[MAX_NUM_VERTEX];
int* DAG_ngb_query_ngb_index[MAX_NUM_VERTEX];
int DAG_ngb_query_size[MAX_NUM_VERTEX];

// bit_vector DAG_ancestor[MAX_NUM_VERTEX];
uint64_t* DAG_ancestor[MAX_NUM_VERTEX];  // 20190708

chrono::high_resolution_clock::time_point GetClock() {
  return chrono::high_resolution_clock::now();
}
double TimeDiffInMilliseconds(
    const chrono::high_resolution_clock::time_point start,
    const chrono::high_resolution_clock::time_point end) {
  chrono::duration<double> durationDiff = end - start;
  return durationDiff.count() * 1000;
}

inline void ReadCFLFormat(const string& graphFile, vector<Graph*>& graphList) {
  ifstream inFile(graphFile);
  if (!inFile) {
    cout << "Failed to open a data graph file " << graphFile << endl;
    exit(1);
  }
  char c;
  int id, m, currLabel, currDegree, labelID, matrixIndex;
  number_type numVertex, sumDegree;
  string line;
  unordered_map<string, int>::const_iterator sIter;
  unordered_map<int, int>::const_iterator iIter;
  Graph* graph = NULL;
  while (!inFile.eof()) {
    inFile >> c >> id >> numVertex >> sumDegree;
    if (inFile.eof()) break;
    graph = new Graph();
    graphList.push_back(graph);
    graph->initVertex(numVertex);
    graph->initEdge(sumDegree / 2);
    matrixIndex = 0;
    for (int i = 0; i < numVertex; ++i) {
      graph->nbrOffset[i] = matrixIndex;
      inFile >> m >> currLabel >> currDegree;
      if (currDegree > graph->maxDegree) graph->maxDegree = currDegree;
      sIter = labelMap.find(to_string(currLabel));
      if (sIter != labelMap.end()) {
        labelID = sIter->second;
      } else {
        labelID = labelMap.size();
        labelMap[to_string(currLabel)] = labelID;
      }
      iIter = graph->labelFrequency.find(labelID);
      if (iIter != graph->labelFrequency.end())
        graph->labelFrequency[labelID] += 1;
      else
        graph->labelFrequency[labelID] = 1;
      graph->label[i] = labelID;
      graph->core[i] = graph->degree[i] = currDegree;
      for (int j = 0; j < currDegree; ++j) {
#ifdef EDGE_LABEL
        inFile >> m;
        graph->nbr.set(matrixIndex, m);
#else
        inFile >> graph->nbr[matrixIndex];
#endif
        ++matrixIndex;
      }
    }
    graph->nbrOffset[numVertex] = matrixIndex;
  }
  inFile.close();
}

inline void ReadGFUFormat(const string& graphFile, vector<Graph*>& graphList) {
  ifstream inFile(graphFile);
  if (!inFile) {
    cout << "Failed to open a data graph file " << graphFile << endl;
    exit(1);
  }

  Graph* graph = NULL;
#ifdef EDGE_LABEL
  vector<tuple<int, int, int>> edges;
#else
  vector<pair<int, int>> edges;
#endif
  string line;
  number_type pos, numVertex, numEdge;
  int fst, snd, labelID;
  int cntVertex = 0, state = Start;
  unordered_map<string, int>::const_iterator sIter;
  unordered_map<int, int>::const_iterator iIter;
  while (!inFile.eof()) {
    switch (state) {
      case Start: {
        inFile >> line;
        if (inFile.eof()) break;
        graph = new Graph();
        graphList.push_back(graph);
        cntVertex = 0;
        edges.clear();
        state = LoadNumVertex;
        break;
      }
      case LoadNumVertex: {
        inFile >> numVertex;
        graph->initVertex(numVertex);
        state = LoadVertex;
        break;
      }
      case LoadVertex: {
        inFile >> line;
        sIter = labelMap.find(line);
        if (sIter != labelMap.end()) {
          labelID = sIter->second;
        } else {
          labelID = labelMap.size();
          labelMap[line] = labelID;
        }
        iIter = graph->labelFrequency.find(labelID);
        if (iIter != graph->labelFrequency.end()) {
          graph->labelFrequency[labelID] += 1;
        } else {
          graph->labelFrequency[labelID] = 1;
        }
        graph->label[cntVertex] = labelID;
        ++cntVertex;
        if (cntVertex >= graph->nVertex) state = LoadNumEdge;
        break;
      }
      case LoadNumEdge: {
        inFile >> numEdge;
        graph->initEdge(numEdge);
        graph->nLabel = labelMap.size();
        edges.reserve(numEdge);
        state = LoadEdge;
        break;
      }
      case LoadEdge: {
#ifdef EDGE_LABEL
        inFile >> fst >> snd >> labelID;
        edges.push_back(make_tuple(fst, snd, labelID));
#else
        inFile >> fst >> snd;
        edges.push_back(make_pair(fst, snd));
#endif
        graph->degree[fst] += 1;
        graph->degree[snd] += 1;
        if (edges.size() >= graph->nEdge) {
          pos = 0;
          for (int i = 0; i < graph->nVertex; ++i) {
            graph->nbrOffset[i] = pos;
            pos += (number_type)graph->degree[i];
            graph->core[i] = graph->degree[i];
            if (graph->degree[i] > graph->maxDegree)
              graph->maxDegree = graph->degree[i];
          }
          graph->nbrOffset[graph->nVertex] = pos;
          memset(graph->degree, 0, sizeof(int) * graph->nVertex);
          for (int i = 0; i < edges.size(); ++i) {
#ifdef EDGE_LABEL
            fst = std::get<0>(edges[i]);
            snd = std::get<1>(edges[i]);
            labelID = std::get<2>(edges[i]);
            graph->nbr.set(graph->nbrOffset[fst] + graph->degree[fst], snd,
                           labelID);
            graph->nbr.set(graph->nbrOffset[snd] + graph->degree[snd], fst,
                           labelID);
#else
            fst = edges[i].first;
            snd = edges[i].second;
            graph->nbr[graph->nbrOffset[fst] + graph->degree[fst]] = snd;
            graph->nbr[graph->nbrOffset[snd] + graph->degree[snd]] = fst;
#endif
            graph->degree[fst] += 1;
            graph->degree[snd] += 1;
          }
          state = Start;
        }
        break;
      }
      default: {
        exit(1);
      }
    }
  }
  inFile.close();
  edges.clear();
  /*
  for(int j = 0; j < graphList.size(); ++j){
    Graph& g = *graphList[j];
    g.print(j);
  }
  */
}

inline void ReadIgraphFormat(const string& graphFile,
                             vector<Graph*>& graphList) {
  ifstream inFile(graphFile);
  if (!inFile) {
    cout << "Failed to open a data graph file " << graphFile << endl;
    exit(1);
  }

  Graph* graph = NULL;
  vector<pair<int, int>> edges;
  vector<int> vertices;
  string line;
  char symbol, state;
  int fst, snd, labelID, cntVertex;
  number_type pos, numVertex, numEdge;
  unordered_map<string, int>::const_iterator sIter;
  unordered_map<int, int>::const_iterator iIter;
  while (!inFile.eof()) {
    inFile >> state;
    if (inFile.eof()) break;
    switch (state) {
      case 't':
        if (graph != NULL) {
          graph->initEdge(edges.size());
          pos = 0;
          for (int i = 0; i < graph->nVertex; ++i) {
            graph->nbrOffset[i] = pos;
            pos += (number_type)graph->degree[i];
            graph->core[i] = graph->degree[i];
            if (graph->degree[i] > graph->maxDegree)
              graph->maxDegree = graph->degree[i];
          }
          graph->nbrOffset[graph->nVertex] = pos;
          memset(graph->degree, 0, sizeof(int) * graph->nVertex);

          for (int i = 0; i < edges.size(); ++i) {
            fst = edges[i].first;
            snd = edges[i].second;
#ifdef EDGE_LABEL
            graph->nbr.set(graph->nbrOffset[fst] + graph->degree[fst], snd);
            graph->nbr.set(graph->nbrOffset[snd] + graph->degree[snd], fst);
#else
            graph->nbr[graph->nbrOffset[fst] + graph->degree[fst]] = snd;
            graph->nbr[graph->nbrOffset[snd] + graph->degree[snd]] = fst;
#endif
            graph->degree[fst] += 1;
            graph->degree[snd] += 1;
          }
        }
        inFile >> symbol >> line;
        graph = new Graph();
        graphList.push_back(graph);
        cntVertex = 0;
        edges.clear();
        vertices.clear();
        break;
      case 'v':
        inFile >> fst >> line;
        sIter = labelMap.find(line);
        if (sIter != labelMap.end()) {
          labelID = sIter->second;
        } else {
          labelID = labelMap.size();
          labelMap[line] = labelID;
        }
        iIter = graph->labelFrequency.find(labelID);
        if (iIter != graph->labelFrequency.end()) {
          graph->labelFrequency[labelID] += 1;
        } else {
          graph->labelFrequency[labelID] = 1;
        }
        vertices.push_back(labelID);
        ++cntVertex;
        break;
      case 'e':
        if (edges.size() == 0) {
          graph->initVertex(cntVertex);
          graph->nLabel = labelMap.size();
          for (int i = 0; i < graph->nVertex; ++i)
            graph->label[i] = vertices[i];
        }
        inFile >> fst >> snd >> labelID;
        // inFile >> fst >> snd;
        graph->degree[fst] += 1;
        graph->degree[snd] += 1;
        edges.push_back(make_pair(fst, snd));
        break;
      default:
        exit(1);
        break;
    }
  }
  if (graph != NULL) {
    graph->initEdge(edges.size());
    pos = 0;
    for (int i = 0; i < graph->nVertex; ++i) {
      graph->nbrOffset[i] = pos;
      pos += (number_type)graph->degree[i];
      graph->core[i] = graph->degree[i];
      if (graph->degree[i] > graph->maxDegree)
        graph->maxDegree = graph->degree[i];
    }
    graph->nbrOffset[graph->nVertex] = pos;
    memset(graph->degree, 0, sizeof(int) * graph->nVertex);

    for (int i = 0; i < edges.size(); ++i) {
      fst = edges[i].first;
      snd = edges[i].second;
#ifdef EDGE_LABEL
      graph->nbr.set(graph->nbrOffset[fst] + graph->degree[fst], snd);
      graph->nbr.set(graph->nbrOffset[snd] + graph->degree[snd], fst);
#else
      graph->nbr[graph->nbrOffset[fst] + graph->degree[fst]] = snd;
      graph->nbr[graph->nbrOffset[snd] + graph->degree[snd]] = fst;
#endif
      graph->degree[fst] += 1;
      graph->degree[snd] += 1;
    }
  }
  inFile.close();
  edges.clear();
  vertices.clear();
  /*
  for(int j = 0; j < graphList.size(); ++j){
    Graph& g = *graphList[j];
    g.print(j);
  }
  */
}

inline void ProcessDataGraphs() {
  int start, num, maxValue;
  unordered_map<int, int>::const_iterator lit;
  nGraph = dataGraph.size();
  nUniqueLabel = labelMap.size();
#ifdef NEIGHBOR_LABEL_FREQUENCY
#else
  // NLFSize = BITS_PER_LABEL * (nUniqueLabel + 1) / INT_SIZE + 1;
  NLFSize = (BITS_PER_LABEL * nUniqueLabel - 1) / UINT64_SIZE + 1;
#endif
  for (int i = 0; i < dataGraph.size(); ++i) {
    Graph& g = *dataGraph[i];
    if (g.maxDegree > maxDegree) maxDegree = g.maxDegree;
    if (g.nVertex > maxNumDataVertex) maxNumDataVertex = g.nVertex;
  }
  positiveLabel = new int[nUniqueLabel];
  cntLabel = new int[nUniqueLabel]();
  int* bin = new int[maxDegree + 1];
  int* pos = new int[maxNumDataVertex];
  int* vert = new int[maxNumDataVertex];
  // cout<<"N(unique label): "<<nUniqueLabel<<endl;
  // cout<<"NLF size (bytes): "<<NLFSize<<endl;
  //  For each data graph g
  for (int i = 0; i < dataGraph.size(); ++i) {
    start = 0;
    maxValue = 0;
    memset(bin, 0, sizeof(int) * (maxDegree + 1));
    Graph& g = *dataGraph[i];
    for (lit = g.labelFrequency.begin(); lit != g.labelFrequency.end(); ++lit) {
      if (lit->second > maxValue) {
        g.mostFrequentLabelID = lit->first;
        maxValue = lit->second;
      }
    }
    g.maxNumSameLabelVertex = maxValue;
    if (maxValue > maxNumCandidate) maxNumCandidate = maxValue;
    // Compute core numbers
    for (int i = 0; i < g.nVertex; i++) bin[g.core[i]]++;

    for (int d = 0; d <= g.maxDegree; d++) {
      num = bin[d];
      bin[d] = start;
      start += num;
    }
    for (int i = 0; i < g.nVertex; i++) {
      pos[i] = bin[g.core[i]];
      vert[pos[i]] = i;
      bin[g.core[i]]++;
    }
    for (int d = g.maxDegree; d > 0; d--) bin[d] = bin[d - 1];
    bin[0] = 0;

    for (int i = 0; i < g.nVertex; i++) {
      int v = vert[i];
      for (number_type j = g.nbrOffset[v]; j < g.nbrOffset[v + 1]; j++) {
        int u = g.nbr[j];
        if (g.core[u] > g.core[v]) {
          int du = g.core[u];
          int pu = pos[u];
          int pw = bin[du];
          int w = vert[pw];
          if (u != w) {  // if not the same node, switch the position of the two
                         // nodes.
            pos[u] = pw;
            pos[w] = pu;
            vert[pu] = w;
            vert[pw] = u;
          }
          bin[du]++;
          g.core[u]--;
        }
      }
    }
    // Compute Neighbor Label Frequency
    g.sortNeighbors();
    g.computeNLF();
  }
  delete[] bin;
  delete[] pos;
  delete[] vert;
}

inline void GetQueryDataStructure(Graph& query) {
  // When the core decomposition begins, core number is the degree
  int* pos = posQuery;    // int pos[nQueryVertex];
  int* vert = vertQuery;  // int vert[nQueryVertex];
  int* bin = binQuery;    // int bin[maxDegreeQuery + 1];
  memset(bin, 0, sizeof(int) * (query.maxDegree + 1));
  for (int i = 0; i < nQueryVertex; i++) bin[query.core[i]]++;

  int num, start = 0;
  for (int d = 0; d <= query.maxDegree; d++) {
    num = bin[d];
    bin[d] = start;
    start += num;
  }
  for (int i = 0; i < nQueryVertex; i++) {
    pos[i] = bin[query.core[i]];
    vert[pos[i]] = i;
    bin[query.core[i]]++;
  }

  for (int d = query.maxDegree; d > 0; d--) bin[d] = bin[d - 1];
  bin[0] = 0;

  for (int i = 0; i < nQueryVertex; i++) {
    int v = vert[i];
    for (number_type j = query.nbrOffset[v]; j < query.nbrOffset[v + 1]; j++) {
      int u = query.nbr[j];
      if (query.core[u] > query.core[v]) {
        int du = query.core[u];
        int pu = pos[u];
        int pw = bin[du];
        int w = vert[pw];
        if (u != w) {
          // if not the same node, switch the position of the two nodes.
          pos[u] = pw;
          pos[w] = pu;
          vert[pu] = w;
          vert[pw] = u;
        }
        bin[du]++;
        query.core[u]--;
      }
    }
  }
  query.computeNLF();
}
