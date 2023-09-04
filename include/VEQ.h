#define CELL_COMPUTATION_OPT
#include <algorithm>
#include <cfloat>
#define __STDC_FORMAT_MACROS

#include <unordered_map>

#include "memory.h"
using namespace std;

// Elapsed time variable used in the main function
bool exceedTimeLimit;
double filteringTime;
double verificationTime;
double totalTime;

// File paths
string dataGraphFile;
string queryGraphFile;
string answerFile;

// Temporary variable
uint64_t* b;

struct dynamicDAG_ancestors {
 private:
  uint64_t dag_ancestors[MAX_NUM_VERTEX][MAX_NUM_VERTEX][MAX_NUM_VERTEX];
  // vertex depth partition
 public:
  void clear(int root) {
    fill(dag_ancestors[root][0], dag_ancestors[root][0] + FAILING_SET_SIZE,
         0ULL);
  }

  void addAncestor(int u, int uP, int uNumPar, int uPNumPar) {
    if (uNumPar == 1) {
      fill(dag_ancestors[u][1], dag_ancestors[u][1] + FAILING_SET_SIZE, 0ULL);
    }
    for (int x = 0; x < FAILING_SET_SIZE; x++) {
      dag_ancestors[u][uNumPar][x] = dag_ancestors[u][uNumPar - 1][x];
      dag_ancestors[u][uNumPar][x] |= dag_ancestors[uP][uPNumPar][x];
    }
    addInFailingSet(uP, dag_ancestors[u][uNumPar]);
  }

  uint64_t getSetPartition(int u, int uNumPar, int y) const {
    return dag_ancestors[u][uNumPar][y];
  }
} dyanc;

inline bool FilterByCount(const Graph& query, const Graph& data) {
  if (query.nVertex > data.nVertex || query.nEdge > data.nEdge ||
      query.maxDegree > data.maxDegree ||
      query.labelFrequency.size() > data.labelFrequency.size())
    return false;
  unordered_map<int, int>::const_iterator qIter, dIter;
  for (qIter = query.labelFrequency.begin();
       qIter != query.labelFrequency.end(); ++qIter) {
    dIter = data.labelFrequency.find(qIter->first);
    if (dIter == data.labelFrequency.end())
      return false;
    else if (qIter->second > dIter->second)
      return false;
  }
  return true;
}

inline int SelectRootVertex(const Graph& query, const Graph& data) {
  double leastRank = DBL_MAX, rank;
  int root = -1, label, degree;
  for (int i = 0; i < nQueryVertex; i++) {
    if (!isTree && query.core[i] < 2) continue;
    label = query.label[i];
    degree = query.degree[i];
    int start = data.vertexOffset[label];
    int end = data.vertexOffset[label + 1];
    while (start < end && degree > data.degree[data.vertex[start]]) start++;
    rank = (double)(end - start) / (double)degree;
    if (rank < leastRank) {
      leastRank = rank;
      root = i;
    }
  }
  return root;
}

inline void FindNEC(Graph& query, int u) {
  pair<int, int> p;
  unordered_map<pair<int, int>, int, pairHash>::iterator sizeIter;
  unordered_map<int, int>::const_iterator freqIter;
  // for all of node i's children
  for (int j = query.nbrOffset[u]; j < query.nbrOffset[u + 1]; j++) {
    int child = query.nbr[j];
    // If child is not in the two-core
    if (query.core[child] < 2) {
      // If child has degree one (i.e., NEC vertex)
      if (query.degree[child] == 1) {
        int label = query.label[child];
        p = make_pair(label, u);
        sizeIter = NECSize.find(p);
        if (sizeIter == NECSize.end()) {
          // child is the representative node
          query.NECElement[query.numNECMapping++] = Element(label, u, child);
          query.NECMap[child] = child;
          NECChild[p] = child;
          NECSize[p] = 1;
          freqIter = query.labelFrequency.find(label);
          if (freqIter != query.labelFrequency.end() && freqIter->second > 1)
            query.isProblemLeaf[child] = true;
        } else {
          query.NECMap[child] = NECChild[p];
          sizeIter->second += 1;
        }
        // the label with parent being i, necCount++
      }
      // If child is in tree
      else {
        // extract the query tree for extra region candidate extraction, based
        // on DFS also give a DFS-based query sequence at the same time
        int* stack = arrForQuery;
        int posStack = 0;
        visitedForQuery[u] =
            1;  // this is the start node's parent node (a marked node)
        visitedForQuery[child] = 1;  // this is the start node
        stack[posStack++] = child;
        while (posStack != 0) {
          int curr = stack[posStack - 1];
          posStack--;
          for (int m = query.nbrOffset[curr]; m < query.nbrOffset[curr + 1];
               m++) {
            int child = query.nbr[m];
            if (!visitedForQuery[child]) {
              visitedForQuery[child] = 1;
              // ======== special treatment here: if a node is a leaf (degree
              // being 1), then put it into nec node set
              if (query.degree[child] == 1) {
                int label = query.label[child];
                p = make_pair(label, curr);
                sizeIter = NECSize.find(p);
                if (sizeIter == NECSize.end()) {
                  // child is the repesentive node
                  query.NECElement[query.numNECMapping++] =
                      Element(label, curr, child);
                  query.NECMap[child] = child;
                  NECChild[p] = child;
                  NECSize[p] = 1;
                  freqIter = query.labelFrequency.find(label);
                  if (freqIter != query.labelFrequency.end() &&
                      freqIter->second > 1)
                    query.isProblemLeaf[child] = true;
                } else {
                  query.NECMap[child] = NECChild[p];
                  sizeIter->second += 1;
                  // increase the NEC size of vertices with 'label' and parent
                  // 'curr'
                }
                continue;
              }
              stack[posStack++] = child;
            }
          }
        }
      }
    }
  }
}

inline void ExtractResidualStructures(Graph& query) {
  NECSize.clear();
  NECChild.clear();
  NECElementOffset.clear();
  query.numNECMapping = 0;
  memset(query.NECMap, -1, sizeof(int) * nQueryVertex);
  memset(query.isProblemLeaf, false, sizeof(bool) * nQueryVertex);
  if (isTree) {
    FindNEC(query, root);
  } else {
    for (int i = 0; i < nQueryVertex; i++) {
      if (query.core[i] < 2) continue;  // not in two-core
      FindNEC(query, i);
    }
  }

  sort(query.NECElement, query.NECElement + query.numNECMapping,
       sortByNECLabel);
  // NECElement: a list of (label: lID, parentID: uP, representative: u,
  // sum: |NEC(u)|) sorted in the ascending order of label ID

  int sum;
  int lastLabel;
  int NECSetSize = 0;
  if (query.numNECMapping == 1) {
    Element& elem = query.NECElement[0];
    int label = elem.label;
    int parentID = elem.parentID;
    int representChild = elem.representative;
    sum = NECSize[make_pair(label, parentID)];
    NECElementOffset.push_back(make_pair(label, NECSetSize));
    query.NECElement[0].sum = sum;
    NECElementOffset.push_back(
        make_pair(-1, query.numNECMapping));  // redundant currE to set the end
  } else {
    for (int i = 0; i < query.numNECMapping; i++) {
      Element& elem = query.NECElement[i];
      int label = elem.label;
      int parentID = elem.parentID;
      int representChild = elem.representative;
      pair<int, int> p = make_pair(label, parentID);
      sum = NECSize[make_pair(label, parentID)];
      if (i == 0) {
        NECElementOffset.push_back(make_pair(label, NECSetSize));
        query.NECElement[NECSetSize++].sum = sum;
        lastLabel = label;
        continue;
      } else if (i == query.numNECMapping - 1) {
        if (label != lastLabel)
          NECElementOffset.push_back(make_pair(label, NECSetSize));
        query.NECElement[NECSetSize++].sum = sum;
        NECElementOffset.push_back(make_pair(
            -1, query.numNECMapping));  // redunant currE to set the end
        continue;
      }
      if (label != lastLabel) {
        NECElementOffset.push_back(make_pair(label, NECSetSize));
        lastLabel = label;
      }
      query.NECElement[NECSetSize++].sum = sum;
    }
  }
  // NECElemtntOffset[i] = (l, pos)
  // l: label ID
  // pos: the starting position of NEC vertices with label l in NECElement if
  // l != -1, the last position of NECElement otherwise
}

inline void BuildDAG(Graph& query) {
  char* uCurr = visitedForQuery;  // visitedForQuery[i] = 1 if i-node have uCurr
                                  // from queue.
  memset(uCurr, 0, sizeof(char) * nQueryVertex);
  int* visited =
      visitedForQueryInt;  // visited[i] = level if i-node had pushed into
                           // queue, where level is it's BFS-level.
  memset(visited, 0, sizeof(int) * nQueryVertex);
  memset(DAG_child_query_size, 0, sizeof(int) * nQueryVertex);
  memset(DAG_parent_query_size, 0, sizeof(int) * nQueryVertex);
  memset(DAG_ngb_query_size, 0, sizeof(int) * nQueryVertex);

  arrForQuery[0] = root;
  int startPtr = 0, endPtr = 1, nextEndPtr = 1;
  uSequenceSize = 0;
  visited[root] = 1;
  // Sort labels in the ascending order of label frequency

  while (true) {
    //<sort by label_freq first, and sort by degree for each label group>
    sort(arrForQuery + startPtr, arrForQuery + endPtr, sortQueryByDegree);
    stable_sort(arrForQuery + startPtr, arrForQuery + endPtr,
                sortQueryByLabelFrequency);
    while (startPtr != endPtr) {
      int currNode = arrForQuery[startPtr++];
      uCurr[currNode] = 1;
      for (int i = query.nbrOffset[currNode]; i < query.nbrOffset[currNode + 1];
           ++i) {
        int childNode = query.nbr[i];
        if (uCurr[childNode] == 0)  // childNode is not currNode's parent
        {
          {  //[1]. child and parent
            DAG_parent_query[childNode][DAG_parent_query_size[childNode]] =
                currNode;
            DAG_child_query_size[currNode]++;
            DAG_parent_query_size[childNode]++;
          }
          {
            int& currIdx = DAG_ngb_query_size[currNode];
            int& childIdx = DAG_ngb_query_size[childNode];

            DAG_ngb_query[currNode][currIdx] = childNode;
            DAG_ngb_query[childNode][childIdx] = currNode;

            DAG_ngb_query_ngb_index[currNode][currIdx] = childIdx;
            DAG_ngb_query_ngb_index[childNode][childIdx] = currIdx;

            ++currIdx;
            ++childIdx;
          }
        }
        if (visited[childNode] == 0)  // this child node has not been visited.
        {
          visited[childNode] = 1;
          arrForQuery[nextEndPtr] = childNode;
          nextEndPtr++;
          if (query.core[childNode] < 2) continue;
        }
      }
      orderDAG[currNode] = uSequenceSize;
      uSequence[uSequenceSize++] = currNode;
    }
    if (nextEndPtr == endPtr) break;
    startPtr = endPtr;
    endPtr = nextEndPtr;
  }
  query.nNonLeaf = 0;
  for (int i = 0; i < nQueryVertex; ++i) {
    if (query.NECMap[i] == -1) ++query.nNonLeaf;
  }
}

inline bool filteringWithDAG(const Graph& query, const Graph& data,
                             direction filterDir, ngbType filterNgb,
                             bool useNbrSafety) {
  unordered_map<pair<int, int>, pair<number_type, number_type>,
                pairHash>::const_iterator iter;

  int orderSt, orderEd, orderAdd;

  if (filterDir == topDown) {
    orderSt = 0;
    orderEd = uSequenceSize;
    orderAdd = 1;
  } else {
    orderSt = uSequenceSize - 1;
    orderEd = -1;
    orderAdd = -1;
  }
  for (int i = orderSt; i != orderEd; i += orderAdd) {
    ++index_vis_adj;
    ++index_vis_color;
    const int query_node_cur = uSequence[i];
    if (query.NECMap[query_node_cur] != -1 &&
        query.NECMap[query_node_cur] != query_node_cur)
      continue;
    const int color_cur = query.label[query_node_cur];
    index_set_color = 0;
    int s = query.nbrOffset[query_node_cur];
    int d = query.degree[query_node_cur];

    // 1.1 cnt the number of color
    for (int idx_node_ngb = s; idx_node_ngb < s + d; idx_node_ngb++) {
      const int query_node_ngb = query.nbr[idx_node_ngb];
      const int color_ngb = query.label[query_node_ngb];

      if (color[color_ngb][0] == index_vis_color) {
        color[color_ngb][1]++;
      } else {
        color[color_ngb][0] = index_vis_color;
        color[color_ngb][1] = 1;
        set_color[index_set_color++] = color_ngb;
      }

      if (query.NECMap[query_node_ngb] != -1 &&
          query.NECMap[query_node_ngb] != query_node_ngb) {
        continue;
      }

      // 1.1.2 see the whole cs of neighbor
      const CandidateSpace& node_unit_ngb = candSpace[query_node_ngb];
      for (int y = 0; y < node_unit_ngb.size; y++) {
        const int cand_ngb = node_unit_ngb.candidates[y];
        if (vis_adj[cand_ngb] != index_vis_adj) {
          vis_adj[cand_ngb] = index_vis_adj;
        }
      }

      // 1.1.3 filtering
    }

    // 1.2 calculate parent existence part
    ++index_ngb_existence;
    int nVisitedNeighbor = 0;
    for (int y = 0; y < query.degree[query_node_cur];
         y++)  // Filtering rule 3.1
    {
      const int query_node_nbr =
          query.nbr[query.nbrOffset[query_node_cur] + y];  // Filtering rule 3.1

      if (filterNgb == parentNgb &&
          orderDAG[query_node_cur] < orderDAG[query_node_nbr]) {
        continue;
      }
      if (filterNgb == childNgb &&
          orderDAG[query_node_cur] > orderDAG[query_node_nbr]) {
        continue;
      }

      if (query.NECMap[query_node_nbr] != -1 &&
          query.NECMap[query_node_nbr] != query_node_nbr)
        continue;

      CandidateSpace& nbr_unit = candSpace[query_node_nbr];
      for (int z = 0; z < nbr_unit.size; z++) {
        int cand_ngb = nbr_unit.candidates[z];
#ifdef HUGE_GRAPH
        for (long long w = data.nbrOffset[cand_ngb];
             w < data.nbrOffset[cand_ngb + 1]; w++)
#else
        iter = data.labelToNbrOffset.find(make_pair(cand_ngb, color_cur));
        if (iter == data.labelToNbrOffset.end()) continue;
        for (int w = iter->second.first;
             w < iter->second.first + iter->second.second; w++)
#endif
        {
          int cand_around = data.nbr[w];
#ifdef HUGE_GRAPH
          if (data.label[cand_around] != color_cur) continue;
#endif
          if (nVisitedNeighbor == 0) {
            num_ngb_existence[cand_around][0] = index_ngb_existence;
            num_ngb_existence[cand_around][1] = 1;
          } else {
            if (num_ngb_existence[cand_around][0] == index_ngb_existence &&
                num_ngb_existence[cand_around][1] == nVisitedNeighbor) {
              num_ngb_existence[cand_around][1]++;
            }
          }
        }
      }
      ++nVisitedNeighbor;
    }

    // 1.3 filtering
    CandidateSpace& node_unit_cur = candSpace[query_node_cur];
    int new_size = node_unit_cur.size;
    int x = 0;
    while (x < new_size) {
      const int cand_cur = node_unit_cur.candidates[x];
      bool iSVaildCand = true;
      // 1.3.1 : filtering by parent

      // Filtering rule 3.1 (in the codition below after ||)
      if ((num_ngb_existence[cand_cur][0] != index_ngb_existence ||
           num_ngb_existence[cand_cur][1] != nVisitedNeighbor) &&
          (nVisitedNeighbor != 0)) {
        iSVaildCand = false;
      }

      // 1.3.2 : filtering by cardinality
      if (iSVaildCand && useNbrSafety) {
#ifdef HUGE_GRAPH
        int t = 0;
        for (long long z = data.nbrOffset[cand_cur];
             z < data.nbrOffset[cand_cur + 1]; ++z) {
          int ngb = data.nbr[z];
          if (vis_adj[ngb] == index_vis_adj) {
            ++t;
          }
        }
        if (t < d) {
          iSVaildCand = false;
        }
#else
        for (int y = 0; y < index_set_color; y++) {
          int t = 0;
          iter = data.labelToNbrOffset.find(make_pair(cand_cur, set_color[y]));
          if (iter == data.labelToNbrOffset.end()) {
            iSVaildCand = false;
            break;
          }
          for (int z = iter->second.first;
               z < iter->second.first + iter->second.second; ++z) {
            int ngb = data.nbr[z];

            if (vis_adj[ngb] == index_vis_adj) {
              ++t;
            }
          }
          if (t < color[set_color[y]][1]) {
            iSVaildCand = false;
            break;
          }
        }
#endif
      }
      // 1.3.3 : remain or remove
      if (!iSVaildCand) {
        int y = node_unit_cur.candidates[new_size - 1];
        node_unit_cur.candidates[x] = y;
        --new_size;
      } else {
        ++x;
      }
    }
    if (new_size == 0) return false;
    node_unit_cur.size = new_size;
  }
  return true;
}

inline bool TopDownInitial(const Graph& query, const Graph& data) {
  // Initialize candidates of root vertex
  CandidateSpace& rootUnit = candSpace[root];
  int label = query.label[root];
  int start = data.vertexOffset[label];
  int end = data.vertexOffset[label + 1];
  unordered_map<pair<int, int>, pair<number_type, number_type>,
                pairHash>::const_iterator iter;
  while (start < end && query.degree[root] > data.degree[data.vertex[start]])
    start++;
  // For each v in G with label
  for (int j = start; j < end; j++) {
    int v = data.vertex[j];  // v
    if (data.core[v] < query.core[root] ||
        data.maxNbrDegree[v] < query.maxNbrDegree[root])
      continue;
    char isAdded = 1;
    for (int pos = NLFSize - 1; pos >= 0; pos--) {
#ifdef HUGE_GRAPH
      if (data.NLF[(long long)v * (long long)NLFSize + (long long)pos] !=
          (query.NLF[root * NLFSize + pos] |
           data.NLF[(long long)v * (long long)NLFSize + (long long)pos]))
#else
      if (data.NLF[v * NLFSize + pos] !=
          (query.NLF[root * NLFSize + pos] | data.NLF[v * NLFSize + pos]))
#endif
      {
        isAdded = 0;
        break;
      }
    }
    if (isAdded) rootUnit.candidates[rootUnit.size++] = v;
  }
  if (rootUnit.size == 0) return false;
  // For each query vertex u of q in a top-down fashion
  for (int i = 1; i < uSequenceSize; i++) {
    // i = 0 is root, which is already processed above.
    int currVertex = uSequence[i];  // u
    if (query.NECMap[currVertex] != -1 &&
        query.NECMap[currVertex] != currVertex)
      continue;
    int currLabel = query.label[currVertex];
    char checkVal = 0;  // Cnt = 0

    // For each parent p of u in dag
    for (int x = 0; x < DAG_parent_query_size[currVertex]; x++) {
      int nbr = DAG_parent_query[currVertex][x];  // p
      CandidateSpace& nbrUnit = candSpace[nbr];
      // For each candidate vertex v' in p.C
      for (int y = 0; y < nbrUnit.size; y++) {
        int parentCand = nbrUnit.candidates[y];  // v'
// For each vertex v adjacent to v' with label
#ifdef HUGE_GRAPH
        for (long long z = data.nbrOffset[parentCand];
             z < data.nbrOffset[parentCand + 1]; ++z) {
          int canID = data.nbr[z];
          if (data.label[canID] != currLabel) continue;
          if (isVisited[canID] == checkVal) {
            isVisited[canID]++;  // v.cnt++;
            if (checkVal == 0) arrayToClean[toCleanIndex++] = canID;
          }
        }
#else
        iter = data.labelToNbrOffset.find(make_pair(parentCand, currLabel));
        if (iter == data.labelToNbrOffset.end()) continue;
        for (int z = iter->second.first;
             z < iter->second.first + iter->second.second; ++z) {
          int canID = data.nbr[z];
          if (isVisited[canID] == checkVal) {
            isVisited[canID]++;  // v.cnt++;
            if (checkVal == 0) arrayToClean[toCleanIndex++] = canID;
          }
        }
#endif
      }
      checkVal++;  // update the check value by one
    }
    // For each vertex v with v.cnt = Cnt, if v passes filters, add candidates
    CandidateSpace& currSet = candSpace[currVertex];
    int candIndex = 0;
    for (int checkIndex = 0; checkIndex < toCleanIndex; checkIndex++) {
      int canID = arrayToClean[checkIndex];  // v
      if (isVisited[canID] == checkVal) {
        // check degree, core, and max_neighbor degree together
        if (data.degree[canID] < query.degree[currVertex] ||
            data.core[canID] < query.core[currVertex] ||
            data.maxNbrDegree[canID] < query.maxNbrDegree[currVertex])
          continue;
        // check lightweight NLF
        char isAdded = 1;
        for (int pos = NLFSize - 1; pos >= 0; pos--) {
#ifdef HUGE_GRAPH
          if (data.NLF[(long long)canID * (long long)NLFSize +
                       (long long)pos] !=
              (query.NLF[currVertex * NLFSize + pos] |
               data.NLF[(long long)canID * (long long)NLFSize +
                        (long long)pos]))
#else
          if (data.NLF[canID * NLFSize + pos] !=
              (query.NLF[currVertex * NLFSize + pos] |
               data.NLF[canID * NLFSize + pos]))
#endif
          {
            isAdded = 0;
            break;
          }
        }
        // lightweight NLF OK
        if (isAdded) {
          currSet.candidates[candIndex] = canID;
          candIndex++;
        }
      }
    }
    currSet.size = candIndex;
    // Reset v.cnt for all v in G such that v.cnt > 0
    while (toCleanIndex != 0) isVisited[arrayToClean[--toCleanIndex]] = 0;
    if (currSet.size == 0) return false;
  }  // In a top-down manner
  return true;
}

inline bool BottomUp(const Graph& query, const Graph& data) {
  unordered_map<pair<int, int>, pair<number_type, number_type>,
                pairHash>::const_iterator iter;
  // For each query vertex u of dag(q) in a bottom-up fashion
  for (int i = uSequenceSize - 1; i >= 0; --i) {
    int currVertex = uSequence[i];  // u
    // No child, No prunning.
    if (DAG_child_query_size[currVertex] == 0) continue;
    int currLabel = query.label[currVertex];
    char checkVal = 0;  // Cnt = 0
    // For each child uc of u in dag(q)
    for (int x = query.nbrOffset[currVertex];
         x < query.nbrOffset[currVertex + 1]; x++) {
      int nbr = query.nbr[x];  // uc
      CandidateSpace& nbrSet = candSpace[nbr];
      // process only represent node when the child is in NEC
      if (query.NECMap[nbr] != -1 && query.NECMap[nbr] != nbr) continue;
      // For each candidate vc in uc.C
      for (int y = 0; y < nbrSet.size; y++) {
        int childCand = nbrSet.candidates[y];  // vc
// For each vertex v adjacent to vc
#ifdef HUGE_GRAPH
        for (long long z = data.nbrOffset[childCand];
             z < data.nbrOffset[childCand + 1]; ++z)
#else
        iter = data.labelToNbrOffset.find(make_pair(childCand, currLabel));
        if (iter == data.labelToNbrOffset.end()) continue;
        for (int z = iter->second.first;
             z < iter->second.first + iter->second.second; ++z)
#endif
        {
          int canID = data.nbr[z];  // v
#ifdef HUGE_GRAPH
          if (data.label[canID] != currLabel) continue;
#endif
          if (isVisited[canID] == checkVal) {
            isVisited[canID]++;
            if (checkVal == 0) arrayToClean[toCleanIndex++] = canID;
          }
        }
      }
      checkVal++;  // update the check value by one
    }
    // For each vertex v in u.C
    CandidateSpace& currSet = candSpace[currVertex];
    for (int candIndex = 0; candIndex < currSet.size; candIndex++) {
      int canID = currSet.candidates[candIndex];  // v
      // if v.cnt != Cnt then Remove v from u.C
      if (isVisited[canID] != checkVal) {
        currSet.candidates[candIndex] = currSet.candidates[--currSet.size];
        --candIndex;
      }
    }
    // Reset v.cnt = 0 for all v in G such that v.cnt > 0
    while (toCleanIndex != 0) isVisited[arrayToClean[--toCleanIndex]] = 0;
    if (currSet.size == 0) return false;
  }  // In a bottom-up manner
  return true;
}

inline bool TopDown(const Graph& query, const Graph& data) {
  unordered_map<pair<int, int>, pair<number_type, number_type>,
                pairHash>::const_iterator iter;
  // For each query vertex u of dag(q) in a top-down fashion
  for (int i = 1; i < uSequenceSize; i++) {
    // root(0) has no parent
    int currVertex = uSequence[i];  // u
    if (query.NECMap[currVertex] != -1 &&
        query.NECMap[currVertex] != currVertex)
      continue;
    int currLabel = query.label[currVertex];
    char checkVal = 0;  // Cnt = 0
    // For each parent up of u in dag(q)
    for (int x = query.nbrOffset[currVertex];
         x < query.nbrOffset[currVertex + 1]; x++) {
      int nbr = query.nbr[x];  // up
      CandidateSpace& nbrUnit = candSpace[nbr];
      if (query.NECMap[nbr] != -1 && query.NECMap[nbr] != nbr) continue;
      // For each candidate vp in up.C
      for (int y = 0; y < nbrUnit.size; y++) {
        int parentCand = nbrUnit.candidates[y];  // vp
// For each vertex v adjacent to vp
#ifdef HUGE_GRAPH
        for (long long z = data.nbrOffset[parentCand];
             z < data.nbrOffset[parentCand + 1]; ++z) {
          int canID = data.nbr[z];  // v
          if (data.label[canID] != currLabel) continue;
          if (isVisited[canID] == checkVal) {
            isVisited[canID]++;
            if (checkVal == 0) arrayToClean[toCleanIndex++] = canID;
          }
        }
#else
        iter = data.labelToNbrOffset.find(make_pair(parentCand, currLabel));
        if (iter == data.labelToNbrOffset.end()) continue;
        for (int z = iter->second.first;
             z < iter->second.first + iter->second.second; ++z) {
          int canID = data.nbr[z];  // v
          if (isVisited[canID] == checkVal) {
            isVisited[canID]++;
            if (checkVal == 0) arrayToClean[toCleanIndex++] = canID;
          }
        }
#endif
      }
      checkVal++;  // update the check value by one
    }
    // For each vertex v in u.C
    CandidateSpace& currSet = candSpace[currVertex];
    for (int candIndex = 0; candIndex < currSet.size; candIndex++) {
      int canID = currSet.candidates[candIndex];  // v
      // if v.cnt != Cnt then Remove v from u.C
      if (isVisited[canID] != checkVal) {
        currSet.candidates[candIndex] = currSet.candidates[--currSet.size];
        --candIndex;
      }
    }
    // Reset v.cnt = 0 for all v in G such that v.cnt > 0
    while (toCleanIndex != 0) isVisited[arrayToClean[--toCleanIndex]] = 0;
    if (currSet.size == 0) return false;
  }  // In a top-down manner
  return true;
}

inline bool ConstructAdjacencyList(const Graph& query, const Graph& data,
                                   bool filter) {
  unordered_map<pair<int, int>, pair<number_type, number_type>,
                pairHash>::const_iterator iter;
  // For each query vertex u of dag(q) in a bottom-up fashion
  for (int i = uSequenceSize - 1; i >= 0; --i) {
    int u = uSequence[i];  // u
    int currLabel = query.label[u];
    if (query.NECMap[u] != -1 && query.NECMap[u] != u) continue;
    CandidateSpace& currSet = candSpace[u];
    // No child, No adjacency list (adjacency list is stored in the child)
    if (filter) {
      // <GlobalLabelFrequency>
      for (int vPos = 0; vPos < currSet.size; ++vPos) {
        int currCand = currSet.candidates[vPos];
        if (isVisited[currCand] == 0) {
          isVisited[currCand] = 1;
          arrayToClean[toCleanIndex++] = currCand;
          ++cntLabel[currLabel];
        }
      }
      // </GlobalLabelFrequency>
    }
    // if (DAG_child_query_size[u] == 0) continue;

    // Preprocess u.C 1) to check (if v in u.C) in constant time 2) to know v's
    // posCand
    for (int vPos = 0; vPos < currSet.size; ++vPos) {
      candToPos[currSet.candidates[vPos]] = vPos;
    }
    // 1.1. for each neighbor un of u in q
    for (int unPos = 0; unPos < DAG_ngb_query_size[u]; ++unPos) {
      int un = DAG_ngb_query[u][unPos];  // uc
      CandidateSpace& nbrSet = candSpace[un];
      int uPos = DAG_ngb_query_ngb_index[u][unPos];

      // Process only represent node when the child is in NEC
      if (query.NECMap[un] != -1 && query.NECMap[un] != un) continue;
      // Initialize candidate space
      // Make sure it won't allocate new array every time
      memset(nbrSet.nAdjacent[uPos], 0, sizeof(int) * currSet.size);
      for (int vPos = 0; vPos < currSet.size; ++vPos) {
        int currSize = data.degree[currSet.candidates[vPos]];
        if (nbrSet.capacity[uPos][vPos] < currSize) {
          nbrSet.capacity[uPos][vPos] = currSize;
          nbrSet.adjacent[uPos][vPos] = new int[currSize];
        }
      }
      // 1.1.1. for each vertex vn in C(un) Most time-consuming in
      // ConstructAdjacencyList
      for (int vnPos = 0; vnPos < nbrSet.size; ++vnPos) {
        int vn = nbrSet.candidates[vnPos];
// For each v adjacent to vc in G with label
#ifdef HUGE_GRAPH
        for (long long z = data.nbrOffset[vn]; z < data.nbrOffset[vn + 1];
             ++z) {
          // For each vertex v in N_G(v_p) with label l_q(u)
          int v = data.nbr[z];  // v
          if (data.label[v] != currLabel) continue;
          int vPos = candToPos[v];
          if (vPos < 0) continue;
          nbrSet.adjacent[uPos][vPos][nbrSet.nAdjacent[uPos][vPos]++] = vnPos;
        }
#else
        iter = data.labelToNbrOffset.find(make_pair(vn, currLabel));
        if (iter == data.labelToNbrOffset.end()) continue;
        for (int z = iter->second.first;
             z < iter->second.first + iter->second.second; ++z) {
          // For each vertex v in N_G(v_p) with label l_q(u)
          int v = data.nbr[z];  // v
          int vPos = candToPos[v];
          if (vPos < 0) continue;
          nbrSet.adjacent[uPos][vPos][nbrSet.nAdjacent[uPos][vPos]++] = vnPos;
        }
#endif
      }
    }
    // Reset v in G such that v in u.C and v.index = v's index in u.C
    for (int vPos = 0; vPos < currSet.size; ++vPos) {
      candToPos[currSet.candidates[vPos]] = -1;
    }
  }
  if (filter) {
    // <GlobalLabelFrequency>
    while (toCleanIndex != 0) isVisited[arrayToClean[--toCleanIndex]] = 0;
    for (unordered_map<int, int>::const_iterator it =
             query.labelFrequency.begin();
         it != query.labelFrequency.end(); ++it) {
      if (cntLabel[it->first] < it->second) return false;
    }
    // </GlobalLabelFrequency>
  }
  return true;
}

inline void GetCell(const Graph& query, int u, int d, int* extCand,
                    int nExtCand) {
  int uPos, nbr, nCand = 0;
  CandidateSpace& currSet = candSpace[u];
  if (u == root) {
    for (int j = 0; j < currSet.size; ++j) candArray[j] = j;
    nCand = currSet.size;
  } else {
    int curr = -1;
    for (int j = 0; j < nExtCand; ++j) {
      if (currSet.cell[extCand[j]] < 0 && curr != extCand[j])
        candArray[nCand++] = extCand[j];
      curr = extCand[j];
    }
  }
  if (nCand == 0) return;
  // cout<<"Bundling up vertices in C(u"<<u<<")"<<endl;
  uI = u;
  nCandOffset = 0;
  candOffset[nCandOffset++] = 0;
  candOffset[nCandOffset++] = nCand;
  int nNbr = DAG_ngb_query_size[u];
  int index = 0;
  if (nNbr > 0) {
    nbr = DAG_ngb_query[u][index];
    while (query.NECMap[nbr] != -1 && query.NECMap[nbr] != nbr) {
      ++index;
      if (index >= nNbr) break;
      nbr = DAG_ngb_query[u][index];
    }
    if (index < nNbr) {
      uPos = DAG_ngb_query_ngb_index[u][index];
      for (int x = 0; x < nCand; ++x) posArray[candArray[x]] = 0;
      int nPrevCandOffset = nCandOffset;
      for (int j = 0; j < nPrevCandOffset - 1; ++j) {
        if (candOffset[j + 1] - candOffset[j] == 1) {
          currSet.cell[candArray[candOffset[j]]] = 0;
          currSet.cellVertex[currSet.nCellVertex++] = candArray[candOffset[j]];
          continue;
        }
        // Compute cells in C(u)
        sortByNeighbor(candArray, index, uPos, nbr, 0, candOffset[j],
                       candOffset[j + 1]);
      }
      // Sort the offset of cells in candArray in the ascending order
      sort(candOffset, candOffset + nCandOffset);
    }
  }
  // Save cell info in candidate space
  for (int j = 0; j < nCandOffset - 1; ++j) {
    if (candOffset[j + 1] - candOffset[j] == 1) {
      if (currSet.cell[candArray[candOffset[j]]] < 0)
        currSet.cellVertex[currSet.nCellVertex++] = candArray[candOffset[j]];
      currSet.cell[candArray[candOffset[j]]] = 0;
      continue;
    }
    Array* a = new Array(currSet.candidates, candArray, candOffset[j],
                         candOffset[j + 1]);
    uIter = cellToID.find(*a);
    int currID;
    if (uIter == cellToID.end()) {
      currID = cellToID[*a] = globalCellID++;
      if (globalCellID >= nMaxCell) {
        int prevSize = CELL_SIZE;
        nMaxCell *= 2;
        CELL_SIZE = ((nMaxCell + 63) >> 6);
        for (int i = 0; i < nQueryVertex; ++i) {
          element[i].buffer = new uint64_t[CELL_SIZE]();
        }
        for (int i = 0; i < d; ++i) {
          for (int j = 0; j < prevSize; ++j) {
            element[i].buffer[j] = element[i].conflictCell[j];
          }
        }
        for (int i = 0; i < nQueryVertex; ++i) {
          element[i].conflictCell = element[i].buffer;
        }
      }
      cellPos[currID] = cellPos[currID - 1];
      for (int k = candOffset[j]; k < candOffset[j + 1]; ++k) {
        cellList[cellPos[currID]++] = currSet.candidates[candArray[k]];
      }
    } else {
      currID = uIter->second;
    }
    for (int k = candOffset[j]; k < candOffset[j + 1]; ++k) {
      if (currSet.cell[candArray[k]] < 0)
        currSet.cellVertex[currSet.nCellVertex++] = candArray[k];
      currSet.cell[candArray[k]] = currID;
    }
  }
}

inline void GetEquivalence(Stack& retE) {
  bool isRoot = (retE.vertex == root);
  CandidateSpace& retSet = candSpace[retE.vertex];
  int initPos = retE.addressPos;
  int idx = isRoot ? initPos : retE.address[initPos];
  int cID = retSet.cell[idx];
  if ((retE.conflictCell[0] & 1ULL) != 0 || cID == 0) return;
  int nCell = 1, cPos = 0, aPos = cellPos[cID - 1];
  int v = retSet.candidates[idx];
  for (int y = cellPos[cID - 1]; y < cellPos[cID]; ++y) {
    arrayToClean[toCleanIndex++] = cellList[y];
    ++isVisited[cellList[y]];
  }
  for (int z = 0; z < CELL_SIZE; ++z) {
    uint64_t arr = retE.conflictCell[z];
    nCell += __builtin_popcountll(arr);
    while (arr) {
      int id = (z << 6) + __builtin_ctzl(arr);
      for (int y = cellPos[id - 1]; y < cellPos[id]; ++y) {
        if (isVisited[cellList[y]] == 1) ++isVisited[cellList[y]];
      }
      arr ^= (arr & -arr);
    }
  }
  for (aPos = initPos; aPos < retE.addressSize; ++aPos) {
    if (retE.isPruned[aPos]) continue;
    idx = isRoot ? aPos : retE.address[aPos];
    v = retSet.candidates[idx];
    if (isVisited[v] == nCell) {
      retE.isPruned[aPos] = true;
      retE.pruned[retE.nPruned++] = aPos;
    }
  }
  while (toCleanIndex != 0) isVisited[arrayToClean[--toCleanIndex]] = 0;
}

inline bool BFS_MM() {
  // BFS function for maximum matching
  /*
   * Note that in this function, we are NOT manipulating the actual query nodes
   * and data nodes. Instead, we are dealing with the index in LLC and the data
   * node index in the sum NEC candidate set.
   */
  int* queue = arrayToClean;
  int queueStart = 0;  // use for popping up
  int queueEnd = 0;    // use for pushing back
  for (int u = 1; u < nodeIndex; u++) {
    if (pairU[u] == 0)  // NULL
    {
      dist[u] = 0;
      queue[queueEnd++] = u;
    } else
      dist[u] = INT_MAX;
  }
  dist[0] = INT_MAX;  // NULL
  while (queueStart != queueEnd) {
    int u = queue[queueStart];
    queueStart++;
    if (dist[u] < dist[0])  // NULL
    {
      int nec = necMapping[u];  // get the actual index in the LLC
      for (int j = 0; j < NECRegionSize[nec]; j++) {
        int v = uCand[nec * maxNumCandidate + j];
        if (dist[pairV[v]] == INT_MAX) {
          dist[pairV[v]] = dist[u] + 1;
          // here use "u" instead of "nec", because "nec" is only used to get
          // the candidate set.
          queue[queueEnd++] = pairV[v];
        }
      }
    }
  }
  return (dist[0] != INT_MAX);  // NULL
}

inline bool DFS_MM(int u) {
  if (u != 0)  // NULL
  {
    int nec = necMapping[u];  // get the actual index in the LLC
    for (int j = 0; j < NECRegionSize[nec]; j++) {
      int v = uCand[nec * maxNumCandidate + j];
      if (dist[pairV[v]] == dist[u] + 1)
        if (DFS_MM(pairV[v])) {
          pairV[v] = u;
          pairU[u] = v;
          return true;
        }
    }
    dist[u] = INT_MAX;
    return false;
  }
  return true;
}

inline void combine(double& cntEmbedding, double resultSoFar) {
  if (cntEmbedding * resultSoFar > 0) return;
  // Step One: find a candidate that is with the max number of NEC nodes'
  // candidate set containing it locate the candidate with the last coverage
  // number
  int maxCand = -1, maxCandSize = 1, maxCandIdx = -1;
  int* localSumNECCand = sumNECCand[labelID];
  int& localSizeSum = nSumNECCand[labelID];
  for (int i = 1; i < localSizeSum; i++) {
    int can = localSumNECCand[i];
    int temp = vCandIndex[can];
    if (mapTo[can] < 0 && temp > maxCandSize) {
      maxCandSize = temp;
      maxCandIdx = i;
    }
  }
  if (maxCandIdx == -1) {
    // in this case, there is no such a cand that has more than one query node
    // having it as a candidate hence, we can perform the combination and
    // permutation to get the result here
    double res = 1;
    for (int i = 0; i < NECCountSetSize; i++) {
      if (NECCountSet[i] != localFlagQuery[i]) {
        // this node is unmatched or still has unmatched NEC nodes
        int rest = NECCountSet[i] - localFlagQuery[i];
        int remainingNECRegionSize = NECRegionSize[i];
        if (remainingNECRegionSize < rest)  // smaller ==> no result
          return;
        if (remainingNECRegionSize == rest)  // equal ==> only one result
          res *= 1;
        else {
          // larger than rest, then perform the combination
          for (int x = remainingNECRegionSize - rest + 1;
               x <= remainingNECRegionSize; x++)
            res *= x;
          res /= factorization(rest);
        }
      }
    }
    res *= PRE_COMPUTED_PERMUTATION;
    cntEmbedding += res;
    return;
  }                                       // end if maxCand = -1
  maxCand = localSumNECCand[maxCandIdx];  // the case that maxCand is not -1
  swapValue(localSumNECCand[localSizeSum - 1],
            localSumNECCand[maxCandIdx]);  // swap to the end
  localSizeSum--;                          // reduce the size by one

  // then we examine the nec nodes that are covered by the max cand
  /*
   * The pruning rule here is that :
   *   for this cand, if there exists multiple query nodes covered by it and
   * those nodes' cand size (the conditional cand size for the nec node) is 1,
   *   then there is no result. Because, one cand cannot be mapped to multiple
   * query node at the same time. But if there only exists one such node, then
   * it is OK to precede and this cand should be mapped to this query node.
   */

  // now we try to map this cand node, and remove it from all candidate regions
  // that contain it this for loop remove the maxCand from the candidate sets
  // that contain it
  for (int i = 0; i < vCandIndex[maxCand]; i++) {
    // for each query node whose cand set containing it
    // qNodeIndex is the position of this maxCand in the nec node set
    int qNodeIndex = vCandInfo[maxCand][i].first;
    int posInRegion = vCandInfo[maxCand][i].second;
    int regionSize = NECRegionSize[qNodeIndex];
    if (posInRegion == regionSize - 1)
      NECRegionSize[qNodeIndex]--;
    else if (posInRegion < regionSize - 1) {
      pair<int, int>& tempP =
          uCandInfo[qNodeIndex * maxNumCandidate + regionSize - 1];
      int lastCand = tempP.first;
      int lastPos = tempP.second;
      swapValue(tempP, uCandInfo[qNodeIndex * maxNumCandidate + posInRegion]);
      vCandInfo[maxCand][i].second = regionSize - 1;
      vCandInfo[lastCand][lastPos].second = posInRegion;
      NECRegionSize[qNodeIndex]--;
    }
  }

  // ====== multiple ======
  // sort accord to size of each query node's effective candidate size
  sort(vCandInfo[maxCand].begin(),
       vCandInfo[maxCand].begin() + vCandIndex[maxCand],
       sortByEffectiveCandSize);
  // the following loop maintains the position info after the sorting
  for (int i = 0; i < vCandIndex[maxCand]; i++) {
    // for each query node containing maxCand
    pair<int, int> tempP = vCandInfo[maxCand][i];
    int qNodeIndex = tempP.first;
    int posInRegion = tempP.second;
    // set the position of "qNodeIndex" in the u_cand set of the candidate
    // "maxCand"
    uCandInfo[qNodeIndex * maxNumCandidate + posInRegion].second = i;
  }
  for (int i = 0; i < vCandIndex[maxCand]; i++) {
    // for each query node containing maxCand
    int qNodeIndex = vCandInfo[maxCand][i].first;
    if ((int)localFlagQuery[qNodeIndex] ==
        NECCountSet[qNodeIndex])  // check whether a nec node has been fully
                                  // mapped or not
      continue;
    mapTo[maxCand] = nQueryVertex;
    localFlagQuery[qNodeIndex]++;  // increase the flag value accordingly
    // here, if the query node "qNodeIndex" has been fully mapped, then we need
    // to delete it from all the data node's candidate set which containing it,
    // such  that this query node will not be counted in the future round of
    // finding the next maxCand
    if (localFlagQuery[qNodeIndex] == NECCountSet[qNodeIndex]) {
      // this query node is fully mapped
      for (int j = 0; j < NECRegionSize[qNodeIndex]; j++) {
        // for each of its candidate, get the candidate from position "j"
        int vNode = uCandInfo[qNodeIndex * maxNumCandidate + j].first;
        vector<pair<int, int>>& vcInfo = vCandInfo[vNode];
        int pos = uCandInfo[qNodeIndex * maxNumCandidate + j].second;
        if (pos == vCandIndex[vNode] - 1)
          vCandIndex[vNode]--;
        else if (pos < vCandIndex[vNode] - 1) {
          int size = vCandIndex[vNode];
          swapValue(uCandInfo[vcInfo[size - 1].first * maxNumCandidate +
                              vcInfo[size - 1].second],
                    uCandInfo[vcInfo[pos].first * maxNumCandidate +
                              vcInfo[pos].second]);
          swapValue(vcInfo[size - 1], vcInfo[pos]);
          vCandIndex[vNode]--;
        }
      }                                  // end for j
    }                                    // end if
    combine(cntEmbedding, resultSoFar);  // recursive function
    if (cntEmbedding * resultSoFar > 0) {
      // clean up and return
      mapTo[maxCand] = -1;
      return;
    }
    // now we have to update the c(v) again, by putting back the last mapped
    // query node (actually here: the qNodeIndex)
    if (localFlagQuery[qNodeIndex] == NECCountSet[qNodeIndex]) {  // NEC equal
      for (int j = 0; j < NECRegionSize[qNodeIndex]; j++) {
        int vNode = uCandInfo[qNodeIndex * maxNumCandidate + j].first;
        if (vNode == maxCand)  // skip "maxCand"
          continue;
        vCandIndex[vNode]++;  // extend the size by one
      }
    }
    localFlagQuery[qNodeIndex]--;  // reduce this query node's flag, because it
                                   // has been unmapped.
  }
  mapTo[maxCand] = -1;
  combine(cntEmbedding, resultSoFar);  // recursive function
  if (cntEmbedding * resultSoFar > 0) {
    // clean up and return
    mapTo[maxCand] = -1;
    return;
  }
  // end =============== (multiple)

  // put back the max cand into each query node's candidate set
  for (int i = 0; i < vCandIndex[maxCand]; i++) {
    int qNodeIndex = vCandInfo[maxCand][i].first;
    NECRegionSize[qNodeIndex]++;
  }
  localSizeSum++;
  mapTo[maxCand] = -1;
}

inline double CountLeafMatch(const Graph& query, double matchFound) {
  double matchResult = 1;
  int NECSetSize = NECElementOffset.size() -
                   1;  // -1 is because the last currE is actually redundant
  int can, uniqueUCandCounter,
      sumNEC;  // the sum number of all nec node of a label
  int uCandCounter, sumCand;
  memset(necMapping, -1, sizeof(int) * (nQueryVertex + 1));  // for MM
  leafCandIndex = 0;
  // === First scan === fastly identify whether all candidates satisfy the
  // (cand < nec) and (sumCand < sum_nec). (the sum cands here are not unique
  // cands)
  for (int i = 0; i < NECSetSize; i++) {
    // for each label
    int c, d, loc;
    int start = NECElementOffset[i].second;    // start position
    int end = NECElementOffset[i + 1].second;  // end position
    sumNEC = 0;
    sumCand = 0;
    for (int j = start; j < end; j++) {
      // basic info about this nec node
      int parentID = query.NECElement[j].parentID;
      int NECNum = query.NECElement[j].sum;
      int uR = query.NECElement[j].representative;
      sumNEC += NECNum;
      uCandCounter = 0;  // record the number of candidate for this nec node
      int parentPos = candPos[parentID];
      CandidateSpace& unit = candSpace[uR];
      order[uR] = query.nNonLeaf;
      GetCell(query, uR, order[uR], unit.adjacent[0][parentPos],
              unit.nAdjacent[0][parentPos]);
      for (int it = 0; it < unit.nAdjacent[0][parentPos]; it++) {
        int idx = unit.adjacent[0][parentPos][it];
        int v = unit.candidates[idx];
        if (mapTo[v] == uR) mapTo[v] = -1;
        if (mapTo[v] < 0) {
          ++uCandCounter;
        } else {
          Stack& ancElem = element[order[mapTo[v]]];
          c = unit.cell[idx];
          ancElem.conflictCell[c >> 6] |= (1ULL << (c & 0x3f));
        }
      }
      if (uCandCounter < NECNum) return 0;
      sumCand += uCandCounter;
    }
    if (sumCand < sumNEC) return 0;
  }
  // === Second scan: extract cands and perform the maximum matching
  for (int i = 0; i < NECSetSize; i++) {       // for each label
    int start = NECElementOffset[i].second;    // start position
    int end = NECElementOffset[i + 1].second;  // end position
    int label = NECElementOffset[i].first;     // the label of this nec set
    // initialization
    posSumNECCand = 1;
    // flagSumNECCand must be all -1 here.
    uniqueUCandCounter = 0;
    sumNEC = 0;
    sumCand = 0;
    nodeIndex = 1;  // for MM
    int* localSumNECCand = sumNECCand[label];
    for (int j = start; j < end; j++) {
      // basic info about this nec node
      int parentID = query.NECElement[j].parentID;
      int NECNum = query.NECElement[j].sum;
      int uR = query.NECElement[j].representative;
      sumNEC += NECNum;
      for (int mm = 0; mm < NECNum; mm++)  // for MM
        necMapping[nodeIndex++] = j;
      uCandCounter = 0;  // record the number of candidate for this nec node
      int parentPos = candPos[parentID];
      CandidateSpace& unit = candSpace[uR];
      int startPos = leafCandIndex;
      for (int it = 0; it < unit.nAdjacent[0][parentPos]; it++) {
        int can = unit.candidates[unit.adjacent[0][parentPos][it]];
        if (mapTo[can] < 0) {
          if (flagSumNECCand[can] ==
              -1) {  // if not already in the candidate set,
            flagSumNECCand[can] =
                posSumNECCand;  // and mark it true and indicate its position in
                                // "sumNECCands"
            localSumNECCand[posSumNECCand++] =
                can;  // then add it into the candidate's "posSumNECCand"
                      // position,
            uniqueUCandCounter++;
          }
          uCand[j * maxNumCandidate + uCandCounter] =
              flagSumNECCand[can];  // for the maximum matching
          uCandCounter++;
          leafCand[leafCandIndex++] =
              can;  // store the candidate for the second stage (combination)
        }
      }
      leafCandInfo[j] = make_pair(startPos, uCandCounter);
      NECRegionSize[j] = uCandCounter;  // set the size of the j-th candidate
                                        // set of this nec node
      sumCand += uCandCounter;
      nSumNECCand[label] = posSumNECCand;  // added on 20:15 09/05/2016
    }
    if (uniqueUCandCounter < sumNEC) {
      // no result
      // reset "flagSumNECCand"
      for (int i = 1; i < posSumNECCand; i++)  // it starts from 1 NOT 0
        flagSumNECCand[localSumNECCand[i]] = -1;
      return 0;
    }
    // after the  "sumNECCands" is initialized, we reset the "flagSumNECCand"
    // for next time use
    for (int i = 1; i < posSumNECCand; i++)  // it starts from 1 NOT 0
      flagSumNECCand[localSumNECCand[i]] = -1;
    // ================= BEGIN checking MAXIMUM MATCHING =================
    // using nec_region as the adjacent list
    memset(pairU, NULL, sizeof(int) * nQueryVertex);  // must reset before using
    memset(pairV, NULL,
           sizeof(int) * currG->nVertex);  // must reset before using
    int matching = 0;
    while (BFS_MM()) {
      for (int u = 1; u < nodeIndex; u++)
        if (pairU[u] == NULL)
          if (DFS_MM(u)) matching++;
    }
    if (matching != nodeIndex - 1) {
      return 0;
    }
    // ================== END checking MAXIMUM MATCHING ==================
    NECRank[i] = make_pair(i, make_pair(sumNEC, sumCand));
  }
  sort(NECRank, NECRank + NECSetSize, sortLLCNode);
  // Deal with the labels one by one
  for (int i = 0; i < NECSetSize; i++) {  // for each nec node label
    if (matchResult > 0) continue;
    double cntLocalMatchLabel = 0;
    int posElem = NECRank[i].first;  // the label index of this nec label set
                                     // ==> this is not the actual label!!!!
    int nodeSum =
        NECRank[i].second.first;  // the number of node in this nec label set
    labelID = NECElementOffset[posElem].first;
    if (nodeSum == 1) {
      // ==== CASE ONE : there is only one node in this nec set with this label
      // we have computed this one in the last step
      cntLocalMatchLabel = (long long)NECRank[i].second.second;
      if (cntLocalMatchLabel == 0) return 0;
    } else {  // ==== CASE TWO : more than one node, and possible more than one
              // start nodes (nec_units)
      memset(localFlagQuery, 0, sizeof(char) * nQueryVertex);
      posSumNECCand = 1;  // omit zero, start from 1 ===> Leave 0 for NULL of
                          // Maximum Matching
      int start = NECElementOffset[posElem].second;
      int end = NECElementOffset[posElem + 1].second;
      int localNECSize = end - start;  // number of nec with posElem
      // ======= sorting here is to decide the processing sequence of the nec
      // nodes in this nec label set
      for (int j = start, x = 0; j < end; j++, x++)
        vNECCount[x] = make_pair(j, query.NECElement[j].sum);
      sort(vNECCount, vNECCount + localNECSize, sortBySecondElement);
      int* toCleanVCandIndex =
          arrayToClean;  // store the candidates that need to be cleaned for the
                         // array "vCandIndex"
      NECCountSetSize = localNECSize;
      int idxToCleanVCandIndex = 0;
      // in this loop, for this point beyond, before any return, "vCandIndex"
      // must be cleaned up using the above two parameters
      for (int j = 0; j < localNECSize; j++) {
        int necIndex = vNECCount[j].first;
        int necCount = vNECCount[j].second;
        NECCountSet[j] = necCount;  // the sum of nodes that the nec
                                    // representative stands for
        uCandCounter = 0;
        pair<int, int> tempP = leafCandInfo[necIndex];
        for (int x = tempP.first; x < tempP.first + tempP.second; x++) {
          can = leafCand[x];
          int temp =
              vCandIndex[can];  // get the position for storing this candidate
          // the temp-th position of vCandInfo[can] stores a candidate "can"
          // which is for "j-th" nec node and stored in the "uCandCounter-th"
          // position in the candidate set
          vCandInfo[can][temp] = make_pair(
              j, uCandCounter);  // pair<j-th nec node, can's pos in the
                                 // candidate set of j-th nec node>
          // store the candidate and its corresponding info ==>  <the candidate
          // can, the candidate's pos in vCandIndex[can]> in the single-array:
          // "j-th" nec node and the "uCandCounter-th" position.
          uCandInfo[j * maxNumCandidate + uCandCounter] = make_pair(can, temp);
          if (temp == 0)  // record for the cleaning purpose
            toCleanVCandIndex[idxToCleanVCandIndex++] = can;
          vCandIndex[can]++;  // update the counter(next available position) for
                              // "vCandIndex[can]"
          uCandCounter++;     // update the candidate counter(next available
                              // position)
        }
        // if the number of candidate for this nec node, is smaller than the nec
        // count of this nec node, then no result can be found.
        NECRegionSize[j] = uCandCounter;
      }  // end for j

      // to reduce the computation cost, we pre-compute the value of the all
      // factorization
      PRE_COMPUTED_PERMUTATION = 1;
      for (int i = 0; i < NECCountSetSize; i++)
        PRE_COMPUTED_PERMUTATION *= factorization(NECCountSet[i]);
      int* mappedDataVertex = flagSumNECCand;
      int posMappedDataVertex = 0;
      // in this loop, for this point beyond, before any return, "mapTo" must be
      // cleaned up using the above two parameters.
      int found = 1;
      while (found > 0) {
        // =============== preprocess before the combination ===============
        found = 0;
        // first, we scan each nec node, to see whether there exists such a node
        // that its candidate size is equal to its unmatched nec size
        for (int i = 0; i < NECCountSetSize; i++) {
          if (NECRegionSize[i] != 0 && NECCountSet[i] > NECRegionSize[i]) {
            while (posMappedDataVertex != 0) {
              mapTo[mappedDataVertex[--posMappedDataVertex]] = -1;
              mappedDataVertex[posMappedDataVertex] = -1;
            }
            while (idxToCleanVCandIndex != 0)  // cleaning up vCandIndex
              vCandIndex[toCleanVCandIndex[--idxToCleanVCandIndex]] = 0;
            return 0;
          }
          // all can be matched one on one
          if (NECCountSet[i] == NECRegionSize[i]) {
            found++;
            localFlagQuery[i] =
                NECCountSet[i];  // mark this nec node all matched
            // update the c(v); there is no neccessary to update the c(v)
            for (int j = 0; j < NECRegionSize[i]; j++) {
              int vNode = uCandInfo[i * maxNumCandidate + j]
                              .first;       // the j-th candidate
              mapTo[vNode] = nQueryVertex;  // set this cand's flag to 1
              mappedDataVertex[posMappedDataVertex++] =
                  vNode;  // recording this cand
              // for each query node that contains this candidate "vNode"
              for (int x = 0; x < vCandIndex[vNode];
                   x++) {  // x is the position in "vCandInfo[vNode]"
                int qNodeIndex =
                    vCandInfo[vNode][x]
                        .first;  // get the original index in the nec node set
                if (qNodeIndex == i) continue;
                // the position of vNode in the corresponding("qNodeIndex-th")
                // candidate set.
                int posInRegion = vCandInfo[vNode][x].second;
                int regionSize = NECRegionSize[qNodeIndex];
                if (posInRegion ==
                    regionSize -
                        1)  // this is the last node in the candidate set
                  NECRegionSize[qNodeIndex]--;  // remove this candidate by
                                                // downsizing the NECRegionSize
                                                // by one
                else if (posInRegion < regionSize - 1) {
                  // normal case: not the last one.
                  // swap the two candidates and their corresponding info, then
                  // downsize the NECRegionSize by one
                  pair<int, int>& tempP =
                      uCandInfo[qNodeIndex * maxNumCandidate + regionSize - 1];
                  int lastCand = tempP.first;
                  int lastPos = tempP.second;
                  swapValue(tempP,
                            uCandInfo[qNodeIndex * maxNumCandidate +
                                      posInRegion]);  // swap the candidate
                  vCandInfo[vNode][x].second = regionSize - 1;
                  vCandInfo[lastCand][lastPos].second = posInRegion;
                  NECRegionSize[qNodeIndex]--;
                }
                if (NECRegionSize[qNodeIndex] < NECCountSet[qNodeIndex]) {
                  while (posMappedDataVertex != 0) {
                    mapTo[mappedDataVertex[--posMappedDataVertex]] = -1;
                    mappedDataVertex[posMappedDataVertex] = -1;
                  }
                  while (idxToCleanVCandIndex != 0)  // cleaning up vCandIndex
                    vCandIndex[toCleanVCandIndex[--idxToCleanVCandIndex]] = 0;
                  return 0;
                }
              }  // end for x
            }    // end for j
            NECRegionSize[i] = 0;
          }  // end if
        }    // end for i
      }      // end PREPROCESSING
      // ======================== END PREPROCESSING ========================
      combine(cntLocalMatchLabel, matchResult);

      // Clean up
      while (posMappedDataVertex != 0) {
        mapTo[mappedDataVertex[--posMappedDataVertex]] = -1;
        mappedDataVertex[posMappedDataVertex] = -1;
      }
      while (idxToCleanVCandIndex != 0)  // cleaning up vCandIndex
        vCandIndex[toCleanVCandIndex[--idxToCleanVCandIndex]] = 0;
      if (cntLocalMatchLabel == 0) return 0;
    }  // end else : end of case TWO
    matchResult *= cntLocalMatchLabel;
  }
  return matchResult;
}

void conflictClass(int uCurr, int uID, int depth) {
  ancestors[depth][uCurr >> 6] |= (1ULL << (uCurr & 0x3f));

  const uint64_t arr = (1ULL << (uID & 0x3f));
  addInFailingSet(uID, currE->failingSet);
  addInFailingSet(uID, ancestors[depth]);
}

inline long long FindProblemLeafMatch(int depth, Stack& elem, int u, int label,
                                      int uP) {
  unordered_map<pair<int, int>, int, pairHash>::iterator sizeIter;
  sizeIter = NECSize.find(make_pair(label, uP));
  int c, d, loc, v, pos, pPos = candPos[uP];
  long long NECNum = sizeIter->second;
  long long uCandCounter = 0;
  CandidateSpace& unit = candSpace[u];
  order[u] = order[uP] + 1;
  GetCell(*currQ, u, order[u], unit.adjacent[0][pPos], unit.nAdjacent[0][pPos]);
  for (int i = 0; i < unit.nAdjacent[0][pPos]; ++i) {
    pos = unit.adjacent[0][pPos][i];
    v = unit.candidates[pos];
    if (mapTo[v] < 0) {
      ++uCandCounter;
    } else {
      Stack& ancElem = element[order[mapTo[v]]];
      c = unit.cell[pos];
      ancElem.conflictCell[c >> 6] |= (1ULL << (c & 0x3f));
    }
  }
  if (uCandCounter < NECNum) {
    for (int i = 0; i < unit.nAdjacent[0][pPos]; ++i) {
      pos = unit.adjacent[0][pPos][i];
      v = unit.candidates[pos];
      if (mapTo[v] < 0) {
        if (useFailingSet) {
          addInFailingSet(u, ancestors[depth]);
        }
      } else {
        if (useFailingSet) {
          int uID = mapTo[v];
          conflictClass(u, uID, depth);
        }
      }
    }
    return 0;
  } else {
    if (uCandCounter == NECNum) {
      elem.problemChild[elem.nProblemChild++] = u;
      for (int i = 0; i < unit.nAdjacent[0][pPos]; ++i) {
        int v = unit.candidates[unit.adjacent[0][pPos][i]];
        if (mapTo[v] < 0) mapTo[v] = u;
      }
      isMapped[u] = true;
    }
    return uCandCounter;
  }
}

inline void ReleaseProblemLeafMatch(Stack& elem, int uP, Stack& childElem) {
  for (int idx = 0; idx < elem.nProblemChild; ++idx) {
    int u = elem.problemChild[idx];
    CandidateSpace& unit = candSpace[u];
    for (int i = 0; i < unit.nAdjacent[0][candPos[uP]]; ++i) {
      int v = unit.candidates[unit.adjacent[0][candPos[uP]][i]];
      if (mapTo[v] == u) mapTo[v] = -1;
    }
    isMapped[u] = false;
  }
  elem.nProblemChild = 0;
}

long long getExtendableCandidates(int u, int uIdx, int up, int candParentIdx,
                                  int numParentsMapped) {
  int c;
  CandidateSpace& unit = candSpace[u];
  const int iecPos = numParentsMapped - 1;
  int uNgbIdx = DAG_ngb_query_ngb_index[up][uIdx];

  iecSize[u][iecPos] = 0;
  if (numParentsMapped == 1) {
    GetCell(*currQ, u, order[up], unit.adjacent[uNgbIdx][candParentIdx],
            unit.nAdjacent[uNgbIdx][candParentIdx]);
    for (int i = 0; i < unit.nAdjacent[uNgbIdx][candParentIdx]; i++) {
      const int candIdx = unit.adjacent[uNgbIdx][candParentIdx][i];
      const int vID = unit.candidates[candIdx];
      iec[u][0][iecSize[u][0]++] = candIdx;
      if (mapTo[vID] < 0) {
        Stack& ancElem = element[order[mapTo[vID]]];
        c = unit.cell[candIdx];
        ancElem.conflictCell[c >> 6] |= (1ULL << (c & 0x3f));
      }
    }
  } else {
    int j = 0, k = 0, candIdx, newCandIdx;
    int indexSize = iecSize[u][iecPos - 1];
    int newIndexSize = unit.nAdjacent[uNgbIdx][candParentIdx];
    GetCell(*currQ, u, order[up], iec[u][iecPos - 1], indexSize);
    while (j < indexSize and k < newIndexSize) {
      candIdx = iec[u][iecPos - 1][j];
      newCandIdx = unit.adjacent[uNgbIdx][candParentIdx][k];
      if (candIdx == newCandIdx) {
        const int vID = unit.candidates[candIdx];
        if (useFailingSet) {
          iec[u][iecPos][iecSize[u][iecPos]++] = candIdx;
          if (mapTo[vID] >= 0) {
            Stack& ancElem = element[order[mapTo[vID]]];
            c = unit.cell[candIdx];
            ancElem.conflictCell[c >> 6] |= (1ULL << (c & 0x3f));
          }
        } else {
          if (mapTo[vID] < 0) {
            iec[u][iecPos][iecSize[u][iecPos]++] = candIdx;
          }
        }
        ++j;
        ++k;
      } else if (candIdx < newCandIdx)
        ++j;
      else
        ++k;
    }
  }
  return (long long)iecSize[u][iecPos];
}

inline void reset(int depth, int ri, int rootCand) {
  while (depth > 0) {
    Stack& tempE = element[depth];
    int uID = tempE.vertex;
    int vID = currMapping[depth];
    isMapped[uID] = false;
    for (int i = 0; i < DAG_ngb_query_size[uID]; ++i) {
      --nMappedParent[DAG_ngb_query[uID][i]];
    }
    queueFllExt.clear_nInserted(depth);
    ReleaseProblemLeafMatch(tempE, uID, element[depth + 1]);
    mapTo[vID] = -1;
    while (tempE.nPruned != 0)
      tempE.isPruned[tempE.pruned[--tempE.nPruned]] = false;
    Stack& higherE = element[depth - 1];
    GetEquivalence(higherE);
    tempE.address = NULL;
    --depth;
  }
  ReleaseProblemLeafMatch(element[0], root, element[1]);
  mapTo[rootCand] = -1;
  queueFllExt.clearQueue();
  if (useFailingSet) isRedundant = false;
  Stack& rootE = element[0];
  GetEquivalence(rootE);
}

void updateExtendableRoot(const Graph& query, int ri, int u, int depth = 1) {
  long long weight;
  for (int j = 0; j < DAG_ngb_query_size[u]; ++j) {
    int rc = DAG_ngb_query[u][j];
    if (query.NECMap[rc] != -1) {
      if (query.isProblemLeaf[rc])
        weight = FindProblemLeafMatch(0, element[0], rc, query.label[rc], root);
    } else {
      weight = getExtendableCandidates(rc, j, u, ri, nMappedParent[rc]);

      if (nMappedParent[rc] == 1) {
        queueFllExt.insertToQueue(rc);
      }
    }
  }  // 3. end for
}

bool updateExtendableU(const Graph& query, int u, int depth, const int uPos,
                       int& uID, int& vCID) {
  bool flag = true;
  long long weight;
  for (int j = 0; j < DAG_ngb_query_size[uCurr]; ++j) {
    int uC = DAG_ngb_query[uCurr][j];

    if (isMapped[uC]) continue;

    if (query.isProblemLeaf[uC])
      weight = FindProblemLeafMatch(depth, *currE, uC, query.label[uC], uCurr);
    else if (query.NECMap[uC] != -1)
      continue;
    else {
      weight = getExtendableCandidates(uC, j, uCurr, uPos, nMappedParent[uC]);
    }
    if (weight == 0) {
      flag = false;
      if (useFailingSet) {
        addInFailingSet(uC, ancestors[depth]);
      }
      break;
    }
    if (query.isProblemLeaf[uC]) continue;
    if (nMappedParent[uC] == 1) {
      queueFllExt.insertToQueue(uC);
      queueFllExt.add_nInserted(depth);
    }
  }
  return flag;
}

void updateReleaseCandidate(int depth) {
  Stack& higherE = element[depth];
  Stack& lowerE = element[depth + 1];
  ReleaseProblemLeafMatch(higherE, higherE.vertex, lowerE);
  mapTo[currMapping[depth]] = -1;
  queueFllExt.removeFromQueue(depth);
  queueFllExt.clear_nInserted(depth);
}

void updateCellInfo(int depth) {
  Stack& lowerE = element[depth];
  while (lowerE.nPruned != 0)
    lowerE.isPruned[lowerE.pruned[--lowerE.nPruned]] = false;

  Stack& higherE = element[depth - 1];
  GetEquivalence(higherE);
}

void updateAncestors(int uCurr) {
  for (int i = 0; i < DAG_ngb_query_size[uCurr]; ++i) {
    int uC = DAG_ngb_query[uCurr][i];
    if (!isMapped[uC]) {
      ++nMappedParent[uC];
      dyanc.addAncestor(uC, uCurr, nMappedParent[uC], nMappedParent[uCurr]);
    }
  }
}

int updateMappingQueryVer(int uCurr, int depth) {
  isMapped[uCurr] = true;
  order[uCurr] = depth;
  updateAncestors(uCurr);
  if (useFailingSet) {
    for (int x = 0; x < FAILING_SET_SIZE; ++x) {
      currE->failingSet[x] = 0;
    }
  }
  currE->address = iec[uCurr][nMappedParent[uCurr] - 1];
  currE->addressSize = iecSize[uCurr][nMappedParent[uCurr] - 1];
}

void updateReleaseQueryVer(int uCurr, int depth) {
  queueFllExt.reinsertToQueue(uCurr, depth);
  isMapped[uCurr] = false;
  for (int i = 0; i < DAG_ngb_query_size[uCurr]; ++i) {
    int uC = DAG_ngb_query[uCurr][i];
    if (!isMapped[uC]) {
      --nMappedParent[uC];
    }
  }
}

int getNextuPosIdx() {
  weight_type currVal, maxVal;
  int optOffset;
  int uPosIdx;
  uPosIdx = currE->addressPos;
  return uPosIdx;
}

void updateFailingSet(int depth, bool isRedundant) {
  for (int x = 0; x < FAILING_SET_SIZE; ++x) {
    if (!isRedundant) {
      uint64_t arr = ancestors[depth][x];
      while (arr) {
        int idx = (x << 6) + __builtin_ctzl(arr);
        for (int y = 0; y < FAILING_SET_SIZE; ++y) {
          currE->failingSet[y] |=
              dyanc.getSetPartition(idx, nMappedParent[idx], y);
        }
        arr ^= (arr & -arr);
      }
    }
    ancestors[depth][x] = 0;
  }
}

void moveUpFailingSet(const Stack& higherE, int depth) {
  int uID = higherE.vertex;
  b = currE->failingSet;
  if (depth != 0) {
    if ((b[uID >> 6] & (1ULL << (uID & 0x3f))) == 0) {
      for (int x = 0; x < FAILING_SET_SIZE; ++x) higherE.failingSet[x] = b[x];
      isRedundant = true;
    } else {
      for (int x = 0; x < FAILING_SET_SIZE; ++x) higherE.failingSet[x] |= b[x];
    }
  }
}

inline bool findAllEmbeddings(int dataGraphID) {
  if (answer[dataGraphID]) return true;
  return false;
}

inline bool cannotMap(int dataGraphID) {
  if (currE->addressPos == currE->addressSize || isRedundant ||
      answer[dataGraphID]) {
    return true;
  }
  return false;
}

inline void Backtrack(const Graph& query, const Graph& data, int dataGraphID) {
  long long weight;
  int uPos, vID, vCID, uID;
  bool flag;
  CandidateSpace *currSet, *childSet;

  dyanc.clear(root);
  queueFllExt.clearQueue();

  memset(isMapped, false, sizeof(bool) * nQueryVertex);
  memset(nMappedParent, 0, sizeof(int) * nQueryVertex);

  isRedundant = false;
  if (useFailingSet) {
    for (int i = 0; i < nQueryVertex; ++i)
      memset(ancestors[i], 0, sizeof(uint64_t) * FAILING_SET_SIZE);
  }
  Stack& rootE = element[0];
  order[root] = 0;
  while (rootE.nPruned != 0)
    rootE.isPruned[rootE.pruned[--rootE.nPruned]] = false;
  for (int i = 0; i < uSequenceSize; ++i) {
    CandidateSpace& currSet = candSpace[uSequence[i]];
    while (currSet.nCellVertex != 0)
      currSet.cell[currSet.cellVertex[--currSet.nCellVertex]] = -1;
  }
  isMapped[root] = true;

  updateAncestors(root);

  rootE.vertex = root;
  rootE.address = candSpace[root].candidates;
  rootE.addressPos = -1;
  rootE.addressSize = candSpace[root].size;
  GetCell(query, root, 0, NULL, -1);
  recursiveCallCountPerGraph = 1;
  ++recursiveCallCount;
  for (int ri = 0; ri < candSpace[root].size; ++ri) {
    if (findAllEmbeddings(dataGraphID)) {
      break;
    }
    // 1. For each data vertex v in r.C
    ++rootE.addressPos;
    int rootCand = candSpace[root].candidates[ri];
    if (rootE.isPruned[ri]) {
      continue;
    } else {
      for (int x = 0; x < CELL_SIZE; ++x) rootE.conflictCell[x] = 0;
    }
    mapTo[rootCand] = root;
    currMapping[0] = rootCand;  // M(r) = v
    candPos[root] = ri;
    char backtrack = 0;
    // must reset this to 1 here
    int depth = 1;  // the current query index of the query sequence, because 0
                    // is the root has already been matched
    // 3. For each child uc of r in dag(q)
    updateExtendableRoot(query, ri, root);
    while (true) {
      if (depth == 0) break;  // "No MATCH found!"
      if (depth == query.nNonLeaf) {
        if (query.numNECMapping != 0) {
          nCurrMatch = CountLeafMatch(query, nMatch);
          if (nCurrMatch > 0 && !answer[dataGraphID]) {
            ++nMatch;
            answer[dataGraphID] = true;
          }
        } else {
          nCurrMatch = 1;
          if (!answer[dataGraphID]) ++nMatch;
          answer[dataGraphID] = true;
        }
        if (findAllEmbeddings(dataGraphID)) {
          element[depth].address = NULL;
          --depth;
          reset(depth, ri, rootCand);
          return;
        }

        updateCellInfo(depth);
        --depth;
        updateReleaseCandidate(depth);
        if (useFailingSet) {
          for (int x = 0; x < FAILING_SET_SIZE; ++x)
            element[depth].failingSet[x] = ~0ULL;
        }
        continue;
      }  // if (depth == query.nNonLeaf)
      currE = &element[depth];
      // Find \intersection_{u_p in u.p} (N^{u_p}_{u}(M(u_p)))
      if (currE->address == NULL) {
        ++recursiveCallCount;
        ++recursiveCallCountPerGraph;
        queueFllExt.popFromQueue(uCurr, depth);

        currE->vertex = uCurr;
        currSet = &candSpace[uCurr];

        updateMappingQueryVer(uCurr, depth);
        if (currE->addressSize == 0) {
          currE->address = NULL;
          isRedundant = false;
          if (depth != 0) {
            if (useFailingSet) {
              updateFailingSet(depth, isRedundant);
            }
            updateReleaseQueryVer(uCurr, depth);
          }         // if (depth != 0)
          --depth;  // roll back one node in the matching sequence
          // updateReleaseCandidate(depth);
          continue;
        }  // if (currE->addressSize == 0 || isRedundant)
        currE->addressPos = 0;
      } else  // currE->address != NULL
      {
        uCurr = currE->vertex;
        currSet = &candSpace[uCurr];
        currE->addressPos++;  // update the index by one

        if (cannotMap(dataGraphID)) {
          if (useFailingSet) {
            updateFailingSet(depth, isRedundant);
            isRedundant = false;
          }
          updateReleaseQueryVer(uCurr, depth);
          updateCellInfo(depth);
          Stack& higherE = element[depth - 1];
          currE->address = NULL;
          --depth;
          updateReleaseCandidate(depth);

          uID = higherE.vertex;
          if (useFailingSet) {
            moveUpFailingSet(higherE, depth);
          }
          continue;
        }  // if (currE->addressPos == currE->addressSize || isRedundant)
      }    // if (currE->address == NULL)
      backtrack = 0;  // actually, this line is not necessary, when processed
                      // here, the backtrack must be false...
      while (true) {
        // Break, until find a mapping for this node
        // or cannot find a mapping after examining all candidates
        // no recursive call count increase here!
        if (answer[dataGraphID]) {
          while (depth != 0) {
            mapTo[currMapping[depth]] = -1;
            if (useFailingSet) {
              for (int x = 0; x < FAILING_SET_SIZE; ++x)
                ancestors[depth][x] = 0;
            }
            element[depth].address = NULL;  // 20170414
            --depth;
          }
          mapTo[rootCand] = -1;
          break;
        }

        int uPosIdx = getNextuPosIdx();
        uPos = currE->address[uPosIdx];
        vID = currSet->candidates[uPos];

        if (mapTo[vID] < 0)  // if v is unvisited
        {
          int cellID = currSet->cell[uPos];
          if (currE->isPruned[uPosIdx]) {
            currE->addressPos++;
            if (currE->addressPos == currE->addressSize) {
              backtrack = 1;
              break;
            }
            continue;
          }
          currMapping[depth] = vID;  // 3. M(u) = v
          candPos[uCurr] = uPos;
          mapTo[vID] = uCurr;  // 3. Mark v as visited
          if (usePathSizeOrder) {
            queueFllExt.set_minPosition(-1);
            queueFllExt.set_optWeight(WEIGHT_MAX);
          }
          for (int x = 0; x < CELL_SIZE; ++x) {
            currE->conflictCell[x] = 0;
          }
          // 4. For each child uc of u in dag(q) s.t.
          // all parents of uc has been mapped
          if (!updateExtendableU(query, uCurr, depth, uPos, uID, vCID)) {
            // 6. Remove just inserted vertices from Q
            updateReleaseCandidate(depth);
            currE->addressPos++;
            if (currE->addressPos == currE->addressSize) {
              backtrack = 1;
              break;
            } else {
              continue;
            }
          }
          break;
        } else {  // If vID is already visited
          uID = mapTo[vID];
          int cellID = currSet->cell[uPos];
          Stack& ancElem = element[order[uID]];
          ancElem.conflictCell[cellID >> 6] |= (1ULL << (cellID & 0x3f));
          if (useFailingSet) {
            conflictClass(uCurr, uID, depth);
          }
          currE->addressPos++;  // not ok, then we need the next result
          if (currE->addressPos ==
              currE->addressSize) {  // No more extendable vertex, so cannot
                                     // find a match for this query node
            backtrack = 1;
            break;
          }
        }
      }  // end while
      if (answer[dataGraphID]) break;
      if (backtrack) {
        backtrack = 0;
        // 7. Reinsert (u, weight) to Q.
        if (useFailingSet) {
          updateFailingSet(depth, isRedundant);
          isRedundant = false;
        }
        updateReleaseQueryVer(uCurr, depth);
        updateCellInfo(depth);
        currE->address = NULL;
        --depth;

        Stack& higherE = element[depth];
        updateReleaseCandidate(depth);
        uID = higherE.vertex;
        if (useFailingSet) {
          moveUpFailingSet(higherE, depth);
        }
      } else  // if(!backtrack)
      {
        depth++;
      }
    }  // while
    reset(depth, ri, rootCand);
  }  // for candSpace[root]
}
