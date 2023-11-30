#include <algorithm>
#include <unordered_map>

#include "config.h"
using namespace std;

int* cntLabel = nullptr;
int* positiveLabel = nullptr;

struct pairHash {
  size_t operator()(const pair<int, int>& p) const {
    return p.first ^ p.second;
  }
};

struct Element {
  int label = 0;
  int parentID = 0;
  int representative = 0;
  int sum = 0;
  Element() {}
  Element(int label, int parentID, int representative) {
    this->label = label;
    this->parentID = parentID;
    this->representative = representative;
  }
};

#ifdef EDGE_LABEL
class PairArray {
 public:
  pair<int, int>* arr = nullptr;
  PairArray() {}
  PairArray(int size) { arr = new pair<int, int>[size]; }
  inline int operator[](const int idx) const { return arr[idx].first; }
  inline pair<int, int>& get(const int idx) { return arr[idx]; }
  inline void set(int idx, int fst, int snd) { arr[idx] = make_pair(fst, snd); }
  inline void set(int idx, int fst) { arr[idx].first = fst; }
};
#endif
class Graph {
 public:
  number_type nVertex = 0;
  number_type nEdge = 0;
  int nLabel = 0;
  int maxDegree = 0;
  int* label = 0;
#ifdef EDGE_LABEL
  PairArray nbr;
#else
  int* nbr = 0;
#endif
  number_type* nbrOffset = nullptr;
  int* degree = nullptr;
  // int* vSortedByLabelDegree;
  // int* degreeArray;
  // pair<int, int>* posLabel;
  int* vertex = nullptr;
  int* vertexOffset = nullptr;
  int* maxNbrDegree = nullptr;
#ifdef NEIGHBOR_LABEL_FREQUENCY
  unordered_map<int, int>* NLF = nullptr;
#else
  // int* NLF;
  uint64_t* NLF = nullptr;
#endif
  unordered_map<int, int> labelFrequency;
#ifdef HUGE_GRAPH
#else
  unordered_map<pair<int, int>, pair<number_type, number_type>, pairHash>
      labelToNbrOffset;
#endif
  // int* labelToNbrOffset;

  // Only for a data graph
  int mostFrequentLabelID = 0;
  int maxNumSameLabelVertex = 0;

  // Only for a query graph
  int nNonLeaf = 0;
  int* core = nullptr;
  int* NECMapping = nullptr;
  Element* NECElement = nullptr;
  int* NECMap = nullptr;
  int numNECMapping = 0;
#ifdef LEAF_ADAPTIVE_MATCHING
  bool* isProblemLeaf = nullptr;
#endif
  // Only for tests
  bool fail = 0;

  Graph() {
    nVertex = 0;
    nEdge = 0;
    nLabel = 0;
    maxDegree = 0;

    label = NULL;
#ifdef EDGE_LABEL
#else
    nbr = NULL;
#endif
    nbrOffset = NULL;
    core = NULL;
    degree = NULL;
    vertex = NULL;        // vertices sorted by labels
    vertexOffset = NULL;  // vertexOffset[l]: the offset of the starting
                          // position of vertices with label l
    maxNbrDegree = NULL;
#ifdef NEIGHBOR_LABEL_FREQUENCY
#else
    NLF = NULL;
#endif
    NECMap = NULL;
    numNECMapping = 0;
    NECMapping = NULL;
    NECElement = NULL;
  }

  inline void initVertex(const number_type n) {
    nVertex = n;
    label = new int[nVertex]();
    nbrOffset = new number_type[nVertex + 1]();
    core = new int[nVertex]();
    degree = new int[nVertex]();
    vertex = new int[nVertex]();
    maxNbrDegree = new int[nVertex]();
  }

  inline void initEdge(const number_type n) {
    nEdge = n;
#ifdef EDGE_LABEL
    nbr = PairArray(nEdge * 2);
#else
    nbr = new int[nEdge * 2]();
#endif
  }

  inline void sortNeighbors() {
    int prev, curr;
    number_type currCnt, currIndex;
    // Sort adjaceny lists in the ascending order of labels
    // For the same labels, sort in the ascending order of degrees
    vertexOffset = new int[nUniqueLabel + 1]();

#ifdef HUGE_GRAPH
    for (int i = 0; i < nVertex; ++i) ++vertexOffset[label[i]];
#else
    for (int i = 0; i < nVertex; ++i) {
      ++vertexOffset[label[i]];
      if (degree[i] == 0) continue;
      sort(
#ifdef EDGE_LABEL
          nbr.arr + nbrOffset[i], nbr.arr + nbrOffset[i + 1],
          [this](const pair<int, int>& p1, const pair<int, int>& p2) -> bool {
            if (label[p1.first] == label[p2.first])
              return (degree[p1.first] < degree[p2.first]);
            else
              return (label[p1.first] < label[p2.first]);
          }
#else
          nbr + nbrOffset[i], nbr + nbrOffset[i + 1],
          [this](const int& v1, const int& v2) -> bool {
            if (label[v1] == label[v2])
              return (degree[v1] < degree[v2]);
            else
              return (label[v1] < label[v2]);
          }
#endif
      );
      currCnt = 1;
      currIndex = nbrOffset[i];
      prev = label[nbr[currIndex]];
      if (degree[i] > 1) {
        number_type j = nbrOffset[i];
        for (j = nbrOffset[i] + 1; j < nbrOffset[i + 1]; ++j) {
          curr = label[nbr[j]];
          if (prev == curr) {
            ++currCnt;
            continue;
          } else {
            labelToNbrOffset[make_pair(i, prev)] =
                make_pair(currIndex, currCnt);
            currIndex += currCnt;
            prev = curr;
            currCnt = 1;
          }
        }
        labelToNbrOffset[make_pair(i, prev)] = make_pair(currIndex, currCnt);
        continue;
      } else {
        labelToNbrOffset[make_pair(i, prev)] = make_pair(currIndex, currCnt);
      }
    }
#endif

    // Calculate label-to-vertex lists
    curr = vertexOffset[0];
    vertexOffset[0] = 0;
    for (int i = 1; i <= nUniqueLabel; ++i) {
      prev = vertexOffset[i];
      vertexOffset[i] = vertexOffset[i - 1] + curr;
      curr = prev;
    }
    for (int i = 0; i < nVertex; ++i) {
      vertex[vertexOffset[label[i]]++] = i;
    }
    for (int i = nUniqueLabel; i >= 1; --i) {
      vertexOffset[i] = vertexOffset[i - 1];
    }
    vertexOffset[0] = 0;
    for (int i = 0; i < nUniqueLabel; ++i) {
      sort(vertex + vertexOffset[i], vertex + vertexOffset[i + 1],
           [this](const int& v1, const int& v2) -> bool {
             return (degree[v1] < degree[v2]);
           });
    }
  }

  inline void computeNLF() {
    int lid, nPositive = 0;
#ifdef NEIGHBOR_LABEL_FREQUENCY
    unordered_map<int, int>::iterator iter;
    NLF = new unordered_map<int, int>[nVertex];
#else
    NLF = new uint64_t[(long long)nVertex * (long long)NLFSize]();
#endif
    for (int j = 0; j < nVertex; ++j) {
      if (degree[j] == 0) continue;
      for (number_type k = nbrOffset[j]; k < nbrOffset[j + 1]; ++k) {
        int neighbor = nbr[k];
        if (degree[neighbor] > maxNbrDegree[j])
          maxNbrDegree[j] = degree[neighbor];
        lid = label[neighbor];
#ifdef NEIGHBOR_LABEL_FREQUENCY
        iter = NLF[j].find(lid);
        if (iter == NLF[j].end()) {
          NLF[j][lid] = 1;
        } else {
          iter->second += 1;
        }
#else
        if (cntLabel[lid] == 0) positiveLabel[nPositive++] = lid;
        if (cntLabel[lid] < BITS_PER_LABEL) {
          int start = lid * BITS_PER_LABEL;
          // int idx = NLFSize - 1 - (start + cntLabel[lid]) / INT_SIZE;
          // int pos = (start + cntLabel[lid]) % INT_SIZE;
          int idx = NLFSize - 1 - (start + cntLabel[lid]) / UINT64_SIZE;
          int pos = (start + cntLabel[lid]) % UINT64_SIZE;
#ifdef HUGE_GRAPH
          NLF[(long long)j * (long long)NLFSize + (long long)idx] |=
              (1ULL << pos);
#else
          NLF[j * NLFSize + idx] |= (1ULL << pos);
#endif
          ++cntLabel[lid];
        }
#endif
      }
#ifdef NEIGHBOR_LABEL_FREQUENCY
#else
      while (nPositive > 0) cntLabel[positiveLabel[--nPositive]] = 0;
#endif
    }
  }
#ifdef PRINT_LOG
  void print(int id = 0) {
    unordered_map<pair<int, int>, pair<int, int>, pairHash>::const_iterator
        iter;
    cout << "[Graph " << id << "]" << endl;
    cout << "|V|: " << nVertex << ". |E|: " << nEdge << ". |L|: " << nLabel
         << ". maxD: " << maxDegree << endl;
    for (int i = 0; i < nVertex; ++i) {
      cout << "(u" << i << ") lv: " << label[i] << ". degree: " << degree[i]
           << ". nbr:";
      for (number_type j = nbrOffset[i]; j < nbrOffset[i + 1]; ++j) {
#ifdef EDGE_LABEL
        pair<int, int> p = nbr.get(j);
        cout << " u" << p.first << "(le: " << p.second << ")";
#else
        cout << " u" << nbr[j];
#endif
      }
      cout << endl;
      for (int j = 0; j < nUniqueLabel; ++j) {
        iter = labelToNbrOffset.find(make_pair(i, j));
        if (iter == labelToNbrOffset.end()) continue;
        cout << "\tN(u: u" << i << ", lv: " << j << "):";
        for (int k = iter->second.first;
             k < iter->second.first + iter->second.second; ++k) {
          cout << " u" << nbr[k];
        }
        cout << endl;
      }
    }
  }
#endif
};

// Variables for buildling CS
struct CandidateSpace {
  int size = 0;         // the size for both path and candidates
  int* candidates = nullptr;  // candidate set
  int*** adjacent =
      NULL;  // adjacent[i][j] = candidates of this unit when the i-th parent,
             // regarding DAG, mapped to the j-th candidate of the parent.
  int** nAdjacent =
      NULL;  // nAdjacent[i][j] = size of back_trak_index[i][j]. That is, the
             // number of candidates of this unit when the i-th parent mapped to
             // the j-th candidate of the parent.

#ifdef N_OPTIMIZATION
  int** capacity = NULL;
  int** capacityNgb = NULL;
#else
  int* nParentCand =
      NULL;  // nParentCand[i] = the number of candidates of the i-th parent.
  int* nNgbCand = NULL;
#endif
  long long* weight = NULL;
#ifdef PRUNING_BY_EQUIVALENCE_SETS
  int* cell = nullptr;
  int* cellVertex = nullptr;
  int nCellVertex = 0;
#endif
};
