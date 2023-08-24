#include <algorithm>

#include "config.h"
using namespace std;

int* cntLabel;
int* positiveLabel;

struct pairHash {
  size_t operator()(const pair<int, int>& p) const {
    return p.first ^ p.second;
  }
};

struct Element {
  int label;
  int parentID;
  int representative;
  int sum;
  Element() {}
  Element(int label, int parentID, int representative) {
    this->label = label;
    this->parentID = parentID;
    this->representative = representative;
  }
};

class Graph {
 public:
  number_type nVertex;
  number_type nEdge;
  int nLabel;
  int maxDegree;
  int* label;
  int* nbr;
  number_type* nbrOffset;
  int* degree;
  int* vertex;
  int* vertexOffset;
  int* maxNbrDegree;
  uint64_t* NLF;
  unordered_map<int, int> labelFrequency;
  unordered_map<pair<int, int>, pair<number_type, number_type>, pairHash>
      labelToNbrOffset;

  // Only for a data graph
  int mostFrequentLabelID;
  int maxNumSameLabelVertex;

  // Only for a query graph
  int nNonLeaf;
  int* core;
  int* NECMapping;
  Element* NECElement;
  int* NECMap;
  int numNECMapping;
  bool* isProblemLeaf;

  // Only for tests
  bool fail;

  Graph() {
    nVertex = 0;
    nEdge = 0;
    nLabel = 0;
    maxDegree = 0;

    label = NULL;
    nbr = NULL;
    nbrOffset = NULL;
    core = NULL;
    degree = NULL;
    vertex = NULL;        // vertices sorted by labels
    vertexOffset = NULL;  // vertexOffset[l]: the offset of the starting
                          // position of vertices with label l
    maxNbrDegree = NULL;
    NLF = NULL;
    NECMap = NULL;
    numNECMapping = 0;
    NECMapping = NULL;
    NECElement = NULL;
  }

  inline void initVertex(const number_type n) {
    nVertex = n;
    label = new int[nVertex];
    nbrOffset = new number_type[nVertex + 1];
    core = new int[nVertex];
    degree = new int[nVertex]();
    // vSortedByLabelDegree = new int[nVertex];
    // degreeArray = new int[nVertex];
    // posLabel = new pair<int, int>[nVertex];
    vertex = new int[nVertex];
    maxNbrDegree = new int[nVertex]();
  }

  inline void initEdge(const number_type n) {
    nEdge = n;
    nbr = new int[nEdge * 2];
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
      sort(nbr + nbrOffset[i], nbr + nbrOffset[i + 1],
           [this](const int& v1, const int& v2) -> bool {
             if (label[v1] == label[v2])
               return (degree[v1] < degree[v2]);
             else
               return (label[v1] < label[v2]);
           });
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
    NLF = new uint64_t[(long long)nVertex * (long long)NLFSize]();
    for (int j = 0; j < nVertex; ++j) {
      if (degree[j] == 0) continue;
      for (number_type k = nbrOffset[j]; k < nbrOffset[j + 1]; ++k) {
        int neighbor = nbr[k];
        if (degree[neighbor] > maxNbrDegree[j])
          maxNbrDegree[j] = degree[neighbor];
        lid = label[neighbor];
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
      }
      while (nPositive > 0) cntLabel[positiveLabel[--nPositive]] = 0;
    }
  }
};

// Variables for buildling CS
struct CandidateSpace {
  int size;         // the size for both path and candidates
  int* candidates;  // candidate set
  int*** adjacent =
      NULL;  // adjacent[i][j] = candidates of this unit when the i-th parent,
             // regarding DAG, mapped to the j-th candidate of the parent.
  int** nAdjacent =
      NULL;  // nAdjacent[i][j] = size of back_trak_index[i][j]. That is, the
             // number of candidates of this unit when the i-th parent mapped to
             // the j-th candidate of the parent.

  int** capacity = NULL;
  int** capacityNgb = NULL;
  long long* weight = NULL;
  int* cell;
  int* cellVertex;
  int nCellVertex;
};
