#include <limits.h>  // for LLONG_MAX which is 9223372036854775807

#include "process.h"
using namespace std;

// Variables for search
int uCurr;
int labelID = -1;
pair<int, pair<int, double>> NECRank[MAX_NUM_VERTEX];
int NECRegionSize[MAX_NUM_VERTEX];
int NECCountSet[MAX_NUM_VERTEX];
int NECCountSetSize;
pair<int, int>
    vNECCount[MAX_NUM_VERTEX];  // used in the map special tree function
double PRE_COMPUTED_PERMUTATION;
int candPos[MAX_NUM_VERTEX];
int currMapping[MAX_NUM_VERTEX];
int nMappedParent[MAX_NUM_VERTEX];
#ifdef QUERY_VERTEX_WITH_MAX
weight_type WEIGHT_MIN = LLONG_MIN;
#else
weight_type WEIGHT_MAX = LLONG_MAX;
#endif
int* iec[MAX_NUM_VERTEX][MAX_QUERY_DEGREE];
int iecSize[MAX_NUM_VERTEX][MAX_QUERY_DEGREE];

// Variables for pruning by failing sets
int FAILING_SET_SIZE = ((MAX_NUM_VERTEX + 63) >> 6);
uint64_t* ancestors[MAX_NUM_VERTEX];
bool isRedundant;
// Variables for pruning by equivalence sets
#ifdef PRUNING_BY_EQUIVALENCE_SETS
int CELL_SIZE;
int nMaxCell;
#endif

// Variables for filtering by neighbor safety
#ifdef FILTERING_BY_NEIGHBOR_SAFETY
using lint = long long;
bool on = false;
int index_vis_color = 0;
int index_set_color = 0;
int index_vis_adj = 0;
int index_ngb_existence = 0;
double s1 = 0;
double s2 = 0;
#ifdef HUGE_GRAPH
int DATAV;
int LABEL;
int** color;              // label 2
int* set_color;           // label
int* vis_adj;             // datav
int** num_ngb_existence;  // datav 2
lint* cnt_included_cs;    // datav
#else
const int DATAV = 4300000;
const int LABEL = 50000;
int color[LABEL][2] = {
    0,
};
int set_color[LABEL] = {
    0,
};
int vis_adj[DATAV] = {
    0,
};
int num_ngb_existence[DATAV][2] = {
    0,
};
lint cnt_included_cs[DATAV] = {
    0,
};
#endif
#endif

inline long long getWeightExtendable(int u) {
  return iecSize[u][nMappedParent[u] - 1];
}

inline bool isInFailingSet(int u, const uint64_t* failingSet) {
  return failingSet[u >> 6] & (1ULL << (u & 0x3f));
}

inline void addInFailingSet(int u, uint64_t* failingSet) {
  failingSet[u >> 6] |= (1ULL << (u & 0x3f));
}

inline void swapValue(int& v1, int& v2) {
  int temp = v1;
  v1 = v2;
  v2 = temp;
}

inline void swapValue(pair<int, int>& v1, pair<int, int>& v2) {
  pair<int, int> temp = v1;
  v1 = v2;
  v2 = temp;
}

inline bool sortLLCNode(pair<int, pair<int, double>> a,
                        pair<int, pair<int, double>>
                            b) {  // the second pair is <node sum, cand sum>
  if (a.second.first == b.second.first)
    return a.second.second < b.second.second;
  else
    return a.second.first < b.second.first;
}
inline bool sortByEffectiveCandSize(pair<int, int> p1, pair<int, int> p2) {
  return NECRegionSize[p1.first] < NECRegionSize[p2.first];
}
inline bool sortBySecondElement(pair<int, int> p1, pair<int, int> p2) {
  return (p1.second < p2.second);
}
inline bool sortByNECLabel(const Element& p1, const Element& p2) {
  return (p1.label < p2.label);
}
inline bool sortQueryByDegree(const int id1, const int id2) {
  return (currQ->degree[id1] > currQ->degree[id2]);
}
inline bool sortQueryByLabelFrequency(const int id1, const int id2) {
  return (currG->labelFrequency[currQ->label[id1]] <
          currG->labelFrequency[currQ->label[id2]]);
}
inline bool sortPairsByWeight(const pair<long long, long long>& p1,
                              const pair<long long, long long>& p2) {
  return p1.second < p2.second;
}

#ifdef SUBGRAPH_MATCHING
inline string doubleToString(double value) {
  string strVal = to_string(value);
  int pos = strVal.find_first_of(".");
  return strVal.substr(0, pos);
}
#endif

inline double factorization(int x) {
  switch (x) {
    case 1:
      return 1;
    case 2:
      return 2;
    case 3:
      return 6;
    case 4:
      return 24;
    case 5:
      return 120;
    case 6:
      return 720;
    case 7:
      return 5040;
    case 8:
      return 40320;
    case 9:
      return 362880;
    case 10:
      return 3628800;
    case 11:
      return 39916800;
    case 12:
      return 479001600;
    case 13:
      return 6227020800;
    case 14:
      return 87178291200;
    case 15:
      return 1307674368000;
    case 16:
      return 20922789888000;
    case 17:
      return 355687428096000;
    case 18:
      return 6402373705728000;
    case 19:
      return 121645100408832000;
    case 20:
      return 2432902008176640000;
  }

  double result = 2432902008176640000;
  for (double i = 21; i <= x; i++) result *= i;

  return result;
}
