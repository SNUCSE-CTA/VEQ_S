#include "util.h"

#ifdef FILTERING_BY_NEIGHBOR_SAFETY
enum direction { topDown, bottomUp };
enum ngbType { parentNgb, childNgb, allNgb };
#endif

#ifdef PRUNING_BY_EQUIVALENCE_SETS
class Array {
 public:
  int* b;
  int size;
  Array(){};
  Array(int* arr, int* pos, int start, int end) {
    size = end - start;
    b = new int[size];
    for (int i = 0; i < size; ++i) b[i] = arr[pos[start + i]];
    sort(b, b + size);
  }

  bool operator==(const Array& other) const {
    if (size != other.size)
      return false;
    else {
      bool ret = true;
      for (int x = 0; x < size; ++x) ret &= (b[x] == other.b[x]);
      return ret;
    }
  }
};

struct ArrayHash {
  size_t operator()(const Array& array) const {
    size_t ret = 0;
    for (int i = 0; i < array.size; ++i) ret = ret ^ array.b[i];
    return ret + array.b[0];
  }
};
#endif

struct Stack {
  int* address = NULL;
  int addressSize;
  int addressPos;
  int vertex;
  // for failing set
  uint64_t* failingSet;
#ifdef PRUNING_BY_EQUIVALENCE_SETS
  uint64_t* buffer;
  uint64_t* conflictCell;
  bool* isPruned;
  int nPruned;
  int* pruned;
#ifdef SUBGRAPH_MATCHING
  uint64_t* auxBuffer;
  uint64_t* nonConflictCell;
  bool* neverVisited;
  double* nFoundMatch;
#endif
#endif
#ifdef LEAF_ADAPTIVE_MATCHING
  int* problemChild;
  int nProblemChild;
#endif
#ifdef CANDIDATE_VERTEX_WITH_MAX
  bool* isVisited;
#endif
};

struct Queue {
 private:
  int posExtendable = 0;
  int minPosition = -1;
  long long optWeight = LLONG_MAX;
  int extendable[MAX_NUM_VERTEX] = {
      0,
  };
  int positions[MAX_NUM_VERTEX] = {
      0,
  };
  int nInserted[MAX_NUM_VERTEX] = {
      0,
  };

 public:
  Queue() { this->posExtendable = 0; }
#ifdef MAPPING_FUNCTION_LOG
  inline void printQueue() const {
    string str = "Q:";
    for (int i = 0; i < posExtendable; i++) {
      str += " u" + to_string(extendable[i]);
    }
    printf("%s\n", str.c_str());
  }
#endif
  inline void reinsertToQueue(int u, int depth) {
    int position = positions[depth];

    extendable[posExtendable] = extendable[position];
    posExtendable++;

    extendable[position] = u;
#ifdef MAPPING_FUNCTION_LOG
    printf("\t(rI) u%d (pos: %d), Min. weight: %ld (|Q|: %d)\n", u, minPosition,
           optWeight, posExtendable);
#endif
  }

  inline void insertToQueue(int u) {
    extendable[posExtendable] = u;
    posExtendable++;
#ifdef MAPPING_FUNCTION_LOG
    printf("\t(I) u%d (pos: %d), Min. weight: %ld (|Q|: %d)\n", u, minPosition,
           optWeight, posExtendable);
#endif
  }

  inline void popFromQueue(int& current, int depth) {
    int optPos = -1;

    weight_type optWeight = WEIGHT_MAX;

    for (int i = 0; i < posExtendable; i++) {
      if (getWeightExtendable(extendable[i]) < optWeight) {
        optPos = i;
        optWeight = getWeightExtendable(extendable[i]);
      }
    }

    current = extendable[optPos];
    positions[depth] = optPos;
    extendable[optPos] = extendable[posExtendable - 1];
    posExtendable--;

#ifdef MAPPING_FUNCTION_LOG
    printf("\t(P) weight(u%d): %ld in pos %d in Q (|Q|: %d)\n", current,
           optWeight, optPos, posExtendable);
#endif
  }

  inline void removeFromQueue(int depth) {
    posExtendable -= nInserted[depth];
#ifdef MAPPING_FUNCTION_LOG
    printf("\t(R) |Q|: %d\n", posExtendable);
#endif
  }

  inline void clearQueue() {
    memset(nInserted, 0, sizeof(int) * nQueryVertex);
    posExtendable = 0;
#ifdef MAPPING_FUNCTION_LOG
    printf("\t(R) |Q|: %d\n", posExtendable);
#endif
  }

  void clear_nInserted(int depth) { this->nInserted[depth] = 0; }

  void add_nInserted(int depth) { this->nInserted[depth]++; }

  int set_optWeight(int optWeight) { this->optWeight = optWeight; }

  int set_minPosition(int minPosition) { this->minPosition = minPosition; }
};
