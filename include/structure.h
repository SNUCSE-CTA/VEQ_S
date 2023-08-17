#include "util.h"

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

struct Stack {
  int* address = NULL;
  int addressSize;
  int addressPos;
  int vertex;
  // for failing set
  uint64_t* failingSet;
  uint64_t* buffer;
  uint64_t* conflictCell;
  bool* isPruned;
  int nPruned;
  int* pruned;
  int* problemChild;
  int nProblemChild;
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
  inline void reinsertToQueue(int u, int depth) {
    int position = positions[depth];
    extendable[posExtendable] = extendable[position];
    posExtendable++;
    extendable[position] = u;
  }

  inline void insertToQueue(int u) {
    extendable[posExtendable] = u;
    posExtendable++;
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
  }

  inline void removeFromQueue(int depth) { posExtendable -= nInserted[depth]; }

  inline void clearQueue() {
    memset(nInserted, 0, sizeof(int) * nQueryVertex);
    posExtendable = 0;
  }

  void clear_nInserted(int depth) { this->nInserted[depth] = 0; }

  void add_nInserted(int depth) { this->nInserted[depth]++; }

  int set_optWeight(int optWeight) { this->optWeight = optWeight; }

  int set_minPosition(int minPosition) { this->minPosition = minPosition; }
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
