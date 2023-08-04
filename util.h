// Prints the DAG representation of a query graph.
//
// This function takes a query graph and prints
// its DAG representation to the standard error stream.
// It iterates through the vertices of the query graph,
// displaying each vertex's label and its children in the DAG, if any.
void printDAG(Graph& query) {
  cerr << "Query DAG" << endl;
  for (int i = 0; i < nQueryVertex; i++) {
    int u = uSequence[i];
    cerr << "[" << i << "] l(u" << u << "): " << query.label[u] << ". Child(u"
         << u << "):";
    for (int j = 0; j < DAG_child_query_size[u]; j++) {
      cerr << " u" << DAG_child_query[u][j];
    }
    cerr << endl;
  }
}

// Prints information about the CS(Candidate Set)
// and adjacency lists for nodes in a query graph.
//
// This function takes a query graph and prints information
// related to the CS and adjacency lists of the graph's nodes.
// It displays the core nodes of the graph, their corresponding candidate sets,
// and the adjacency lists of the nodes.
void printCS(Graph& query) {
  cerr << "Core:";
  for (int i = 0; i < nQueryVertex; i++) {
    if (query.core[i] > 1) cerr << " u" << i;
  }
  cerr << endl;
  cerr << "CS" << endl;
  for (int i = 0; i < nQueryVertex; i++) {
    int currNode = i;
    // NEC boost
    if (query.NECMap[currNode] != -1 && query.NECMap[currNode] != currNode) {
      cerr << "C(u" << currNode << "): same as C(u" << query.NECMap[currNode]
           << "(NEC)" << endl;
      continue;
    }
    CandidateSpace& currSet = candSpace[currNode];
    cerr << "C(u" << currNode << ")(|C(u)|: ";
    int size = 0;
    for (int j = 0; j < currSet.size; j++) {
      if (currSet.candidates[j] != -1) size++;
    }
    cerr << size << "):";
    for (int j = 0; j < currSet.size; j++) {
      if (currSet.candidates[j] == -1) continue;
#ifdef PRUNING_BY_EQUIVALENCE_SETS
      cerr << " v" << currSet.candidates[j] << "(" << j << ", "
           << currSet.cell[j] << ")";
#else
      cerr << " v" << currSet.candidates[j] << "(" << j << ")";
#endif
    }
    cerr << endl;
  }
  cerr << "Adjacency list" << endl;
  for (int i = 0; i < nQueryVertex; i++) {
    int currNode = i;
    // NEC boost
    if (query.NECMap[currNode] != -1 && query.NECMap[currNode] != currNode) {
      int parentNode = DAG_parent_query[currNode][0];
      cerr << "N^{u" << parentNode << "}_{u" << currNode << "]: same as N^{u"
           << parentNode << "}_{u" << query.NECMap[currNode] << "} (NEC)"
           << endl;
      continue;
    }
    for (int parentIndex = 0; parentIndex < DAG_parent_query_size[currNode];
         parentIndex++) {
      int parentNode = DAG_parent_query[currNode][parentIndex];
      CandidateSpace& currSet = candSpace[currNode];
      CandidateSpace& parentSet = candSpace[parentNode];
      for (int j = 0; j < parentSet.size; j++) {
        if (parentSet.candidates[j] == -1) continue;
        cerr << "N^{u" << parentNode << "}_{u" << currNode << "} (v"
             << parentSet.candidates[j] << "): ";
        for (int k = 0; k < currSet.nAdjacent[parentIndex][j]; k++) {
          cerr << "v"
               << currSet.candidates[currSet.adjacent[parentIndex][j][k]];
          if (k != currSet.nAdjacent[parentIndex][j] - 1) cerr << ", ";
        }
        cerr << endl;
      }
    }
  }
}
