#include "libVEQ.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

void run(int argc, char *argv[]) {
  setvbuf(stdout, NULL, _IOLBF, 0);
  int queryFormat = GFU_FORMAT;
  int dataFormat = GFU_FORMAT;
  for (int i = 1; i < argc; ++i) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
        case 'd':
          switch (argv[i][2]) {
            case 'i':
              dataFormat = IGRAPH_FORMAT;
              break;
            case 'g':
              dataFormat = GFU_FORMAT;
              break;
          }
          dataGraphFile = argv[i + 1];
          break;
        case 'q':
          switch (argv[i][2]) {
            case 'c':
              queryFormat = CFL_FORMAT;
              break;
            case 'i':
              queryFormat = IGRAPH_FORMAT;
              break;
            case 'g':
              queryFormat = GFU_FORMAT;
              break;
          }
          queryGraphFile = argv[i + 1];
          break;
        case 'o':
          answerFile = argv[i + 1];
          break;
      }
    }
  }
  cout << "Data file: " << dataGraphFile << endl;
  cout << "Query file: " << queryGraphFile << endl;
  cout << "Output file: " << answerFile << endl;
  if (dataFormat == IGRAPH_FORMAT)
    ReadIgraphFormat(dataGraphFile, dataGraph);
  else
    ReadGFUFormat(dataGraphFile, dataGraph);

  ProcessDataGraphs();
  AllocateForDataGraph();

  if (queryFormat == CFL_FORMAT)
    ReadCFLFormat(queryGraphFile, queryGraph);
  else if (queryFormat == IGRAPH_FORMAT)
    ReadIgraphFormat(queryGraphFile, queryGraph);
  else
    ReadGFUFormat(queryGraphFile, queryGraph);

  ProcessQuery();
}

PYBIND11_MODULE(libVEQ, m) {
  m.doc() = "Test";
  m.def("run", [](vector<string> args) {
    vector<char *> cstrs;
    cstrs.reserve(args.size());
    for (auto &s : args) cstrs.push_back(const_cast<char *>(s.c_str()));
    return run(cstrs.size(), cstrs.data());
  });
}
