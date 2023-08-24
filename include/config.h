typedef long long weight_type;
typedef unsigned short degree_type;
#ifdef HUGE_GRAPH
typedef long long number_type;
#else
typedef int number_type;
#endif
// Constants
#ifdef HUGE_GRAPH
const int BITS_PER_LABEL = 1;
#else
const int BITS_PER_LABEL =
    4;  // the maximum number of labels per label that we count as NLF
#endif
const int Start = 0;
const int LoadNumVertex = 1;
const int LoadVertex = 2;
const int LoadNumEdge = 3;
const int LoadEdge = 4;
const int GFU_FORMAT = 5;
const int IGRAPH_FORMAT = 6;
const int CFL_FORMAT = 7;
const int MAX_NUM_VERTEX = 257;
const int MAX_QUERY_DEGREE = 100;
const int MAX_NUM_LABEL = 48000;
const int UINT64_SIZE = sizeof(uint64_t) * 8;
// Execution option flags
const bool usePathSizeOrder = false;
const bool useFailingSet = true;
bool isTree;

int NLFSize = -1;
int nUniqueLabel = 0;
