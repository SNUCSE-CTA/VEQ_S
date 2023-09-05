# VEQs
## 1. Background and Objectives
With the global attention on big data, related research and technologies are rapidly growing. Among the methods for storing big data, graphs are one of the prominent ones. In particular, the subgraph query processing problem involves finding all data graphs that contain subgraphs isomorphic to a given query graph among multiple data graphs and one query graph. This problem has practical applications in various fields such as social network analysis and chemical compound search.

The aim of this project is to refine the source code of the latest algorithm for the subgraph query processing problem, VEQs, and make it available as a library for convenient use.

## 2. Development Environment and Programming Language
Programming Language:
- C++ 14

Development Environment:
- Operating System: macOS Ventura
- Integrated Development Environment (IDE): Visual Studio Code

## 3. System Configuration and Architecture
This project includes headers inside the `include` directory, `main.cpp` inside the `src` directory, and other graph files and test data. `VEQ.h` provides the core algorithm of the project, `config.h` contains necessary settings for program execution, `memory.h`, `structure.h`, and `util.h` are responsible for memory allocation, data structure declarations, and declaration of other inline functions in preprocessing stages. `graph.h` contains the data structure forms for graphs and candidate spaces.

During program execution and testing, necessary functionalities are provided by `run.h` and `compare.h.` The figure below illustrates the essential steps of the algorithm.

![VEQS_architecture](https://github.com/SNUCSE-CTA/VEQ_S/assets/69783927/88921009-f343-40d2-a0d2-e9d64ca05e5c)
*Figure 1. ProcessQuery execution process in main and associated header file structure*

## 4. Key Features of the Project
The VEQs algorithm utilizes static equivalence, dynamic equivalence, and neighbor-safety to significantly enhance query processing performance and enable efficient search on large-scale data graphs.

### 4.1. Subgraph Query Processing Problem
The subgraph query processing problem involves finding all data graphs that contain subgraphs isomorphic to a given query graph among multiple data graphs and one query graph. For example, in Figure 2, given a set of data graphs $D=\{G_1, G_2\}$ and a query graph $q$, $G_1$ contains a subgraph isomorphic to $q$, but $G_2$ does not. Therefore, the answer is $\{G_1\}$.

![Example of Subgraph Query Processing Problem](https://github.com/SNUCSE-CTA/VEQ_S/assets/69783927/d3de2e66-08d2-4241-a2a0-712317177b2d)
*Figure 2. Example of the Subgraph Query Processing Problem*

### 4.2. Overview of the VEQs Algorithm
VEQs [1, 2] follows the filtering-verification approach and consists of three stages: (1) Building a query DAG, (2) Building Candidate Space (CS), and (3) Matching.

#### 4.2.1. Building a query DAG
In this stage, a query DAG $q_d$ is created by assigning directions to edges of the query graph $q$. The vertex in $q$ with an infrequent label and a large degree is selected as the root, and a breadth-first search is performed starting from the root. During this process, neighbor equivalence classes are found for every vertex with degree 1, and vertices belonging to the same neighbor equivalence class are merged into one vertex.

#### 4.2.2. CS Construction Stage
For each data graph and query graph, an auxiliary data structure called CS is constructed. CS consists of candidate sets $C(u)$ for each vertex $u$ in the query graph and edges between vertices in the candidate vertex set. While DAF, a previous subgraph isomorphism algorithm, proposed CS [3], VEQ reduces the candidate set size using neighbor safety (Figure 3). If the candidate set for a vertex in the query graph is empty, there is no isomorphic subgraph to the query graph. Therefore, the algorithm proceeds to the next data graph. If not, the algorithm moves on to the exploration stage.

![Reducing Candidate Sets](https://github.com/SNUCSE-CTA/VEQ_S/assets/69783927/d8e8d1fa-5633-41fa-abca-c4d78fd4c974)
*Figure 3. Reducing Candidate Sets*

#### 4.2.3. Exploration Stage
For each data graph's CS, the algorithm checks if an isomorphic subgraph to the query graph exists. Each vertex $u$ in the query graph is mapped to vertices in $C(u)$, and during this process, techniques such as failing sets and dynamic equivalence are used to eliminate duplicate subgraphs. If at least one isomorphic subgraph to the query graph is found for a data graph, it is added to the answer set, and the algorithm proceeds to the next data graph.

The project repository supports extensive experiments on well-known real datasets using various tests and test data. It confirms that the algorithm outperforms existing algorithms in multiple aspects.

Additionally, continuous integration (CI) is used to execute scripts that check if the program correctly reads graphs, providing automated testing for this script, enhancing stability and reliability.

## 5. Expected Impact and Application Areas
In this project, the query execution time improvement for the subgraph query processing problem is an astounding 41741%. These research results have garnered inquiries and requests for code sharing from researchers at over 30 prestigious overseas universities, including the University of Edinburgh, New York University, Arizona State University, Hong Kong University of Science and Technology, Chinese University of Hong Kong, Peking University, Fudan University, Osaka University, University of Verona, University of Salerno, Eindhoven University of Technology, and University of Sydney. Furthermore, requests for code sharing have also been received from NTT, a telecommunications company in Japan, and Inc. AIgenDrug, a startup in the field of drug development.

Currently, big data is being generated rapidly, and services utilizing it are advancing significantly. The generation and analysis of graph big data are also becoming active, increasing the demand for efficient algorithms for large-scale graph problems. The technology from this project can be used for searching specific patterns in big data graphs, detecting cyberattacks, and aiding in drug development. The use of this graph analysis technology is expected to facilitate toxicological analysis, a critical step in drug development, potentially benefiting the field of drug development.

## 6. References
[1] H. Kim, Y. Choi, K. Park, X. Lin, S.-H. Hong, and W.-S. Han. 2021. Versatile equivalences: Speeding up subgraph query processing and subgraph matching. In Proceedings of ACM SIGMOD.

[2] H. Kim, Y. Choi, K. Park, X. Lin, S.-H. Hong, and W.-S. Han. 2023. Fast subgraph query processing and subgraph matching via static and dynamic equivalences. The VLDB Journal, 32.

[3] M. Han, H. Kim, G. Gu, K. Park, and W.-S. Han. 2019. Efficient subgraph matching: Harmonizing dynamic programming, adaptive matching order, and failing set together. In Proceedings of ACM SIGMOD.
