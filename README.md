# VEQ
## 1. Background and Objectives
Our research team has developed an algorithm that is over 400 times faster than the previous state-of-the-art algorithms for the subgraph query processing problem. We presented this algorithm at SIGMOD, which is one of the top computer science conferences. As a result of this research, we have received inquiries and requests for code sharing from researchers at prestigious universities overseas. Additionally, a startup in the field of drug development, Inc. AIgenDrug, has also requested code sharing.  
  
Therefore, the goal of this project is to refine the source code of VEQ, which is the world's top-performing algorithm for the subgraph query processing and subgraph matching, and make it publicly available as a library for more convenient use.

## 2. Development Environment and Programming Language
Programming Language:
- C++ 14

Development Environment:
- Operating System: macOS Ventura
- Integrated Development Environment (IDE): Visual Studio Code

## 3. System Configuration and Architecture
This project includes headers inside the `include` directory, `main.cpp` inside the `src` directory, and other graph files and test data. `VEQ.h` provides the core algorithm of the project, `config.h` contains necessary settings for program execution, `memory.h`, `structure.h`, and `util.h` are responsible for memory allocation, data structure declarations, and declaration of other inline functions in preprocessing stages. `graph.h` contains the data structure forms for graphs and candidate spaces.

During program execution and testing, necessary functionalities are provided by `run.h` and `compare.h.` The figure below illustrates the essential steps of the algorithm.

![VEQS_architecture](https://github.com/SNUCSE-CTA/VEQ_S/assets/83649602/76422cab-ce5e-42b5-b02f-27186859b445)  
*Figure 1. ProcessQuery execution process in main and associated header file structure*

## 4. Key Features of the Project
The VEQ algorithm utilizes static equivalence, dynamic equivalence, and neighbor-safety to significantly enhance query processing performance and enable efficient search on large-scale data graphs.

### 4.1. Subgraph Query Processing Problem
The subgraph query processing problem involves finding all data graphs that contain subgraphs isomorphic to a given query graph among multiple data graphs and one query graph. For example, in Figure 2, given a set of data graphs $D=\{G_1, G_2\}$ and a query graph $q$, $G_1$ contains a subgraph isomorphic to $q$, but $G_2$ does not. Therefore, the answer is $\{G_1\}$.

![Example of Subgraph Query Processing Problem](https://github.com/SNUCSE-CTA/VEQ_S/assets/83649602/fb8b4d4d-bc01-44d5-9204-b630bf2efd87)  
*Figure 2. Example of the Subgraph Query Processing Problem*

### 4.2. Overview of the VEQ Algorithm
VEQ [1, 2] follows the filtering-verification approach and consists of three stages: (1) Building a query DAG, (2) Building Candidate Space (CS), and (3) Matching.

#### 4.2.1. Building a query DAG
In this stage, a query DAG $q_d$ is created by assigning directions to edges of the query graph $q$. The vertex in $q$ with an infrequent label and a large degree is selected as the root, and a breadth-first search is performed starting from the root. During this process, neighbor equivalence classes are found for every vertex with degree 1, and vertices belonging to the same neighbor equivalence class are merged into one vertex.

#### 4.2.2. CS Construction Stage
For each data graph and query graph, an auxiliary data structure called CS is constructed. CS consists of candidate sets $C(u)$ for each vertex $u$ in the query graph and edges between vertices in the candidate vertex set. While DAF, a previous subgraph isomorphism algorithm, proposed CS [3], VEQ reduces the candidate set size using neighbor safety (Figure 3). If the candidate set for a vertex in the query graph is empty, there is no isomorphic subgraph to the query graph. Therefore, the algorithm proceeds to the next data graph. If not, the algorithm moves on to the exploration stage.

![Reducing Candidate Sets](https://github.com/SNUCSE-CTA/VEQ_S/assets/83649602/27b65097-93cc-44e9-a33a-07d1ca6fddc8)  
*Figure 3. Reducing Candidate Sets*

#### 4.2.3. Exploration Stage
For each data graph's CS, the algorithm checks if an isomorphic subgraph to the query graph exists. Each vertex $u$ in the query graph is mapped to vertices in $C(u)$, and during this process, techniques such as failing sets and dynamic equivalence are used to eliminate duplicate subgraphs. If at least one isomorphic subgraph to the query graph is found for a data graph, it is added to the answer set, and the algorithm proceeds to the next data graph.

The project repository supports extensive experiments on well-known real datasets using various tests and test data. It confirms that the algorithm outperforms existing algorithms in multiple aspects.

Additionally, continuous integration (CI) is used to execute scripts that check if the program correctly reads graphs, providing automated testing for this script, enhancing stability and reliability.

## 5. How to run
Compile and make an executable binary
```bash
make
```
Run VEQs and VEQm
```bash
./VEQ_S -dg [data graph file] -qg [query graph file] -o [output file]
./VEQ_M -dg [data graph file] -qg [query graph file]
```
(Example)
```bash
./VEQ_S -dg graph/search/data/COLLAB.gfu -qg graph/search/query/COLLAB/bfs/8/q0.gfu
./VEQ_M -dg graph/matching/data/yeast.gfu -qg graph/matching/query/yeast/sparse/50/q0.gfu
```

Run VEQs and VEQm using example data and query
```bash
make run1
make run2
```
Run tests
```bash
make test
```
Remove all binaries and object files
```bash
make clean
```
## 6. Library Build and Python Binding

## 6.1. Requirements
- pybind11 is required for Python binding. pybind11 installation guide is available on the [official document](https://pybind11.readthedocs.io/en/stable/installing.html).
- Adding the path for pybind11 to the include path may be required. This can be done by executing the following command. `export CPLUS_INCLUDE_PATH=<path-to-pybind11>:$CPLUS_INCLUDE_PATH`. Alternatively, the following argument can be added to the compile command in `build_library.sh`. `-I"<path-to-pybind11>"`. The argument `<path-to-pybind11>` depends on the system. In the case of the development environment, it was `/opt/homebrew/Cellar/pybind11/2.11.1/include`.
- Adding python related headers to the include path may also be required. This can be done by executing the following command. `export CPLUS_INCLUDE_PATH=<path-to-required-headers>:$CPLUS_INCLUDE_PATH`. Alternatively, the following argument can be added to the compile command in `build_library.sh`. `-I"<path-to-required-headers>"`. The argument `<path-to-required-headers>` depends on the system. In the case of the development environment, it was `/Library/Developer/CommandLineTools/Library/Frameworks/Python3.framework/Versions/3.8/Headers`.

## 6.2. Building the library
Build the library
```bash
./build_library.sh
```
Alternatively, the library can also be built using the following method.
```bash
cd libVEQ
make
cd ..
```
Through the provided `run` function within the `libVEQ.so` library, the same functionality as the `VEQ_S` binary can be tested. This feature is available in `libVEQ/binding.py`, and its usage is identical to `VEQ_S`. Below is an example execution command.
```bash
python3 libVEQ/binding.py -dg graph/data/COLLAB.gfu -qg graph/query/COLLAB/randomwalk/8/q30.gfu
```
## 7. Expected Impact and Application Areas
Currently, big data is being generated at a rapid pace, and services utilizing this data are advancing significantly. The generation and analysis of graph big data are also becoming active, leading to an increasing demand for efficient algorithms for large-scale graph problems. The technology developed in this project can be used for tasks such as searching for specific patterns in big data graphs and aiding in drug development. This graph analysis technology is expected to facilitate toxicological analysis, a crucial step in drug development. In fact, the startup Inc. AIgenDrug has successfully utilized the technology from this project to quickly address toxicity analysis of compound data.

## 8. References
[1] H. Kim, Y. Choi, K. Park, X. Lin, S.-H. Hong, and W.-S. Han. 2021. Versatile equivalences: Speeding up subgraph query processing and subgraph matching. In Proceedings of ACM SIGMOD.

[2] H. Kim, Y. Choi, K. Park, X. Lin, S.-H. Hong, and W.-S. Han. 2023. Fast subgraph query processing and subgraph matching via static and dynamic equivalences. The VLDB Journal, 32.

[3] M. Han, H. Kim, G. Gu, K. Park, and W.-S. Han. 2019. Efficient subgraph matching: Harmonizing dynamic programming, adaptive matching order, and failing set together. In Proceedings of ACM SIGMOD.
