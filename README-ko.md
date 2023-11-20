# VEQs
## 1. 개발배경 및 목적
본 연구팀은 부분그래프 질의 처리 문제에 대해서 이전 최고 성능 알고리즘에 비해 400배 이상 빠 른 알고리즘을 개발하여, 컴퓨터 분야 최우수학술대회인 SIGMOD에서 발표하였다. 이러한 연구 결과에 대해 해외 유수 대학의 연구진들로부터 알고리즘에 대한 문의와 코드 공유 요청을 받았고, 또한 신약 개발 관련 스타트업 (주)AIgenDrug으로부터도 코드 공유 요청을 받았다.  
  
이에 본 프로젝트는 부분 그래프 질의 문제에 대한 세계 최고 성능의 알고리즘인 VEQs의 소스코드 를 다듬어 공개하고, 이를 라이브러리화하여 보다 편리하게 사용할 수 있도록 만드는 것을 목표로 한 다.

## 2. 개발환경 및 개발언어
개발 언어
- C++ 14

개발 환경
- 운영체제: macOS Ventura
- 통합 개발 환경(IDE): Visual Studio Code

## 3. 시스템 구성 및 아키텍처
본 프로젝트는 include 디렉토리 내부의 헤더와 src 디렉토리 내부의 main.cpp, 그 외 그래프 파일과 테스트 데이터를 포함한다. VEQ.h는 프로젝트의 핵심 알고리즘을 제공하며, config.h는 프로그램 수행에 필요한 설정, memory.h, structure.h 및 util.h는 각각 전처리 과정에 필요한 메모리 할당, 자료구조의 선언, 기타 inline function의 선언을 담당한다. 또한 graph.h는 그래프 및 candidate space의 자료구조 형태를 담고 있다.

프로그램의 실행과 테스트 수행 시에는 각각 run.h과 compare.h에서 필요한 기능을 제공한다. 아래 그림은 핵심적인 알고리즘의 수행 단계를 보여준다.

![VEQS_architecture](https://github.com/SNUCSE-CTA/VEQ_S/assets/83649602/76422cab-ce5e-42b5-b02f-27186859b445)  
*그림 1. main에서 수행되는 ProcessQuery의 과정 및 연관된 헤더 파일 구조*

## 4. 프로젝트 주요기능
VEQs 알고리즘은 static equivalence와 dynamic equivalence, 그리고 이웃 안정성 (neighbor-safety)을 활용해 query 처리 성능을 혁신적으로 향상시키고, 대규모 data graph에 대한 효율적인 검색을 가능하게 한다.

### 4.1. 부분 그래프 질의 문제
부분 그래프 질의 문제는 다수의 데이터 그래프와 하나의 쿼리 그래프가 주어졌을 때 쿼리 그래프와 동형(isomorphic)인 부분 그래프를 포함하는 모든 데이터 그래프를 찾아내는 문제이다. 예를 들어, 그림 2에서 데이터 그래프 집합 $D=\{G_1, G_2\}$와 쿼리 그래프 $q$가 주어졌을 때, $G_1$은 $q$와 동형인
부분 그래프를 가지지만 $G_2$는 그렇지 않으므로, 정답은 $\{G_1\}$이 된다.

![부분 그래프 질의 문제의 예시](https://github.com/SNUCSE-CTA/VEQ_S/assets/83649602/fb8b4d4d-bc01-44d5-9204-b630bf2efd87)  
*그림 2. 부분 그래프 질의 문제의 예시*

### 4.2. VEQs 알고리즘 개요
VEQs [1, 2]는 filtering-verification 방법을 따르는 알고리즘으로, (1) 쿼리 DAG 구축, (2) CS(Candidate Space) 구축, (3) 탐색 단계로 구성되어 있다.

#### 4.2.1. 쿼리 DAG 구축 단계
이 단계에서는 쿼리 그래프 $q$의 간선들에 방향을 부여하여 쿼리 DAG $q_d$를 만든다. 이때 $q$에서 가장 희소한 레이블(label)을 가지면서 가장 큰 차수(degree)를 가진 정점을 루트(root)로 하고, 그로부터 너비 우선 탐색(breadth first search)을 수행한다. 이 과정에서 차수가 1인 모든 정점에 대해 이웃 동치류(neighbor equivalence class)를 찾고 동일한 이웃 동치류에 속하는 정점들은 하나의 정점으로 묶는다.

#### 4.2.2. CS 구축 단계
각 데이터 그래프와 쿼리 그래프에 대해 보조 자료구조인 CS를 만든다. CS는 쿼리 그래프의 각 정점 $u$에 대한 후보 집합(candidate set) $C(u)$와 후보 정점 집합의 정점들 사이의 간선으로 구성된다. CS는 이전에 부분 그래프 동형 알고리즘인 DAF가 제안하였으나 [3], VEQ는 이웃 안전성을 사용해 DAF의 것보다 작은 후보 집합을 얻는다 (그림 3). 쿼리 그래프의 어느 한 정점에 대해서라 도 그 후보 집합이 비어 있다면 동형인 쿼리 그래프와 동형인 부분 그래프가 없으므로 다음 데이터 그래프로 진행하고, 그렇지 않다면 탐색하기 단계로 넘어간다.

![후보 집합 줄이기](https://github.com/SNUCSE-CTA/VEQ_S/assets/83649602/27b65097-93cc-44e9-a33a-07d1ca6fddc8)  
*그림 3. 후보 집합 줄이기*

#### 4.2.3. 탐색 단계
각 데이터 그래프의 CS에 대해 쿼리 그래프와 동형인 부분 그래프가 존재하는지 확인한다. 쿼리 그래프의 각 정점 $u$를 $C(u)$에 속한 정점에 대응시키며, 그 과정에서 실패 집합(failing sets) 기법과 동적 동치(dynamic equivalence)를 통해 중복된 부분 트리를 제거하는 기법을 사용한다. 쿼리 그래프와 동형인 부분 그래프가 하나라도 발견된다면 해당 데이터 그래프를 정답 집합에 추가하고, 다음 데이터 그래프로 진행한다.

프로젝트 레포지토리에서는 여러 test 및 testdata를 통해 잘 알려진 real dataset에 대해 광범위한 실험을 수행하도록 지원하며, 기존의 알고리즘을 여러 방면에서 능가함을 확인할 수 있다.

또 CI(continuous integration)를 사용하여 로컬 테스트와는 별개로 그래프를 잘 읽어오는지 확인하는 스크립트를 실행시킨다. 이는 해당 스크립트에 대한 자동화된 테스트를 제공하며, 안정성 및 신뢰성을 향상시키는 데 도움을 준다.

## 5. 기대효과 및 활용분야
현재 빠른 속도로 빅 데이터가 생성되고 있고 이를 활용한 서비스가 크게 발전하고 있다. 그래프 빅데이터의 생성과 분석도 활성화되어 대규모 그래프 문제에 대한 효율적인 알고리즘에 대한 수요가 증가하고 있다. 본 프로젝트의 기술은 빅 데이터 그래프에서 특정한 패턴 검색, 신약 개발 등에 사용될 수 있다. 이 그래프 분석 기술을 사용하면 신약개발의 중요 과정인 독성 분석이 용이해지고 이로 인해 신약개발에도 도움이 될 것으로 기대된다. 실제로, 신약 개발 관련 스타트업 (주)AIgenDrug이 본 프로젝트의 기술을 사용하여 화합물 데이터의 독성 분석을 빠른 시간에 해결할 수 있었다.

## 6. 참고문헌
[1] H. Kim, Y. Choi, K. Park, X. Lin, S.-H. Hong, and W.-S. Han. 2021. Versatile equivalences: Speeding up subgraph query processing and subgraph matching. In Proceedings of ACM SIGMOD.

[2] H. Kim, Y. Choi, K. Park, X. Lin, S.-H. Hong, and W.-S. Han. 2023. Fast subgraph query processing and subgraph matching via static and dynamic equivalences. The VLDB Journal, 32.

[3] M. Han, H. Kim, G. Gu, K. Park, and W.-S. Han. 2019. Efficient subgraph matching: Harmonizing dynamic programming, adaptive matching order, and failing set together. In Proceedings of ACM SIGMOD.

