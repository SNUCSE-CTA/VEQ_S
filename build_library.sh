cd libVEQ
g++ -std=c++14 -march=native -flto -O3 -w -DNDEBUG -shared -fPIC -Wl,-undefined,dynamic_lookup -o libVEQ.so libVEQ.cpp
