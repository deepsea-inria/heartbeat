#include <iostream>
#include <chrono>
#include <array>
#include <fstream>
#include <iostream>
#include <istream>
#include <ostream>
#include <sstream>
#include <algorithm>

#include "graph.h"

#ifndef _HEARTBEAT_GRAPHIO_H_
#define _HEARTBEAT_GRAPHIO_H_

namespace pbbs {

template <class Item>
struct read_from_file_struct {
  Item operator()(std::ifstream& in) const {
    Item memory;
    in.read(reinterpret_cast<char*>(&memory), sizeof(Item));
    return memory;
  }
};
  
template <class intT>
struct read_from_file_struct<graph::graph<intT>> {
  graph::graph<intT> operator()(std::ifstream& in) {
    intT n, m;
    in.read(reinterpret_cast<char*>(&n), sizeof(intT));
    in.read(reinterpret_cast<char*>(&m), sizeof(intT));
    intT* degree = new intT[n];
    in.read(reinterpret_cast<char*>(degree), sizeof(intT) * n);
    intT* e = (intT*)malloc(sizeof(intT) * m);
    in.read(reinterpret_cast<char*>(e), sizeof(intT) * m);
    graph::vertex<intT>* v = (graph::vertex<intT>*)malloc(sizeof(graph::vertex<intT>) * n);
    int offset = 0;
    for (int i = 0; i < n; i++) {
      v[i] = graph::vertex<intT>(e + offset, degree[i]);
      offset += degree[i];
    }
    delete [] degree;
    return graph::graph<intT>(v, n, m, e);
  }
};

template <class Item>
Item read_from_file(std::ifstream& in) {
  return read_from_file_struct<Item>()(in);
}
  
template <class Item>
Item read_from_file(std::string file) {
  std::ifstream in(file, std::ifstream::binary);
  return read_from_file<Item>(in);
}

} // namespace

#endif
