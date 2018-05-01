
#include <string>
#include <cstring>
#include <fstream>

#include "geometry.hpp"
#include "graph.hpp"

#ifndef _PBBS_SPTL_READINPUTBINARY
#define _PBBS_SPTL_READINPUTBINARY

namespace sptl {

template <class Item>
Item read_from_file(std::ifstream& in);

template <class Item>
struct read_from_file_struct {
  Item operator()(std::ifstream& in) const {
    Item memory;
    in.read(reinterpret_cast<char*>(&memory), sizeof(Item));
    return memory;
  }
};

template <>
struct read_from_file_struct<std::string> {
  std::string operator()(std::ifstream& in) const {
    int size = 0;
    in.read(reinterpret_cast<char*>(&size), sizeof(int));
    std::string answer;
    answer.resize(size);
    in.read(&answer[0], size);
    return answer;
  }
};

template <class Item>
struct read_from_file_struct<Item*> {
  Item* operator()(std::ifstream& in, long size) const {
    Item* result = (Item*)malloc(sizeof(Item) * size);
    in.read(reinterpret_cast<char*>(result), sizeof(Item) * size);
    return result;
  }
};

template <class Item>
struct read_from_file_struct<parray<Item>> {
  parray<Item> operator()(std::ifstream& in) const {
    long size = 0;
    in.read(reinterpret_cast<char*>(&size), sizeof(long));
    parray<Item> result(size);
    in.read(reinterpret_cast<char*>(result.begin()), sizeof(Item) * size);
    return result;
  }
};

template <>
struct read_from_file_struct<parray<char*>> {
  parray<char*> operator()(std::ifstream& in) const {
    long size = 0;
    in.read(reinterpret_cast<char*>(&size), sizeof(long));
    parray<char*> result(size);
    int* len = new int[size];
    in.read(reinterpret_cast<char*>(len), sizeof(int) * size);
    for (int i = 0; i < size; i++) {
      result[i] = new char[len[i] + 1];
      result[i][len[i]] = 0;
      in.read(&result[i][0], sizeof(char) * len[i]);
    }
    delete [] len;
    return result;
  }
};

template <>
struct read_from_file_struct<parray<std::pair<char*, int>*>> {
  parray<std::pair<char*, int>*> operator()(std::ifstream& in) const {
    long size = 0;
    in.read(reinterpret_cast<char*>(&size), sizeof(long));
    int* len = new int[size];
    in.read(reinterpret_cast<char*>(len), sizeof(int) * size);
    parray<std::pair<char*, int>*> result(size);
    for (int i = 0; i < size; i++) {
      char* f = new char[len[i] + 1];
      f[len[i]] = 0;
      in.read(f, sizeof(char) * len[i]);
      int s = 0;
      in.read(reinterpret_cast<char*>(&s), sizeof(int));
      result[i] = new std::pair<char*, int>(f, s);
    }
    delete [] len;
    return result;
  }
};

template <class Point>
struct read_from_file_struct<triangles<Point>> {
  triangles<Point> operator()(std::ifstream& in) const {
    triangles<Point> t;
    in.read(reinterpret_cast<char*>(&t.num_points), sizeof(long));
    t.p = read_from_file_struct<Point*>()(in, t.num_points);
    in.read(reinterpret_cast<char*>(&t.num_triangles), sizeof(long));
    t.t = read_from_file_struct<triangle*>()(in, t.num_triangles); return t;
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

class ray_cast_test {
public:
  static
  ray_cast_test foo;
  sptl::parray<point3d> points;
  sptl::parray<triangle> triangles;
  sptl::parray<ray<point3d>> rays;

  ray_cast_test(): points(0), triangles(0), rays(0) { }

  ray_cast_test(const ray_cast_test& other) {
    points = other.points;
    triangles = other.triangles;
    rays = other.rays;
  }

  ray_cast_test& operator=(const ray_cast_test& other) {
    if (&other == this) {
      return *this;
    }
    this->points = other.points;
    this->triangles = other.triangles;
    this->rays = other.rays;
    return foo;
  }
};

ray_cast_test ray_cast_test::foo;

template <>
struct read_from_file_struct<ray_cast_test> {
  ray_cast_test operator()(std::ifstream& in) {
    ray_cast_test test;
  
    test.points = read_from_file<parray<point3d>>(in);
    test.triangles = read_from_file<parray<triangle>>(in);
    test.rays = read_from_file<parray<ray<point3d>>>(in);
  
    return test;
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
  
} // end namespace

#endif
