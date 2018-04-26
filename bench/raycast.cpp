#include "heartbeatbench.hpp"
#include "readinputbinary.hpp" 
#include "kdtree.hpp"
#include "kdTree.h"

namespace heartbeatbench {

template <class Item>
using parray = sptl::parray<Item>;

pbbs::_point3d<double> to_pbbs(sptl::_point3d<double> point) {
  return pbbs::_point3d<double>(point.x, point.y, point.z);
}

pbbs::_vect3d<double> to_pbbs(sptl::_vect3d<double> vect) {
  return pbbs::_vect3d<double>(vect.x, vect.y, vect.z);
}

pbbs::triangle to_pbbs(sptl::triangle t) {
  return pbbs::triangle(t.vertices[0], t.vertices[1], t.vertices[2]);
}

pbbs::ray<pbbs::_point3d<double>> to_pbbs(sptl::ray<sptl::_point3d<double>> ray) {
  return pbbs::ray<pbbs::_point3d<double>>(to_pbbs(ray.o), to_pbbs(ray.d));
}

template <class Item1, class Item2>
parray<Item1> to_pbbs(parray<Item2>& a) {
  parray<Item1> result(a.size());
  for (int i = 0; i < a.size(); i++) {
    result[i] = to_pbbs(a[i]);
  }
  return result;
}

void benchmark(sptl::ray_cast_test& x) {
  deepsea::cmdline::dispatcher d;
  parray<pbbs::_point3d<double>> points = to_pbbs<pbbs::_point3d<double>>(x.points);
  parray<pbbs::triangle> triangles = to_pbbs<pbbs::triangle>(x.triangles);
  parray<pbbs::ray<pbbs::_point3d<double>>> rays = to_pbbs<pbbs::ray<pbbs::_point3d<double>>>(x.rays);
  pbbs::triangles<pbbs::_point3d<double>> tri(points.size(), triangles.size(), points.begin(), triangles.begin());
  d.add("heartbeat", [&] {
    intT* dest;
    heartbeat::launch_interpreter<heartbeatbench::rayCast>(tri, rays.begin(), rays.size(), &dest);
  });
  d.add("pbbs", [&] {
    heartbeatbench::run_and_report_elapsed_time([&] {
      pbbs::rayCast(tri, rays.begin(), rays.size());
    });
  });
  d.dispatch("algorithm");
}

} // end namespace

int main(int argc, char** argv) {
  heartbeatbench::initialize(argc, argv);
  std::string infile = deepsea::cmdline::parse_or_default_string("infile", "");
  sptl::ray_cast_test x = sptl::read_from_file<sptl::ray_cast_test>(infile);
  heartbeatbench::benchmark(x);
  return 0;
}
