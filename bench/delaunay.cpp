
#include <math.h>
#include <functional>
#include <stdlib.h>

#include "heartbeatbench.hpp"
//#include "geometrydata.hpp"
#include "readinputbinary.hpp" 
#include "delaunay.hpp"
#undef blocked_for
#include "delaunay.h"

namespace sptl {

parray<pbbs::_point2d<double>> to_pbbs(parray<sptl::_point2d<double>>& points) {
  parray<pbbs::_point2d<double>> result(points.size());
  for (int i = 0; i < points.size(); i++) {
    result[i] = pbbs::_point2d<double>(points[i].x, points[i].y);
  }
  return result;
} 

void benchmark() {
  std::string infile = deepsea::cmdline::parse_or_default_string("infile", "");
  parray<sptl::_point2d<double>> x = sptl::read_from_file<parray<sptl::_point2d<double>>>(infile);
  std::string algorithm = deepsea::cmdline::parse<std::string>("algorithm");
  deepsea::cmdline::dispatcher d;
  heartbeatbench::triangles<_point2d<double>> res;
  int n = x.size();
  d.add("heartbeat", [&] {
    heartbeat::launch_interpreter<heartbeatbench::delaunay>(&x[0], n, &res);
  });
  d.add("pbbs", [&] {
    auto y = to_pbbs(x);
    heartbeatbench::run_and_report_elapsed_time([&] {
      pbbs::delaunay(&y[0], (int)y.size());
    });
  });
  d.dispatch("algorithm");
}
  
} // end namespace

int main(int argc, char** argv) {
  heartbeatbench::initialize(argc, argv);
  sptl::benchmark();
  return 0;
}


