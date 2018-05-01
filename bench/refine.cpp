
#include "heartbeatbench.hpp"
//#include "geometrydata.hpp"
#include "refine.hpp"
#include "readinputbinary.hpp" 
#include "refine.h"

namespace sptl {

pbbs::triangles<pbbs::_point2d<double>> to_pbbs(sptl::triangles<sptl::_point2d<double>>& x) {
  pbbs::triangles<pbbs::_point2d<double>> result;
  result.numPoints = x.num_points;
  result.P = (pbbs::_point2d<double>*) malloc(sizeof(pbbs::_point2d<double>) * x.num_points);
  for (int i = 0; i < x.num_points; i++) {
    result.P[i] = pbbs::_point2d<double>(x.p[i].x, x.p[i].y);
  }
  result.numTriangles = x.num_triangles;
  result.T = (pbbs::triangle*) malloc(sizeof(pbbs::triangle) * x.num_triangles);
  for (int i = 0; i < x.num_triangles; i++) {
    result.T[i] = pbbs::triangle(x.t[i].vertices[0], x.t[i].vertices[1], x.t[i].vertices[2]);
  }
  return result;
}

void benchmark() {
  std::string infile = deepsea::cmdline::parse_or_default_string("infile", "");
  auto x = sptl::read_from_file<sptl::triangles<sptl::_point2d<double>>>(infile);
  std::string algorithm = deepsea::cmdline::parse<std::string>("algorithm");
  deepsea::cmdline::dispatcher d;
  heartbeatbench::triangles<_point2d<double>> res;
  d.add("heartbeat", [&] {
    heartbeat::launch_interpreter<heartbeatbench::refine>(x, &res);
  });
  d.add("pbbs", [&] {
    auto y = to_pbbs(x);
    heartbeatbench::run_and_report_elapsed_time([&] {
      pbbs::refine(y);
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


