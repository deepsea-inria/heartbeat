// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include <iostream>
#include <chrono>
#include <array>
#include <fstream>
#include <iostream>
#include <istream>
#include <ostream>
#include <sstream>
#include <algorithm>

#include "heartbeatbench.hpp"
#include "speculative_for.hpp"
#include "graph.h"
#include "readinputbinary.hpp" 
#include "graphio.hpp"
//#include "graphutils.hpp"

#include "spanning.h"
#include "gettime.h"
#undef blocked_for
#include "spanning.hpp"

namespace pbbs {

pbbs::graph::edgeArray<int> to_pbbs(sptl::graph::edgeArray<int>& g) {
  pbbs::graph::edge<int>* e = (pbbs::graph::edge<int>*) malloc(sizeof(pbbs::graph::edge<int>) * g.nonZeros);
  for (int i = 0; i < g.nonZeros; i++) {
    e[i] = pbbs::graph::edge<int>(g.E[i].u, g.E[i].v);
  }
  return pbbs::graph::edgeArray<int>(e, g.numRows, g.numCols, g.nonZeros);
}
  
void benchmark(std::string infile) {
  int grain = deepsea::cmdline::parse_or_default_int("speculative_for_grain", 50);
  sptl::graph::graph<int> x = sptl::read_from_file<sptl::graph::graph<int>>(infile);  
  sptl::graph::edgeArray<int> edges = sptl::graph::to_edge_array(x);
  auto edges2 = to_pbbs(edges);
  std::pair<intT*, intT> result;
  deepsea::cmdline::dispatcher d;
  d.add("heartbeat", [&] {
    heartbeat::launch_interpreter<heartbeatbench::st>(edges2, grain, &result);
  });
  d.add("pbbs", [&] {
    heartbeatbench::run_and_report_elapsed_time([&] {
      result = spanningTree(edges2);
    });
  });
  d.dispatch("algorithm");
  if (deepsea::cmdline::parse_or_default_bool("check", false)) {
    auto result2 = spanningTree(edges2);
    assert(result.second == result2.second);
    auto n = result2.second;
    auto res1 = result.first;
    auto res2 = result2.first;
    for (auto i = 0; i != n; i++) {
      if (res1[i] != res2[i]) {
        auto c1 = res1[i];
        auto c2 = res2[i];
        std::cout << "bogus i=" << i << " result1[i]= " << c1 << " result2[i]= " << c2 << std::endl;
      }
      assert(res1[i] == res2[i]);
    }
    free(res2);
  }
  assert(result.first != nullptr);
  free(result.first);
}

} // end namespace

int main(int argc, char** argv) {
  heartbeatbench::initialize(argc, argv);
  sequence::initialize();
  std::string infile = deepsea::cmdline::parse_or_default_string("infile", "");
  if (infile == "") {
    std::cout << "bogus input file name" << std::endl;
    exit(0);
  }
  pbbs::benchmark(infile);
  return 0;
}
