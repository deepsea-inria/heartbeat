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
#include "bfs.hpp"
#include "bfs.h"
#include "graphio.hpp"

namespace pbbs {

void benchmark(std::string infile) {
  graph::graph<int> g = read_from_file<graph::graph<int>>(infile);
  intT source = deepsea::cmdline::parse_or_default("source", 0);
  std::pair<intT,intT> result;
  deepsea::cmdline::dispatcher d;
  d.add("heartbeat", [&] {
    heartbeat::launch_interpreter<heartbeatbench::bfs>(source, g, &result);
  });
  d.add("pbbs", [&] {
    heartbeatbench::run_and_report_elapsed_time([&] {
      result = BFS(source, g);
    });
  });
  d.dispatch("algorithm");
  std::cout << "nb_visited " << result.first << std::endl;
  std::cout << "bfs_max_dist " << result.second << std::endl;
}

} // end namespace

int main(int argc, char** argv) {
  heartbeatbench::initialize(argc, argv);
  std::string infile = deepsea::cmdline::parse_or_default_string("infile", "");
  if (infile == "") {
    std::cout << "bogus input file name" << std::endl;
    exit(0);
  }
  pbbs::benchmark(infile);
  return 0;
}
