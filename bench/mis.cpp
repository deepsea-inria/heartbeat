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
#include "sequence.hpp"
#include "speculative_for.hpp"
#include "graph.h"

namespace heartbeatbench {

// **************************************************************
//    MAXIMAL INDEPENDENT SET
// **************************************************************

// For each vertex:
//   Flags = 0 indicates undecided
//   Flags = 1 indicates chosen
//   Flags = 2 indicates a neighbor is chosen
struct MISstep {
  char flag;
  char *Flags;  pbbs::graph::vertex<intT>*G;
  MISstep(char* _F, pbbs::graph::vertex<intT>* _G) : Flags(_F), G(_G) {}
  
  bool reserve(intT i) {
    intT d = G[i].degree;
    flag = 1;
    for (intT j = 0; j < d; j++) {
      intT ngh = G[i].Neighbors[j];
      if (ngh < i) {
        if (Flags[ngh] == 1) { flag = 2; return 1;}
        // need to wait for higher priority neighbor to decide
        else if (Flags[ngh] == 0) flag = 0;
      }
    }
    return 1;
  }
  
  bool commit(intT i) { return (Flags[i] = flag) > 0;}
};

class maximalIndependentSet : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  pbbs::graph::graph<intT> GS; char** dest;
  char* Flags; intT tmp; char z = 0;
  
  maximalIndependentSet() { }
    
  maximalIndependentSet(pbbs::graph::graph<intT> GS, char** dest)
  : GS(GS), dest(dest) { }
  
  heartbeat_dc_declare(heartbeat::edsl, maximalIndependentSet, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        s.Flags = malloc_array<char>(s.GS.n);
        return sequence::fill3(st, pt, s.Flags, s.Flags + s.GS.n, &(s.z));
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        MISstep mis(s.Flags, s.GS.V);
        return heartbeatbench::speculative_for4(st, pt, mis, 0, s.GS.n, 20, &s.tmp);
      }),
      dc::stmt([] (sar& s, par& p) {
        *s.dest = s.Flags;
      })
    });
  }
  
};

heartbeat_pcfg_allocate(maximalIndependentSet, get_cfg)

} // end namespace

#include "sequence.h"
#include "parallel.h"
#include "speculativefor.h"
#include "gettime.h"
//#include "loaders.hpp"
#include "graphio.hpp"

namespace pbbs {

int granularity = 20;

char* maximalIndependentSet(pbbs::graph::graph<intT> GS) {
  intT n = GS.n;
  graph::vertex<intT>* G = GS.V;
  char* Flags = newArray(n, (char) 0);
  heartbeatbench::MISstep mis(Flags, G);
  speculative_for(mis,0,n,granularity);
  return Flags;
}
  
void benchmark(std::string infile) {
  granularity = deepsea::cmdline::parse_or_default_bool("speculative_for_grain", granularity);
  graph::graph<int> g = read_from_file<graph::graph<int>>(infile);
  char* flags = nullptr;
  deepsea::cmdline::dispatcher d;
  d.add("heartbeat", [&] {
    heartbeat::launch_interpreter<heartbeatbench::maximalIndependentSet>(g, &flags);
  });
  d.add("pbbs", [&] {
    heartbeatbench::run_and_report_elapsed_time([&] {
      flags = maximalIndependentSet(g);
    });
  });
  d.dispatch("algorithm");
  if (deepsea::cmdline::parse_or_default_bool("check", false)) {
    auto flags2 = maximalIndependentSet(g);
    for (auto i = 0; i != g.n; i++) {
      if (flags[i] != flags2[i]) {
        char c1 = flags[i];
        char c2 = flags2[i];
        std::cout << "bogus i=" << i << " flags[i]= " << c1 << " flags2[i]= " << c2 << std::endl;
      }
      assert(flags[i] == flags2[i]);
    }
    free(flags2);
  }
  assert(flags != nullptr);
  free(flags);
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
