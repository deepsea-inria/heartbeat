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
#include <limits.h>
#include "sequence.h"
#include "parallel.h"
#include "graph.h"
#include "spanning.h"
#include "speculativefor.h"
#include "union.h"
#include "unionfind.hpp"

#ifndef SPANNING_INCLUDED_H
#define SPANNING_INCLUDED_H

namespace heartbeatbench {

using namespace std;

template <class T>
using edge = pbbs::graph::edge<T>;

template <class T>
using edgeArray = pbbs::graph::edgeArray<T>;

template <class intT>
struct unionFindStep {
  intT u;  intT v;  
  edge<intT> *E;  reservation *R;  unionFind UF;
  unionFindStep() { }
  unionFindStep(edge<intT>* _E, unionFind _UF, reservation* _R)
    : E(_E), R(_R), UF(_UF) {} 

  bool reserve(intT i) {
    u = UF.find(E[i].u);
    v = UF.find(E[i].v);
    if (u > v) {intT tmp = u; u = v; v = tmp;}
    if (u != v) {
      R[v].reserve(i);
      return 1;
    } else return 0;
  }

  bool commit(intT i) {
    if (R[v].check(i)) { UF.link(v, u); return 1; }
    else return 0;
  }
};

struct notMax { bool operator() (intT i) {return i < INT_T_MAX;}};

int neg1 = -1;
  
class st : public heartbeat::edsl::pcfg::shared_activation_record {
public:

  edgeArray<intT> G; pair<intT*,intT>* dest; 
  intT m; intT n; unionFind UF; reservation *R; unionFindStep<intT> UFStep;
  intT tmpi; _seq<intT> stIdx; int grain;

  st(edgeArray<intT> G, int grain, pair<intT*,intT>* dest)
    : G(G), grain(grain), dest(dest) { }

  heartbeat_private_activation_record_begin(heartbeat::edsl, st, 1)
    int lo; int hi;
  heartbeat_private_activation_record_end(heartbeat::edsl, st, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par& p) {
        s.m = s.G.nonZeros;
        s.n = s.G.numRows;
        new (&s.UF) unionFind(s.n);
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        return sequence::fill3(st, pt, s.UF.parents, s.UF.parents + s.n, &neg1);
      }),
      dc::stmt([] (sar& s, par& p) {
        s.R = malloc_array<reservation>(s.n);
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.n; },
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto R = s.R;
        for (auto i = lo; i != hi; i++) {
          new ((void*) (R+i)) reservation;
        }
      }),
      dc::stmt([] (sar& s, par& p) {
        new (&s.UFStep) unionFindStep<intT>(s.G.E, s.UF, s.R); 
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        return speculative_for4(st, pt, s.UFStep, 0, s.m, s.grain, &s.tmpi);
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        return sequence::filter4(st, pt, (intT*) s.R, s.n, notMax(), &s.stIdx);
      }),
      dc::stmt([] (sar& s, par& p) {
        s.UF.del(); free(s.R);
        *s.dest = std::make_pair(s.stIdx.A, s.stIdx.n);
      })
    });
  }

};

heartbeat_pcfg_allocate(st, get_cfg)

} // end namespace

#endif
