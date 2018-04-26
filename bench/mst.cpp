/*!
 * \file matching_bench.cpp
 * \brief Benchmarking script for parallel minimal spanning tree
 * \date 2017
 * \copyright COPYRIGHT (c) 2015 Umut Acar, Arthur Chargueraud, and
 * Michael Rainey. All rights reserved.
 * \license This project is released under the GNU Public License.
 *
 */

#include <math.h>
#include <functional>
#include <stdlib.h>

#include "heartbeatbench.hpp"
#include "speculative_for.hpp"
#include "graph.h"
#include "samplesort.hpp"
#include "unionfind.hpp"

namespace heartbeatbench {

struct indexedEdge {
  intT u;
  intT v;
  intT id;
  indexedEdge(intT _u, intT _v, intT _id) : u(_u), v(_v), id(_id) {}
};

struct UnionFindStep {
  intT u;  intT v;  
  indexedEdge *E;  reservation *R;  unionFind UF;  bool *inST;
  UnionFindStep(indexedEdge* _E, unionFind _UF, reservation* _R, bool* ist) 
    : E(_E), R(_R), UF(_UF), inST(ist) {}

  bool reserve(intT i) {
    u = UF.find(E[i].u);
    v = UF.find(E[i].v);
    if (u != v) {
      R[v].reserve(i);
      R[u].reserve(i);
      return 1;
    } else return 0;
  }

  bool commit(intT i) {
    if (R[v].check(i)) {
      R[u].checkReset(i); 
      UF.link(v, u); 
      inST[E[i].id] = 1;
      return 1;}
    else if (R[u].check(i)) {
      UF.link(u, v); 
      inST[E[i].id] = 1;
      return 1; }
    else return 0;
  }
};

template <class E, class F>
class almostKth : public heartbeat::edsl::pcfg::shared_activation_record {
public:

  E* A; E* B; intT k; intT n; F f;
  intT ssize; intT stride; intT km; E* T; E p; bool *flags; intT l;
  intT* dest; _seq<E> tmp;

  almostKth(E* A, E* B, intT k, intT n, F f, intT* dest)
    : A(A), B(B), k(k), n(n), f(f), dest(dest) { }

  heartbeat_private_activation_record_begin(heartbeat::edsl, almostKth, 2)
    int lo; int hi;
  heartbeat_private_activation_record_end(heartbeat::edsl, almostKth, sar, par, dc, get_dc)

  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par& p) {
        s.ssize = std::min<intT>(1000,s.n);
        auto ssize = s.ssize;
        s.stride = s.n/ssize;
        s.km = (intT) (s.k * ((double) ssize) / s.n);
        s.T = malloc_array<E>(ssize);
        auto T = s.T;
        auto A = s.A;
        auto stride = s.stride;
        for (intT i = 0; i < ssize; i++) {
          T[i] = A[i*stride];
        }
        std::sort(T,T+ssize,s.f);
      }),
      dc::stmt([] (sar& s, par& p) {
        s.p = s.T[s.km];
        s.flags = malloc_array<bool>(s.n);
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.n; },
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto flags = s.flags;
        auto A = s.A;
        auto f = s.f;
        auto _p = s.p;
        for (auto i = lo; i != hi; i++) {
          flags[i] = f(A[i],_p);
        }
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        return sequence::pack5(st, pt, s.A,s.B,s.flags,s.n, &s.tmp);
      }),
      dc::stmt([] (sar& s, par& p) {
        s.l = s.tmp.n;
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.n; },
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto flags = s.flags;
        for (auto i = lo; i != hi; i++) {
          flags[i] = ! flags[i];
        }
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        return sequence::pack5(st, pt, s.A,s.B+s.l,s.flags,s.n, &s.tmp);
      }),
      dc::stmt([] (sar& s, par& p) {
        free(s.flags);
        free(s.T);
        *s.dest = s.l;
      })
    });
  }
  
};

template <class E, class F>
typename almostKth<E,F>::cfg_type almostKth<E,F>::cfg = almostKth<E,F>::get_cfg();

typedef std::pair<double,intT> ei;

struct edgeLess {
  bool operator() (ei a, ei b) { 
    return (a.first == b.first) ? (a.second < b.second) 
      : (a.first < b.first);}};

int neg1 = -1;
int zero = 0;

reservation restmp;
  
class mst : public heartbeat::edsl::pcfg::shared_activation_record {
public:

  pbbs::graph::wghEdgeArray<intT> G; std::pair<intT*, intT>* dest;
  pbbs::graph::wghEdge<intT> *E; ei* x;
  intT l; ei* y; unionFind UF; reservation *R; indexedEdge* z;
  bool *mstFlags; bool *flags; intT k; intT* _mst; intT nInMst;
  _seq<intT> tmp; intT tmpi; _seq<ei> tmp22;

  mst(pbbs::graph::wghEdgeArray<intT> G, std::pair<intT*, intT>* dest)
    : G(G), dest(dest) { }

  heartbeat_private_activation_record_begin(heartbeat::edsl, mst, 4)
    int lo; int hi;
  heartbeat_private_activation_record_end(heartbeat::edsl, mst, sar, par, dc, get_dc)

  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par& p) {
        s.E = s.G.E;
        s.x = malloc_array<ei>(s.G.m);
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.G.m; },
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto x = s.x;
        auto E = s.E;
        for (auto i = lo; i != hi; i++) {
          x[i] = ei(E[i].weight,i);
        }
      }),
      dc::stmt([] (sar& s, par& p) {
        s.l = std::min<intT>(4*s.G.n/3,s.G.m);
        s.y = malloc_array<ei>(s.G.m);
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        auto el = edgeLess();
        return heartbeat_call<almostKth<ei,decltype(el)>>(st, pt, s.x, s.y, s.l, s.G.m, el, &s.l);
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        return sampleSort3(st, pt, s.y, s.l, edgeLess());
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        new (&s.UF) unionFind(s.G.n);
        return sequence::fill3(st, pt, s.UF.parents, s.UF.parents + s.G.n, &neg1);
      }),
      dc::stmt([] (sar& s, par& p) {
        s.R = malloc_array<reservation>(s.G.n);
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        return sequence::fill3(st, pt, s.R, s.R + s.G.n, &restmp);
      }),
      dc::stmt([] (sar& s, par& p) {
        s.z = malloc_array<indexedEdge>(s.G.m);
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.l; },
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto y = s.y;
        auto z = s.z;
        auto E = s.E;
        for (auto i = lo; i != hi; i++) {
          intT j = y[i].second;
          z[i] = indexedEdge(E[j].u,E[j].v,j);
        }
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        s.mstFlags = malloc_array<bool>(s.G.m);
        return sequence::fill3(st, pt, s.mstFlags, s.mstFlags + s.G.m, &zero);
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        UnionFindStep UFStep(s.z, s.UF, s.R,  s.mstFlags);
        return speculative_for4(st, pt, UFStep, 0, s.l, 100, &s.tmpi);
      }),
      dc::stmt([] (sar& s, par& p) {
        free(s.z);
        s.flags = malloc_array<bool>(s.G.m-s.l);
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.G.m-s.l; },
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto y = s.y;
        auto& UF = s.UF;
        auto flags = s.flags;
        auto E = s.E;
        auto l = s.l;
        for (auto i = lo; i != hi; i++) {
          intT j = y[i+l].second;
          intT u = UF.find(E[j].u);
          intT v = UF.find(E[j].v);
          if (u != v) flags[i] = 1;
          else flags[i] = 0;
        }
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        return sequence::pack5(st, pt, s.y+s.l, s.x, s.flags, s.G.m-s.l, &s.tmp22);
      }),
      dc::stmt([] (sar& s, par& p) {
        s.k = s.tmp22.n;
        free(s.flags);
        free(s.y);
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        return sampleSort3(st, pt, s.x, s.k, edgeLess());
      }),
      dc::stmt([] (sar& s, par& p) {
        s.z = malloc_array<indexedEdge>(s.k);
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.k; },
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto x = s.x;
        auto E = s.E;
        auto z = s.z;
        for (auto i = lo; i != hi; i++) {
          intT j = x[i].second;
          z[i] = indexedEdge(E[j].u,E[j].v,j);
        }
      }),
      dc::stmt([] (sar& s, par& p) {
        free(s.x);
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        UnionFindStep UFStep(s.z, s.UF, s.R, s.mstFlags);
        return speculative_for4(st, pt, UFStep, 0, s.k, 20, &s.tmpi);
      }),
      dc::stmt([] (sar& s, par& p) {
        free(s.z);
        s._mst = malloc_array<intT>(s.G.m);
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        return sequence::packIndex(st, pt, s._mst, s.mstFlags, s.G.m, &s.tmp);
      }),
      dc::stmt([] (sar& s, par& p) {
        s.nInMst = s.tmp.n;
        free(s.mstFlags);
        s.UF.del();
        free(s.R);
        *s.dest = std::make_pair(s._mst, s.nInMst);
      }),
    });
  }
  
};  

heartbeat_pcfg_allocate(mst, get_cfg)

} // end namespace

#include "readinputbinary.hpp" 
#include "graphio.hpp"
//#include "graphutils.hpp"

#include "mst.h"
#include "gettime.h"
#undef blocked_for

namespace pbbs {
  
graph::wghEdgeArray<int> to_pbbs(sptl::graph::wghEdgeArray<int>& g) {
  graph::wghEdge<int>* e = (graph::wghEdge<int>*) malloc(sizeof(graph::wghEdge<int>) * g.m);
  for (int i = 0; i < g.m; i++) {
    e[i] = graph::wghEdge<int>(g.E[i].u, g.E[i].v, g.E[i].weight);
  }
  return graph::wghEdgeArray<int>(e, g.n, g.m);
}

void benchmark(std::string infile) {
  sptl::graph::graph<int> x = sptl::read_from_file<sptl::graph::graph<int>>(infile);
  sptl::graph::wghEdgeArray<int> edges = to_weighted_edge_array(x);
  auto edges2 = to_pbbs(edges);
  std::pair<intT*, intT> result;
  deepsea::cmdline::dispatcher d;
  d.add("heartbeat", [&] {
    heartbeat::launch_interpreter<heartbeatbench::mst>(edges2, &result);
  });
  d.add("pbbs", [&] {
    heartbeatbench::run_and_report_elapsed_time([&] {
      result = mst(edges2);
    });
  });
  d.dispatch("algorithm");
  if (deepsea::cmdline::parse_or_default_bool("check", false)) {
    auto result2 = mst(edges2);
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
  free(result.first);
}

} // end namespace

int main(int argc, char** argv) {
  heartbeatbench::initialize(argc, argv);
  std::string infile = deepsea::cmdline::parse_or_default_string("infile", "");
  pbbs::benchmark(infile);
  return 0;
}
