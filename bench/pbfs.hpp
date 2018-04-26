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

#include "utils.h"
#include "graph.h"
#include "sequence.hpp"

#ifndef _HEARTBEAT_PBFS_H_
#define _HEARTBEAT_PBFS_H_

namespace heartbeatbench {

// **************************************************************
//    Non-DETERMINISTIC BREADTH FIRST SEARCH
// **************************************************************

// **************************************************************
//    THE NON-DETERMINISTIC BSF
//    Updates the graph so that it is the BFS tree (i.e. the neighbors
//      in the new graph are the children in the bfs tree)
// **************************************************************

struct nonNegF{bool operator() (intT a) {return (a>=0);}};

int zero = 0;

class pbfs : public heartbeat::edsl::pcfg::shared_activation_record {
public:

  intT start; pbbs::graph::graph<intT> GA; pbbs::graph::vertex<intT>* G;
  std::pair<intT,intT>* result;

  intT numVertices; intT numEdges; intT* Frontier; intT* Visited;
  intT* FrontierNext; intT* Counts; intT totalVisited; intT round;
  intT frontierSize; intT nr; intT v; intT o;

  pbfs(intT start, pbbs::graph::graph<intT> GA, std::pair<intT,intT>* result) :
    start(start), GA(GA), result(result) { }

  heartbeat_private_activation_record_begin(heartbeat::edsl, pbfs, 3)
    int s2; int e2;
    int s3; int e3;
    int s4; int e4;
  heartbeat_private_activation_record_end(heartbeat::edsl, pbfs, sar, par, dc, get_dc)

  static
  dc get_dc() {
    return dc::stmts({
       dc::stmt([] (sar& s, par& p) {
         s.numVertices = s.GA.n;
         s.numEdges = s.GA.m;
         s.G = s.GA.V;
         s.Frontier = malloc_array<intT>(s.numEdges);
         s.Visited = malloc_array<intT>(s.numVertices);
         s.FrontierNext = malloc_array<intT>(s.numEdges);
         s.Counts = malloc_array<intT>(s.numVertices);
       }),
       dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
         return sequence::fill3(st, pt, s.Visited, s.Visited + s.numVertices, &zero);
       }),
       dc::stmt([] (sar& s, par& p) {
         s.Frontier[0] = s.start;
         s.frontierSize = 1;
         s.Visited[s.start] = 1;
         s.totalVisited = 0;
         s.round = 0;
       }),
       dc::sequential_loop([] (sar& s, par&) { return s.frontierSize > 0; }, dc::stmts({
         dc::stmt([] (sar& s, par& p) {
           s.round++;
           s.totalVisited += s.frontierSize;
         }),
         dc::parallel_for_loop([] (sar& s, par& p) { p.s2 = 0; p.e2 = s.frontierSize; },
                               [] (par& p) { return std::make_pair(&p.s2, &p.e2); },
                               [] (sar& s, par& p, int lo, int hi) {
           auto Counts = s.Counts;
           auto G = s.G;
           auto Frontier = s.Frontier;
           for (auto i = lo; i != hi; i++) {
             Counts[i] = G[Frontier[i]].degree;
           }
         }),
         dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
           return sequence::scan6(st, pt, s.Counts, s.Counts, s.frontierSize, pbbs::utils::addF<intT>(),(intT)0, &(s.nr));
         }),
         dc::stmt([] (sar& s, par& p) {
           p.s3 = 0;
           p.e3 = s.frontierSize;
         }),
         dc::parallel_for_loop([] (par& p) { return std::make_pair(&p.s3, &p.e3); }, dc::stmts({
           dc::parallel_for_loop([] (sar& s, par& p) {
             auto i = p.s3;
             auto v = s.v = s.Frontier[i];
             auto o = s.o = s.Counts[i];
             p.s4 = 0;
             p.e4 = s.G[v].degree;
           }, [] (par& p) { return std::make_pair(&p.s4, &p.e4);
           }, [] (sar& s, par& p, int lo2, int hi2) {
             auto Frontier = s.Frontier;
             auto Counts = s.Counts;
             auto G = s.G;
             auto Visited = s.Visited;
             auto FrontierNext = s.FrontierNext;
             auto v = s.v;
             auto o = s.o;
             for (auto j = lo2; j != hi2; j++) {
               intT ngh = G[v].Neighbors[j];
               if (Visited[ngh] == 0 && !__sync_val_compare_and_swap(&Visited[ngh], 0, 1)) {
                 FrontierNext[o+j] = ngh;
               } else {
                 FrontierNext[o+j] = -1;
               }
             }
           }),
           dc::stmt([] (sar& s, par& p) {
             p.s3++;
           })
         })),
         // Filter out the empty slots (marked with -1)
         dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
           return sequence::filterDPS5(st, pt, s.FrontierNext, s.Frontier, s.nr, nonNegF(), &(s.frontierSize));
         }),
       })),
       dc::stmt([] (sar& s, par& p) {
         *s.result = std::make_pair(s.totalVisited, s.round);
         free(s.FrontierNext); free(s.Frontier); free(s.Counts); free(s.Visited);
       })
     });
   }
  
};

heartbeat_pcfg_allocate(pbfs, get_cfg)
  
} // end namespace

#endif
