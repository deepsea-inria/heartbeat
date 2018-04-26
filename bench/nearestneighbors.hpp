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

// This is a parallel version of the algorithm described in
//  Juha Karkkainen and Peter Sanders.
//  Simple linear work suffix array construction.
//  Proc. ICALP 2003.  pp 943
// It includes code for finding the LCP
//   Written by Guy Blelloch and Julian Shun

#include <iostream>
#include <limits>

#include "sequence.hpp"
#include "blockradixsort.hpp"
#include "utils.h"
#include "heartbeat.hpp"
#include "octtree.hpp"
#include "gettime.h"

#ifndef _HEARTBEAT_NEARESTNEIGHBORS_H_
#define _HEARTBEAT_NEARESTNEIGHBORS_H_

namespace heartbeatbench {

// A k-nearest neighbor structure
// requires vertexT to have pointT and vectT typedefs
template <class vertexT, int maxK>
struct kNearestNeighbor {
  typedef vertexT vertex;
  typedef typename vertexT::pointT point;
  typedef typename point::vectT fvect;

  typedef gTreeNode<point,fvect,vertex> qoTree;
  typedef gTreeNodeConstructor2<point,fvect,vertex, typename qoTree::ndty> qoTreeCons;
  qoTree *tree;

  void del() {tree->del();}

  struct kNN {
    vertex *ps;  // the element for which we are trying to find a NN
    vertex *pn[maxK];  // the current k nearest neighbors (nearest last)
    double rn[maxK]; // radius of current k nearest neighbors
    int quads;
    int k;
    kNN() {}

    // returns the ith smallest element (0 is smallest) up to k-1
    vertex* operator[] (const int i) { return pn[k-i-1]; }

    kNN(vertex *p, int kk) {
      if (kk > maxK) {std::cout << "k too large in kNN" << std::endl; abort();}
      k = kk;
      quads = (1 << (p->pt).dimension());
      ps = p;
      for (int i=0; i<k; i++) {
        pn[i] = (vertex*) NULL; 
        rn[i] = std::numeric_limits<double>::max();
      }
    }

    // if p is closer than pn then swap it in
    void update(vertex *p) { 
      //inter++;
      point opt = (p->pt);
      fvect v = (ps->pt) - opt;
      double r = v.length();
      if (r < rn[0]) {
        pn[0]=p; rn[0] = r;
        for (int i=1; i < k && rn[i-1]<rn[i]; i++) {
          std::swap(rn[i-1],rn[i]); std::swap(pn[i-1],pn[i]); }
      }
    }

    // looks for nearest neighbors in boxes for which ps is not in
    void nearestNghTrim(qoTree *T) {
      if (!(T->center).out_of_box(ps->pt, (T->size/2)+rn[0])) {
        if (T->IsLeaf())
          for (int i = 0; i < T->count; i++) update(T->vertices[i]);
        else 
          for (int j=0; j < quads; j++) nearestNghTrim(T->children[j]);
      }
    }

    // looks for nearest neighbors in box for which ps is in
    void nearestNgh(qoTree *T) {
      if (T->IsLeaf())
        for (int i = 0; i < T->count; i++) {
          vertex *pb = T->vertices[i];
          if (pb != ps) update(pb);
        }
      else {
        int i = T->findQuadrant(ps);
        nearestNgh(T->children[i]);
        for (int j=0; j < quads; j++) 
          if (j != i) nearestNghTrim(T->children[j]);
      }
    }
  };

  vertex* nearest(vertex *p) {
    kNN nn(p,1);
    nn.nearestNgh(tree); 
    return nn[0];
  }

  // version that writes into result
  void kNearest(vertex *p, vertex** result, int k) {
    kNN nn(p,k);
    nn.nearestNgh(tree); 
    for (int i=0; i < k; i++) result[i] = 0;
    for (int i=0; i < k; i++) result[i] = nn[i];
  }

  // version that allocates result
  vertex** kNearest(vertex *p, int k) {
    vertex** result = newA(vertex*,k);
    kNearest(p,result,k);
    return result;
  }

};

// find the k nearest neighbors for all points in tree
// places pointers to them in the .ngh field of each vertex
template <int maxK, class vertexT>
class ANN : public heartbeat::edsl::pcfg::shared_activation_record {
public:

  typedef kNearestNeighbor<vertexT,maxK> kNNT;

  vertexT** v; int n; int k;
  kNNT T; vertexT** vr;

  ANN(vertexT** v, int n, int k)
    : v(v), n(n), k(k) { } 
  
  heartbeat_private_activation_record_begin(heartbeat::edsl, ANN, 1)
    int s; int e;
  heartbeat_private_activation_record_end(heartbeat::edsl, ANN, sar, par, dc, get_dc)

  static
  dc get_dc() {
    return dc::stmts({
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        using gtc = typename kNNT::qoTreeCons;
        new (&s.T) kNNT;
        return heartbeat_call<gtc>(st, pt, s.v, s.n, &(s.T.tree));
      }),
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        // this reorders the vertices for locality          
        return flatten(st, pt, s.T.tree, &(s.vr));
      }),
      // find nearest k neighbors for each point
      dc::parallel_for_loop([] (sar& s, par& p) { p.s = 0; p.e = s.n;},
                            [] (par& p) { return std::make_pair(&p.s, &p.e); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto& T = s.T;
        auto vr = s.vr;
        auto k = s.k;
        for (auto i = lo; i != hi; i++) {
          T.kNearest(vr[i], vr[i]->ngh, k);
        }
      }),
      dc::stmt([] (sar& s, par& p) {
        delete [] s.vr;
        s.T.del();
      })
    });
  }
  
};

template <int maxK, class vertexT>
typename ANN<maxK,vertexT>::cfg_type ANN<maxK,vertexT>::cfg = ANN<maxK,vertexT>::get_cfg();

template <class PT, int KK>
struct vertexNN {
  typedef PT pointT;
  int identifier;
  pointT pt;         // the point itself
  vertexNN* ngh[KK];    // the list of neighbors
  vertexNN(pointT p, int id) : pt(p), identifier(id) {}
};

  	pbbs::timer fnn1;
  
template <int maxK, class pointT>
class findNearestNeighbors : public heartbeat::edsl::pcfg::shared_activation_record {
public:

  typedef vertexNN<pointT, maxK> vertex;
  
  pointT* p; int n; int k;
  vertex** v; vertex* vv; intT* result;

  findNearestNeighbors(pointT* p, int n, int k, intT* result)
    : p(p), n(n), k(k), result(result) { } 
  
  heartbeat_private_activation_record_begin(heartbeat::edsl, findNearestNeighbors, 2)
    int s; int e;
  heartbeat_private_activation_record_end(heartbeat::edsl, findNearestNeighbors, sar, par, dc, get_dc)

  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par& p) {
        s.v = malloc_array<vertex*>(s.n);
        s.vv = malloc_array<vertex>(s.n);
	fnn1.start();
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.s = 0; p.e = s.n; },
                            [] (par& p) { return std::make_pair(&p.s, &p.e); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto v = s.v;
        auto vv = s.vv;
        auto _p = s.p;
        for (auto i = lo; i != hi; i++) {
          v[i] = new (&vv[i]) vertex(_p[i],i);
        }
      }),
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
	  fnn1.stop();
        return heartbeat_call<ANN<maxK, vertex>>(st, pt, s.v, s.n, s.k);
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.s = 0; p.e = s.n; },
                            [] (par& p) { return std::make_pair(&p.s, &p.e); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto result = s.result;
        auto v = s.v;
        auto k = s.k;
        auto n = s.n;
        for (auto i = lo; i != hi; i++) {
          for (int j = 0; j < std::min(k, n - 1); j++) {
            result[i * k + j] = v[i]->ngh[j]->identifier;
          }
        }
      })
    });
  }
  
};

template <int maxK, class pointT>
typename findNearestNeighbors<maxK,pointT>::cfg_type findNearestNeighbors<maxK,pointT>::cfg = findNearestNeighbors<maxK,pointT>::get_cfg();

template <int maxK, class pointT>
stack_type findNearestNeighbors3(stack_type st, plt_type pt, pointT* p, int n, int k, intT** result) {
  *result = malloc_array<intT>(k * n);
  return sequence::heartbeat_call<findNearestNeighbors<maxK,pointT>>(st, pt, p, n, k, *result);
}

  
} // end namespace

#endif
