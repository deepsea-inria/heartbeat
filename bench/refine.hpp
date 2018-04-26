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
#include <vector>

#include "atomic.hpp"
#include "sequence.hpp"
#include "utils.h"
#include "logging.hpp"
#include "topology.hpp"
#include "deterministichash.hpp"
#include "logging.hpp"

#ifndef _HEARTBEAT_REFINE_H_
#define _HEARTBEAT_REFINE_H_

namespace heartbeatbench {

using namespace std;

  namespace logging = heartbeat::logging;

using vect2d = sptl::vect2d;

// *************************************************************
//   PARALLEL HASH TABLE TO STORE WORK QUEUE OF SKINNY TRIANGLES
// *************************************************************

struct hashTriangles {
  typedef tri* eType;
  typedef tri* kType;
  eType empty() {return NULL;}
  kType getKey(eType v) { return v;}
  uintT hash(kType s) { return pbbs::utils::hash(s->id); }
  int cmp(kType s, kType s2) {
    return (s->id > s2->id) ? 1 : ((s->id == s2->id) ? 0 : -1);
  }
  bool replaceQ(eType s, eType s2) {return 0;}
};

typedef Table<hashTriangles,intT> TriangleTable;
  //TriangleTable makeTriangleTable(intT m) {return TriangleTable(m,hashTriangles());}

// *************************************************************
//   THESE ARE TAKEN FROM delaunay.C
//   Perhaps should be #included
// *************************************************************

// Holds vertex and simplex queues used to store the cavity created 
// while searching from a vertex between when it is initially searched 
// and later checked to see if all corners are reserved.
struct Qs {
  vector<vertex*> vertexQ;
  vector<simplex> simplexQ;
};

// Recursive routine for finding a cavity across an edge with
// respect to a vertex p.
// The simplex has orientation facing the direction it is entered.
//
//         a
//         | \ --> recursive call
//   p --> |T c 
// enter   | / --> recursive call
//         b
//
//  If p is in circumcircle of T then 
//     add T to simplexQ, c to vertexQ, and recurse
void findCavity(simplex t, vertex *p, Qs *q) {
  if (t.inCirc(p)) {
    q->simplexQ.push_back(t);
    t = t.rotClockwise();
    findCavity(t.across(), p, q);
    q->vertexQ.push_back(t.firstVertex());
    t = t.rotClockwise();
    findCavity(t.across(), p, q);
  }
}

// Finds the cavity for v and tries to reserve vertices on the 
// boundary (v must be inside of the simplex t)
// The boundary vertices are pushed onto q->vertexQ and
// simplices to be deleted on q->simplexQ (both initially empty)
// It makes no side effects to the mesh other than to X->reserve
void reserveForInsert(vertex *v, simplex t, Qs *q) {
  // each iteration searches out from one edge of the triangle
  for (int i=0; i < 3; i++) {
    q->vertexQ.push_back(t.firstVertex());
    findCavity(t.across(), v, q);
    t = t.rotClockwise();
  }
  // the maximum id new vertex that tries to reserve a boundary vertex 
  // will have its id written.  reserve starts out as -1
  for (intT i = 0; i < q->vertexQ.size(); i++)
    pbbs::utils::writeMax(&((q->vertexQ)[i]->reserve), v->id);
}

// *************************************************************
//   DEALING WITH THE CAVITY
// *************************************************************

inline bool skinnyTriangle(tri *t) {
  double minAngle = 30;
  if (minAngleCheck(t->vtx[0]->pt, t->vtx[1]->pt, t->vtx[2]->pt, minAngle))
    return 1;
  return 0;
}

inline bool obtuse(simplex t) {
  int o = t.o;
  point2d p0 = t.t->vtx[(o+1)%3]->pt;
  vect2d v1 = t.t->vtx[o]->pt - p0;
  vect2d v2 = t.t->vtx[(o+2)%3]->pt - p0;
  return (v1.dot(v2) < 0.0);
}

inline point2d circumcenter(simplex t) {
  if (t.isTriangle())
    return triangleCircumcenter(t.t->vtx[0]->pt, t.t->vtx[1]->pt, t.t->vtx[2]->pt);
  else { // t.isBoundary()
    point2d p0 = t.t->vtx[(t.o+2)%3]->pt;
    point2d p1 = t.t->vtx[t.o]->pt;
    return p0 + (p1-p0)/2.0;
  }
}

// this side affects the simplex by moving it into the right orientation
// and setting the boundary if the circumcenter encroaches on a boundary
inline bool checkEncroached(simplex& t) {
  if (t.isBoundary()) return 0;
  int i;
  for (i=0; i < 3; i++) {
    if (t.across().isBoundary() && (t.farAngle() > 45.0)) break;
    t = t.rotClockwise();
  }
  if (i < 3) return t.boundary = 1;
  else return 0;
}

bool findAndReserveCavity(vertex* v, simplex& t, Qs* q) {
  t = simplex(v->badT,0);
  if (t.t == NULL) {cout << "refine: nothing in badT" << endl; abort();}
  if (t.t->bad == 0) return 0;
  //t.t->bad = 0;

  // if there is an obtuse angle then move across to opposite triangle, repeat
  if (obtuse(t)) t = t.across();
  while (t.isTriangle()) {
    int i;
    for (i=0; i < 2; i++) {
      t = t.rotClockwise();
      if (obtuse(t)) { t = t.across(); break; } 
    }
    if (i==2) break;
  }

  // if encroaching on boundary, move to boundary
  checkEncroached(t);

  // use circumcenter to add (if it is a boundary then its middle)
  v->pt = circumcenter(t);
  reserveForInsert(v, t, q);
  return 1;
}

// checks if v "won" on all adjacent vertices and inserts point if so
// returns true if "won" and cavity was updated
bool addCavity(vertex *v, simplex t, Qs *q, TriangleTable TT) {
  bool flag = 1;
  for (intT i = 0; i < q->vertexQ.size(); i++) {
    vertex* u = (q->vertexQ)[i];
    if (u->reserve == v->id) u->reserve = -1; // reset to -1
    else flag = 0; // someone else with higher priority reserved u
  }
  if (flag) {
    tri* t0 = t.t;
    tri* t1 = v->t;  // the memory for the two new triangles
    tri* t2 = t1 + 1;  
    t1->initialized = 1;
    if (t.isBoundary()) t.splitBoundary(v, t1);
    else {
      t2->initialized = 1;
      t.split(v, t1, t2);
    }
    
    // update the cavity
    for (intT i = 0; i<q->simplexQ.size(); i++) 
      (q->simplexQ)[i].flip();
    q->simplexQ.push_back(simplex(t0,0));
    q->simplexQ.push_back(simplex(t1,0));
    if (!t.isBoundary()) q->simplexQ.push_back(simplex(t2,0));
    
    for (intT i = 0; i<q->simplexQ.size(); i++) {
      tri* t = (q->simplexQ)[i].t;
      if (skinnyTriangle(t)) {
        TT.insert(t); 
        t->bad = 1;}
      else t->bad = 0;
    }
    v->badT = NULL;
  } 
  q->simplexQ.clear();
  q->vertexQ.clear();
  return flag;
}

// *************************************************************
//    MAIN REFINEMENT LOOP
// *************************************************************

class addRefiningVertices : public heartbeat::edsl::pcfg::shared_activation_record {
public:

  vertex** v; intT n; intT nTotal; TriangleTable TT; intT* dest;
  intT maxR; Qs *qqs; Qs **qs; simplex *t; bool *flags; vertex** h;
  intT cnt; vertex** vv; intT k; _seq<vertex*> tmp; int lo; int hi;
 
  intT top; intT failed;  intT size;
  
  addRefiningVertices(vertex** v, intT n, intT nTotal, TriangleTable TT, intT* dest)
  : v(v), n(n), nTotal(nTotal), TT(TT), dest(dest) { }
  
  heartbeat_private_activation_record_begin(heartbeat::edsl, addRefiningVertices, 2)
    int lo; int hi;
  heartbeat_private_activation_record_end(heartbeat::edsl, addRefiningVertices, sar, par, dc, get_dc)

  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par& p) { 
        auto maxR = s.maxR = (intT) (s.nTotal/500) + 1; // maximum number to try in parallel
        s.qqs = malloc_array<Qs>(maxR);
        s.qs = malloc_array<Qs*>(maxR);
      }),
      dc::mk_if([] (sar& s, par&) { return s.maxR < 500; },
        dc::stmt([] (sar& s, par& p) {
          auto maxR = s.maxR;
          auto qqs = s.qqs;
          auto qs = s.qs;
          for (intT i=0; i < maxR; i++) qs[i] = new (&qqs[i]) Qs;
        }),
        dc::sequential_loop([] (sar& s, par&) { s.lo = 0; s.hi = s.maxR; },
                            [] (sar& s, par&) { return std::make_pair(&s.lo, &s.hi); },
                            [] (sar& s, par&, int lo, int hi) {
          auto qs = s.qs;
          auto qqs = s.qqs;
          for (int i = lo; i != hi; i++) {
            qs[i] = new (&qqs[i]) Qs;
          }
        })),
      dc::stmt([] (sar& s, par& p) { 
        s.t = malloc_array<simplex>(s.maxR);
        s.flags = malloc_array<bool>(s.maxR);
        s.h = malloc_array<vertex*>(s.maxR);
 
        s.top = s.n; s.failed = 0;
        s.size = s.maxR;
      }),
      dc::sequential_loop([] (sar& s, par&) { return s.top > 0; }, dc::stmts({
        dc::stmt([] (sar& s, par& p) {
          s.cnt = min<intT>(s.size,s.top);
          s.vv = s.v+s.top-s.cnt;
        }),
        dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.cnt;},
                              [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                              [] (sar& s, par& p, int lo, int hi) {
          auto flags = s.flags;
          auto vv = s.vv;
          auto t = s.t;
          auto qs = s.qs;
          for (auto j = lo; j != hi; j++) {
            flags[j] = findAndReserveCavity(vv[j],t[j],qs[j]);
          }
        }),
        dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.cnt;},
                              [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                              [] (sar& s, par& p, int lo, int hi) {
          auto flags = s.flags;
          auto vv = s.vv;
          auto t = s.t;
          auto qs = s.qs;
          auto TT = s.TT;
          for (auto j = lo; j != hi; j++) {
            flags[j] = flags[j] && !addCavity(vv[j], t[j], qs[j], TT);
          }
        }),
        dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
          return sequence::pack5(st, pt, s.vv, s.h, s.flags, s.cnt, &s.tmp);
        }),
        dc::stmt([] (sar& s, par& p) {
          s.k = s.tmp.n;
        }),
        dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
          return sequence::copy3(st, pt, s.h, s.h + s.k, s.vv);
        }),
        dc::stmt([] (sar& s, par& p) {
          s.failed += s.k;
          s.top = s.top-s.cnt+s.k;            
        })
      })),
      dc::stmt([] (sar& s, par& p) {
        free(s.qqs); free(s.qs); free(s.t);
        free(s.flags); free(s.h); 
        *s.dest = s.failed;
      }) 
    });
  }

};

heartbeat_pcfg_allocate(addRefiningVertices, get_cfg)

// *************************************************************
//    DRIVER
// *************************************************************

class refine : public heartbeat::edsl::pcfg::shared_activation_record {
public:

  triangles<point2d> Tri; triangles<point2d>* dest;
  int expandFactor = 4;
  intT n; intT m; intT extraVertices; intT totalVertices; intT totalTriangles;
  intT numBad; intT num_triangs; intT num_points;
  bool* flags; _seq<intT> I; point2d* rp; triangle* rt;
  intT nO; intT* II;

  vertex** v; vertex* vv; tri* Triangs;
  TriangleTable workQ; _seq<tri*> badTT; _seq<tri*> badT;
  typedef typename hashTriangles::eType eType;
  eType empty; intT tmpi;

  refine(triangles<point2d> Tri, triangles<point2d>* dest)
    : Tri(Tri), dest(dest) { }

  heartbeat_private_activation_record_begin(heartbeat::edsl, refine, 10)
    int lo; int hi;
  heartbeat_private_activation_record_end(heartbeat::edsl, refine, sar, par, dc, get_dc)

  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par& p) {
        s.n = s.Tri.num_points;
        s.m = s.Tri.num_triangles;
        s.extraVertices = s.expandFactor*s.n;
        s.totalVertices = s.n + s.extraVertices;
        s.totalTriangles = s.m + 2 * s.extraVertices;

        s.v = malloc_array<vertex*>(s.extraVertices);
        s.vv = malloc_array<vertex>(s.totalVertices);
        s.Triangs = malloc_array<tri>(s.totalTriangles);
        new (&s.empty) eType(hashTriangles().empty());
      }),
      dc::stmt([] (sar& s, par& p) {
        logging::push_event(logging::algo_phase);
      }),
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        return heartbeat_call<topologyFromTriangles>(st, pt, s.Tri, &s.vv, &s.Triangs);
      }),
      dc::stmt([] (sar& s, par& p) {
        logging::push_event(logging::algo_phase);
      }),
      //  set up extra triangles
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = s.m; p.hi = s.totalTriangles;},
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto Triangs = s.Triangs;
        for (auto i = lo; i != hi; i++) {
          Triangs[i].id = i;
          Triangs[i].initialized = 0;
        }
      }),
      dc::stmt([] (sar& s, par& p) {
        logging::push_event(logging::algo_phase);
      }), 
     //  set up extra vertices
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.totalVertices-s.n;},
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto Triangs = s.Triangs;
        auto vv = s.vv;
        auto n = s.n;
        auto m = s.m;
        auto v = s.v;
        for (auto i = lo; i != hi; i++) {
          v[i] = new (&vv[i+n]) vertex(point2d(0,0), i+n);
          // give each one a pointer to two triangles to use
          v[i]->t = Triangs + m + 2*i;
        }
      }),
      dc::stmt([] (sar& s, par& p) {
        // these will increase as more are added
        s.num_triangs = s.m;
        s.num_points = s.n;
        s.workQ = TriangleTable(s.num_triangs, hashTriangles());
      }),
      dc::stmt([] (sar& s, par& p) {
        logging::push_event(logging::algo_phase);
      }),
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        return sequence::fill3(st, pt, s.workQ.TA, s.workQ.TA + s.workQ.m, &s.empty);
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.num_triangs;},
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto Triangs = s.Triangs;
        auto& workQ = s.workQ;
        for (auto i = lo; i != hi; i++) {
          if (skinnyTriangle(&Triangs[i])) {
            workQ.insert(&Triangs[i]);
            Triangs[i].bad = 1;
          }
        }
      }),
      dc::stmt([] (sar& s, par& p) {
        logging::push_event(logging::algo_phase);
      }),
      // Each iteration processes all bad triangles from the workQ while
      // adding new bad triangles to a new queue
      dc::sequential_loop([] (sar& s, par&) { return true; }, dc::stmts({
        dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
            using hash_type = hashTriangles;
          return heartbeat_call<Table_entries<hash_type, intT>>(st, pt, &s.workQ, &s.badTT);
        }),
        dc::stmt([] (sar& s, par& p) {
          s.workQ.del();
          // packs out triangles that are no longer bad
          s.flags = malloc_array<bool>(s.badTT.n);
        }),
        dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.badTT.n;},
                              [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                              [] (sar& s, par& p, int lo, int hi) {
          auto flags = s.flags;
          auto& badTT = s.badTT;
          for (auto i = lo; i != hi; i++) {
            flags[i] = badTT.A[i]->bad;
          }
        }),
        dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
          return sequence::pack4(st, pt, s.badTT.A, s.flags, (intT)s.badTT.n, &s.badT);
        }),
        dc::stmt([] (sar& s, par& p) {
          free(s.flags);
          s.badTT.del();
          s.numBad = s.badT.n;
        }),
        dc::mk_if([] (sar& s, par& p) { return s.numBad == 0; },
          dc::exit_loop()),
        dc::mk_if([] (sar& s, par& p) { return s.num_points + s.numBad > s.totalVertices; },
          dc::stmt([] (sar& s, par& p) {
              heartbeat::atomic::die("ran out of vertices");
          })),
        // allocate 1 vertex per bad triangle and assign triangle to it
        dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.numBad;},
                              [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                              [] (sar& s, par& p, int lo, int hi) {
          auto& badT = s.badT;
          auto v = s.v;
          auto n = s.n;
          auto num_points = s.num_points;
          for (auto i = lo; i != hi; i++) {
            badT.A[i]->bad = 2; // used to detect whether touched
            v[i + num_points - n]->badT = badT.A[i];
          }
        }),
        dc::stmt([] (sar& s, par& p) {
          // the new work queue
          s.workQ = TriangleTable(s.numBad, hashTriangles());
        }),
        dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
          return sequence::fill3(st, pt, s.workQ.TA, s.workQ.TA + s.workQ.m, &s.empty);
        }),
        dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
          return heartbeat_call<addRefiningVertices>(st, pt, s.v + s.num_points - s.n, s.numBad, s.num_points, s.workQ, &s.tmpi);
        }),
        dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.numBad;},
                              [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                              [] (sar& s, par& p, int lo, int hi) {
          auto& badT = s.badT;
          auto& workQ = s.workQ;
          for (auto i = lo; i != hi; i++) {
            if (badT.A[i]->bad==2) workQ.insert(badT.A[i]);
          }
        }),
        dc::stmt([] (sar& s, par& p) {
          s.badT.del();
          s.num_points += s.numBad;
          s.num_triangs += 2*s.numBad;            
        })
      })),
      dc::stmt([] (sar& s, par& p) {
        logging::push_event(logging::algo_phase);
      }),
      dc::stmt([] (sar& s, par& p) {
        // Extract Vertices for result
        s.flags = malloc_array<bool>(s.num_triangs);
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.num_points;},
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto flags = s.flags;
        auto vv = s.vv;
        for (auto i = lo; i != hi; i++) {
          flags[i] = (vv[i].badT == NULL);
        }
      }),
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        return sequence::packIndex(st, pt, s.flags, s.num_points, &s.I);
      }),
      dc::stmt([] (sar& s, par& p) {
        logging::push_event(logging::algo_phase);
      }),
      dc::stmt([] (sar& s, par& p) {
        s.nO = s.I.n;
        s.II = s.I.A;
        s.rp = malloc_array<point2d>(s.nO);
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.nO;},
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto vv = s.vv;
        auto rp = s.rp;
        auto II = s.II;
        for (auto i = lo; i != hi; i++) {
          vv[II[i]].id = i;
          rp[i] = vv[II[i]].pt;
        }
      }),
      dc::stmt([] (sar& s, par& p) {
        s.I.del();
      }),
      dc::stmt([] (sar& s, par& p) {
        logging::push_event(logging::algo_phase);
      }),
      // Extract Triangles for result
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.num_triangs;},
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto flags = s.flags;
        auto Triangs = s.Triangs;
        for (auto i = lo; i != hi; i++) {
          flags[i] = Triangs[i].initialized;
        }
      }),
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        return sequence::packIndex(st, pt, s.flags, s.num_triangs, &s.I);
      }),
      dc::stmt([] (sar& s, par& p) {
        logging::push_event(logging::algo_phase);
      }),
      dc::stmt([] (sar& s, par& p) {
        s.rt = malloc_array<triangle>(s.I.n);
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.I.n;},
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto Triangs = s.Triangs;
        auto I = s.I;
        auto rt = s.rt;
        for (auto i = lo; i != hi; i++) {
          tri t = Triangs[I.A[i]];
          rt[i] = triangle(t.vtx[0]->id, t.vtx[1]->id, t.vtx[2]->id);
        }
      }),
      dc::stmt([] (sar& s, par& p) {
        logging::push_event(logging::algo_phase);
      }),
      dc::stmt([] (sar& s, par& p) {
        s.I.del();  free(s.flags);  free(s.Triangs);  free(s.v);  free(s.vv);
        *s.dest = triangles<point2d>(s.nO, s.I.n, s.rp, s.rt);
      })
    });
  }
  
  
};

heartbeat_pcfg_allocate(refine, get_cfg)
  
} // end namespace

#endif
