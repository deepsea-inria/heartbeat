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

#include "sequence.hpp"
#include "utils.h"
#include "logging.hpp"
#include "nearestneighbors.hpp"
#include "topology.hpp"

#ifndef _HEARTBEAT_DELAUNAY_H_
#define _HEARTBEAT_DELAUNAY_H_

namespace heartbeatbench {

using namespace std;

// if on verifies the Delaunay is correct 
#define CHECK 0

// *************************************************************
//    ROUTINES FOR FINDING AND INSERTING A NEW POINT
// *************************************************************

// Finds a vertex (p) in a mesh starting at any triangle (start)
// Requires that the mesh is properly connected and convex
simplex find(vertex *p, simplex start) {
  simplex t = start;
  while (1) {
    int i;
    for (i=0; i < 3; i++) {
      t = t.rotClockwise();
      if (t.outside(p)) {t = t.across(); break;}
    }
    if (i==3) return t;
    if (!t.valid()) return t;
  }
}

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

// checks if v "won" on all adjacent vertices and inserts point if so
bool insert(vertex *v, simplex t, Qs *q) {
  bool flag = 0;
  for (intT i = 0; i < q->vertexQ.size(); i++) {
    vertex* u = (q->vertexQ)[i];
    if (u->reserve == v->id) u->reserve = -1; // reset to -1
    else flag = 1; // someone else with higher priority reserved u
  }
  if (!flag) {
    tri* t1 = v->t;  // the memory for the two new triangles
    tri* t2 = t1 + 1;  
    // the following 3 lines do all the side effects to the mesh.
    t.split(v, t1, t2);
    for (intT i = 0; i<q->simplexQ.size(); i++) {
      (q->simplexQ)[i].flip();
    }
  }
  q->simplexQ.clear();
  q->vertexQ.clear();
  return flag;
}

// *************************************************************
//    CHECKING THE TRIANGULATION
// *************************************************************
#if 0
void checkDelaunay(tri *triangs, intT n, intT boundarySize) {
  intT *bcount = newA(intT,n);
  /* cilk_for */ for (intT j=0; j<n; j++) bcount[j] = 0;
  /*cilk_for*/ for (intT i=0; i<n; i++) {
    if (triangs[i].initialized >= 0) {
      simplex t = simplex(&triangs[i],0);
      for (int i=0; i < 3; i++) {
        simplex a = t.across();
        if (a.valid()) {
          vertex* v = a.rotClockwise().firstVertex();
          if (!t.outside(v)) {
            cout << "Inside Out: "; v->pt.print(); t.print();}
          if (t.inCirc(v)) {
            cout << "In Circle Violation: "; v->pt.print(); t.print(); }
        } else bcount[i]++;
        t = t.rotClockwise();
      }
    }
  }
  /*
  if (boundarySize != sequence::plusReduce(bcount,n))
    cout << "Wrong boundary size: should be " << boundarySize 
    << " is " << bcount << endl; */
  free(bcount);
}
#endif
// *************************************************************
//    CREATING A BOUNDING CIRCULAR REGION AND FILL WITH INITIAL SIMPLICES
// *************************************************************

struct minpt {
  point2d operator() (point2d a, point2d b) {return a.min_coord(b);}};
struct maxpt {
  point2d operator() (point2d a, point2d b) {return a.max_coord(b);}};

// P is the set of points to bound and n the number
// bCount is the number of points to put on the boundary
// v is an array to put the new points
// t is an array to put the new triangles
class generateBoundary : public heartbeat::edsl::pcfg::shared_activation_record {
public:

  point2d* P; intT n; intT bCount; vertex* v; tri* t;
  simplex* dst; int lo; int hi;

  double size; double radius; point2d center; point2d* boundaryP;
  simplex s; point2d minP; point2d maxP;

  generateBoundary(point2d* P, intT n, intT bCount, vertex* v, tri* t, simplex* dst)
    : P(P), n(n), bCount(bCount), v(v), t(t), dst(dst) { }

  heartbeat_dc_declare(heartbeat::edsl, generateBoundary, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({ 
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        return sequence::reduce4(st, pt, s.P, s.n, minpt(), &s.minP);
      }),
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        return sequence::reduce4(st, pt, s.P, s.n, maxpt(), &s.maxP);
      }), 
      dc::stmt([] (sar& s, par& p) {
        s.size = (s.maxP-s.minP).length();
        double stretch = 10.0;
        s.radius = stretch*s.size;
        s.center = s.maxP + (s.maxP-s.minP)/2.0;
        s.boundaryP = malloc_array<point2d>(s.bCount);
      }),
      // Generate the bounding points on a circle far outside the bounding box
      dc::sequential_loop([] (sar& s, par&) { s.lo = 0; s.hi = s.bCount; },
                          [] (sar& s, par&) { return std::make_pair(&s.lo, &s.hi); },
                          [] (sar& s, par&, int lo, int hi) {
        double pi = 3.14159;
        auto radius = s.radius;
        auto bCount = s.bCount;
        auto center = s.center;
        auto boundaryP = s.boundaryP;
        auto n = s.n;
        auto v = s.v;
        for (auto i = lo; i != hi; i++) {
          double x = radius * cos(2*pi*((float) i)/((float) bCount));
          double y = radius * sin(2*pi*((float) i)/((float) bCount));
          boundaryP[i] = center + vect2d(x,y);
          v[i] = vertex(boundaryP[i], i + n);
        }
      }),
      dc::stmt([] (sar& s, par& p) {        
        s.s = simplex(&(s.v[0]), &(s.v[1]), &(s.v[2]), s.t);
      }),
      // Fill with simplices (bCount - 2  total simplices)
      dc::sequential_loop([] (sar& s, par&) { s.lo = 3; s.hi = s.bCount; },
                          [] (sar& s, par&) { return std::make_pair(&s.lo, &s.hi); },
                          [] (sar& s, par&, int lo, int hi) {
        auto v = s.v;
        auto t = s.t;
        auto _s = s.s;
        for (auto i = lo; i != hi; i++) {
          _s = _s.extend(&v[i], t+i-2);
        }
        s.s = _s;
      }),
      dc::stmt([] (sar& s, par& p) {        
        *s.dst = s.s;
      })
    });
  }
      
};

heartbeat_pcfg_allocate(generateBoundary, get_cfg)
  
// *************************************************************
//    MAIN LOOP
// *************************************************************

class incrementallyAddPoints : public heartbeat::edsl::pcfg::shared_activation_record {
public:

  vertex** v; intT n; vertex* start;

  typedef kNearestNeighbor<vertex,1> kNN;

  intT maxR; Qs *qqs; Qs **qs; simplex *t; bool *flags; vertex** h;
  kNN knn; int multiplier; intT nextNN; intT top; intT rounds; intT failed;
  vertex **vv; intT cnt; _seq<vertex*> tmp; intT k;
  int lo; int hi;

  incrementallyAddPoints(vertex** v, intT n, vertex* start)
    : v(v), n(n), start(start) { }

  heartbeat_private_activation_record_begin(heartbeat::edsl, incrementallyAddPoints, 3)
    int lo; int hi;
  heartbeat_private_activation_record_end(heartbeat::edsl, incrementallyAddPoints, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par& p) {
        // various structures needed for each parallel insertion
        s.maxR = (intT) (s.n/100) + 1; // maximum number to try in parallel
        s.qqs = malloc_array<Qs>(s.maxR);
        s.qs = malloc_array<Qs*>(s.maxR);
      }),
      dc::sequential_loop([] (sar& s, par&) { s.lo = 0; s.hi = s.maxR; },
                          [] (sar& s, par&) { return std::make_pair(&s.lo, &s.hi); },
                          [] (sar& s, par&, int lo, int hi) {
        auto qs = s.qs;
        auto qqs = s.qqs;
        for (auto i = lo; i != hi; i++) {
          qs[i] = new (&qqs[i]) Qs;
        }
      }),
      dc::stmt([] (sar& s, par& p) {
        s.t = malloc_array<simplex>(s.maxR);
        s.flags = malloc_array<bool>(s.maxR);
        s.h = malloc_array<vertex*>(s.maxR);
      }),
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        using cons = typename kNN::qoTreeCons;
        new (&s.knn) kNN;
        return heartbeat_call<cons>(st, pt, &s.start, 1, &s.knn.tree);
      }),
      dc::stmt([] (sar& s, par& p) {
        s.multiplier = 8; // when to regenerate
        s.nextNN = s.multiplier;
        s.top=s.n;
        s.rounds = 0;
        s.failed = 0;
      }),
      dc::sequential_loop([] (sar& s, par&) { return s.top > 0; }, dc::stmts({
        // every once in a while create a new point location
        // structure using all points inserted so far
        dc::mk_if([] (sar& s, par& p) {
            return (s.n-s.top)>=s.nextNN && (s.n-s.top) < s.n/s.multiplier;
          }, dc::stmts({
            dc::stmt([] (sar& s, par& p) {
              s.knn.del();
            }),
            dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
              using cons = typename kNN::qoTreeCons;
              new (&s.knn) kNN;
              return heartbeat_call<cons>(st, pt, s.v+s.top, s.n-s.top, &s.knn.tree);
            }),
            dc::stmt([] (sar& s, par& p) {
              s.nextNN = s.nextNN*s.multiplier;
            })
        })),
        dc::stmt([] (sar& s, par& p) {
          // determine how many vertices to try in parallel
          s.cnt = 1 + (s.n-s.top)/100;  // 100 is pulled out of a hat
          s.cnt = (s.cnt > s.maxR) ? s.maxR : s.cnt;
          s.cnt = (s.cnt > s.top) ? s.top : s.cnt;
          s.vv = s.v+s.top-s.cnt;
        }),
        // for trial vertices find containing triangle, determine cavity 
        // and reserve vertices on boundary of cavity
        dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.cnt;},
                              [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                              [] (sar& s, par& p, int lo, int hi) {
          auto& knn = s.knn;
          auto t = s.t;
          auto vv = s.vv;
          auto qs = s.qs;
          for (auto j = lo; j != hi; j++) {
            vertex *u = knn.nearest(vv[j]);
            t[j] = find(vv[j],simplex(u->t,0));
            reserveForInsert(vv[j],t[j],qs[j]);
          }
        }),
        // Pack failed vertices back onto Q and successful
        // ones up above (needed for point location structure)
        dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.cnt;},
                              [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                              [] (sar& s, par& p, int lo, int hi) {
          auto flags = s.flags;
          auto t = s.t;
          auto vv = s.vv;
          auto qs = s.qs;
          for (auto j = lo; j != hi; j++) {
            flags[j] = insert(vv[j],t[j],qs[j]);
          }
        }),
        dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
          return sequence::pack5(st, pt, s.vv, s.h, s.flags, s.cnt, &s.tmp);
        }),
        dc::stmt([] (sar& s, par& p) {
          s.k = s.tmp.n;
        }),
        dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.cnt;},
                              [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                              [] (sar& s, par& p, int lo, int hi) {
          auto flags = s.flags;
          for (auto j = lo; j != hi; j++) {
            flags[j] = !flags[j];
          }
        }),
        dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
          return sequence::pack5(st, pt, s.vv, s.h+s.k, s.flags, s.cnt, &s.tmp);
        }),
        dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
          return sequence::copy3(st, pt, s.h, s.h + s.cnt, s.vv);
        }),
        dc::stmt([] (sar& s, par& p) {
          s.failed += s.k;
          s.top = s.top-s.cnt+s.k; // adjust top, accounting for failed vertices
          s.rounds++;
        })
      })),
      dc::stmt([] (sar& s, par& p) {
        s.knn.del();
        free(s.qqs); free(s.qs); free(s.t); free(s.flags); free(s.h);
      })
    });
  }

};
  
heartbeat_pcfg_allocate(incrementallyAddPoints, get_cfg)  


// *************************************************************
//    DRIVER
// *************************************************************

// A structure for generating a pseudorandom permutation
struct hashID {

private:
  long k;
  intT n;
  intT GCD(intT a, intT b) {
    while( 1 ) {
      a = a % b;
      if( a == 0 ) return b;
      b = b % a;
      if( b == 0 ) return a;
    }
  }

public:
  hashID(intT nn) : n(nn) {
    k = pbbs::utils::hash(nn)%nn;
    while (GCD(k,nn) > 1) k = (k + 1) % nn;     
  }
  intT get(intT i) { return (i*k)%n; }
};

class delaunay : public heartbeat::edsl::pcfg::shared_activation_record {
public:

  point2d* P; intT n; triangles<point2d>* dst;

  hashID hash; intT numVertices;
  vertex** v; vertex* vv; simplex sBoundary; vertex* v0;
  triangle* rt; intT *M; point2d* rp;
  tri* Triangs; intT numTriangles;
  static constexpr intT boundarySize = 10;
  
  delaunay(point2d* P, intT n, triangles<point2d>* dst)
    : P(P), n(n), dst(dst), hash(n) { }

  heartbeat_private_activation_record_begin(heartbeat::edsl, delaunay, 6)
    int lo; int hi;
  heartbeat_private_activation_record_end(heartbeat::edsl, delaunay, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par& p) {
        // allocate space for vertices
        s.numVertices = s.n + s.boundarySize;
        s.v = malloc_array<vertex*>(s.n); // don't need pointers to boundary
        s.vv = malloc_array<vertex>(s.numVertices);
      }),
      // The points are psuedorandomly permuted 
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.n;},
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto v = s.v;
        auto vv = s.vv;
        auto& hash = s.hash;
        auto P = s.P;
        for (auto i = lo; i != hi; i++) {
          v[i] = new (&vv[i]) vertex(P[hash.get(i)], i);
        }
      }),
      dc::stmt([] (sar& s, par& p) {
        // allocate all the triangles needed
        s.numTriangles = 2 * s.n + (s.boundarySize - 2);
        s.Triangs = malloc_array<tri>(s.numTriangles); 
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.n;},
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto v = s.v;
        auto Triangs = s.Triangs;
        for (auto i = lo; i != hi; i++) {
          v[i]->t = Triangs + 2*i;
        }
      }),
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        return heartbeat_call<generateBoundary>(st, pt, s.P, s.n, s.boundarySize, s.vv + s.n, s.Triangs + 2*s.n, &s.sBoundary);
      }),
      dc::stmt([] (sar& s, par& p) {
        s.v0 = s.sBoundary.t->vtx[0];
      }),
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        return heartbeat_call<incrementallyAddPoints>(st, pt, s.v, s.n, s.v0);
      }),
      dc::stmt([] (sar& s, par& p) {
        free(s.v);
        s.rt = malloc_array<triangle>(s.numTriangles);
        // Since points were permuted need to translate back to
        // original coordinates
        s.M = malloc_array<intT>(s.numVertices);
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.n;},
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto M = s.M;
        auto& hash = s.hash;
        for (auto i = lo; i != hi; i++) {
          M[i] = hash.get(i);
        }
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.numVertices;},
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto M = s.M;
        for (auto i = lo; i != hi; i++) {
          M[i] = i;
        }
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.numVertices;},
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto M = s.M;
        auto Triangs = s.Triangs;
        auto rt = s.rt;
        for (auto i = lo; i != hi; i++) {
          vertex** vtx = Triangs[i].vtx;
          rt[i] = triangle(M[vtx[0]->id], M[vtx[1]->id], M[vtx[2]->id]);
        }
      }),
      dc::stmt([] (sar& s, par& p) {
        free(s.M);
        s.rp = malloc_array<point2d>(s.numVertices);
      }),
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        return sequence::copy3(st, pt, s.P, s.P + s.n, s.rp);
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.numVertices;},
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto rp = s.rp;
        auto vv = s.vv;
        for (auto i = lo; i != hi; i++) {
          rp[i] = vv[i].pt;
        }
      }),
      dc::stmt([] (sar& s, par& p) {
        free(s.Triangs);
        free(s.vv);
        *s.dst = triangles<point2d>(s.numVertices, s.numTriangles, s.rp, s.rt);
      })
    });
  }

};
 
heartbeat_pcfg_allocate(delaunay, get_cfg)
  
} // end namespace

#endif
