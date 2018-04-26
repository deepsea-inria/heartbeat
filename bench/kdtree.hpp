#include <algorithm>
#include <float.h>
#include "samplesort.hpp"
#include "rayTriangleIntersect.h"
#include "topology.h"

namespace heartbeatbench {
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

using namespace std;
  
// Stores coordinate of event along with index to its triangle and type
// Stores type of event (START or END) in lowest bit of index
struct event {
  float v;
  intT p;
  event(float value, intT index, bool type) 
    : v(value), p((index << 1) + type) {}
  event() {}
};
#define START 0
#define IS_START(_event) (!(_event.p & 1))
#define END 1
#define IS_END(_event) ((_event.p & 1))
#define GET_INDEX(_event) (_event.p >> 1)

struct cmpVal { bool operator() (event a, event b) {return a.v < b.v || (a.v == b.v && GET_INDEX(a) < GET_INDEX(b));}};

struct range {
  float min;
  float max;
  range(float _min, float _max) : min(_min), max(_max) {}
  range() {}
};

  /*
typedef range* Boxes[3];
typedef event* Events[3];
typedef range BoundingBox[3];
  */
using Boxes = range**;
using Events = event**;
using BoundingBox = range*;

static std::ostream& operator<<(std::ostream& os, const BoundingBox B) {
 return os << B[0].min << ":" << B[0].max << " + " 
	   << B[1].min << ":" << B[1].max << " + " 
	   << B[2].min << ":" << B[2].max;
}

struct cutInfo {
  float cost;
  float cutOff;
  intT numLeft;
  intT numRight;
cutInfo(float _cost, float _cutOff, intT nl, intT nr) 
: cost(_cost), cutOff(_cutOff), numLeft(nl), numRight(nr) {}
  cutInfo() {}
};

struct treeNode {
  treeNode *left;
  treeNode *right;
  /*BoundingBox*/ range box[3];
  int cutDim;
  float cutOff;
  intT* triangleIndices;
  intT n;
  intT leaves;
  
  bool isLeaf() {return (triangleIndices != NULL);}

  treeNode(treeNode* L, treeNode* R, 
	   int _cutDim, float _cutOff, BoundingBox B) 
    : left(L), right(R), triangleIndices(NULL), cutDim(_cutDim), 
      cutOff(_cutOff) {
    for (int i=0; i < 3; i++) box[i] = B[i];
    n = L->n + R->n;
    leaves = L->leaves + R->leaves;
  }

  treeNode(Events E, intT _n, BoundingBox B)
    : left(NULL), right(NULL) {

    event* events = E[0];

    // extract indices from events
    triangleIndices = malloc_array<intT>(_n/2);
    intT k = 0;
    for (intT i = 0; i < _n; i++) 
      if (IS_START(events[i]))
        triangleIndices[k++] = GET_INDEX(events[i]);

    n = _n/2;
    leaves = 1;
    for (int i=0; i < 3; i++) {
      box[i] = B[i];
      free(E[i]);
    }
  }
  
};

class treeNode_del : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  treeNode* T;
  
  treeNode_del(treeNode* T)
    : T(T) { }

  heartbeat_dc_declare(heartbeat::edsl, treeNode_del, sar, par, dc, get_dc)

  static
  dc get_dc() {
    return dc::stmts({
      dc::mk_if([] (sar& s, par&) { return s.T->isLeaf(); },
        dc::stmt([] (sar& s, par&) { free(s.T->triangleIndices); }),
        dc::spawn2_join(
            [] (sar& s, par&, plt pt, stt st) {
              return heartbeat_call<treeNode_del>(st, pt, s.T->left); },
            [] (sar& s, par&, plt pt, stt st) {
              return heartbeat_call<treeNode_del>(st, pt, s.T->right);
      })),
      dc::stmt([] (sar& s, par&) { free(s.T); })
    });
  }
  
};

heartbeat_pcfg_allocate(treeNode_del, get_cfg)

int CHECK = 0;  // if set checks 10 rays against brute force method
int STATS = 0;  // if set prints out some tree statistics

// Constants for deciding when to stop recursion in building the KDTree
float CT = 6.0;
float CL = 1.25;
float maxExpand = 1.6;
int maxRecursionDepth = 25;

// Constant for switching to sequential versions
int minParallelSize = 500000;

typedef double floatT;
typedef pbbs::_point3d<floatT> pointT;
typedef pbbs::_vect3d<floatT> vectT;
typedef pbbs::triangles<pointT> trianglesT;
typedef pbbs::ray<pointT> rayT;

float boxSurfaceArea(BoundingBox B) {
  float r0 = B[0].max-B[0].min;
  float r1 = B[1].max-B[1].min;
  float r2 = B[2].max-B[2].min;
  return 2*(r0*r1 + r1*r2 + r0*r2);
}

float epsilon = .0000001;
range fixRange(float minv, float maxv) {
  if (minv == maxv) return range(minv,minv+epsilon);
  else return range(minv,maxv);
}

inline float inBox(pointT p, BoundingBox B) {
  return (p.x >= (B[0].min - epsilon) && p.x <= (B[0].max + epsilon) &&
	  p.y >= (B[1].min - epsilon) && p.y <= (B[1].max + epsilon) &&
	  p.z >= (B[2].min - epsilon) && p.z <= (B[2].max + epsilon));
}

// sequential version of best cut
cutInfo bestCutSerial0(event* E, range r, range r1, range r2, intT n) {
  if (r.max - r.min == 0.0) return cutInfo(FLT_MAX, r.min, n, n);
#ifdef TIME_MEASURE
//    cutTimer.start();
#endif
  float area = 2 * (r1.max-r1.min) * (r2.max-r2.min);
  float diameter = 2 * ((r1.max-r1.min) + (r2.max-r2.min));

  // calculate cost of each possible split
  intT inLeft = 0;
  intT inRight = n/2;
  float minCost = FLT_MAX;
  intT k = 0;
  intT rn = inLeft;
  intT ln = inRight;
  for (intT i=0; i <n; i++) {
    float cost;
    if (IS_END(E[i])) inRight--;
    float leftLength = E[i].v - r.min;
    float leftArea = area + diameter * leftLength;
    float rightLength = r.max - E[i].v;
    float rightArea = area + diameter * rightLength;
    cost = (leftArea * inLeft + rightArea * inRight);
    if (cost < minCost) {
      rn = inRight;
      ln = inLeft;
      minCost = cost;
      k = i;
    }
    if (IS_START(E[i])) inLeft++;
  }//std::cerr << "Best k: " << k << std::endl;
#ifdef TIME_MEASURE
//    cutTimer.stop();
#endif
  return cutInfo(minCost, E[k].v, ln, rn);
}

class bestCutSerial : public heartbeat::edsl::pcfg::shared_activation_record {
public:

  event* E; range r; range r1; range r2; intT n; cutInfo* dest;
  float area; float diameter; intT inLeft; intT inRight; float minCost;
  intT k; intT rn; intT ln; int lo; int hi;
  
  bestCutSerial(event* E, range r, range r1, range r2, intT n, cutInfo* dest)
    : E(E), r(r), r1(r1), r2(r2), n(n), dest(dest) { }
  
  heartbeat_dc_declare(heartbeat::edsl, bestCutSerial, sar, par, dc, get_dc)

  static
  dc get_dc() {
    return dc::stmts({
      dc::mk_if([] (sar& s, par&) { return s.r.max - s.r.min == 0.0; }, dc::stmts({
        dc::stmt([] (sar& s, par&) {
          *s.dest = cutInfo(FLT_MAX, s.r.min, s.n, s.n);
        }),
        dc::exit_function()
      })),
      dc::stmt([] (sar& s, par&) {
        s.area = 2 * (s.r1.max-s.r1.min) * (s.r2.max-s.r2.min);
        s.diameter = 2 * ((s.r1.max-s.r1.min) + (s.r2.max-s.r2.min));

        // calculate cost of each possible split
        s.inLeft = 0;
        s.inRight = s.n/2;
        s.minCost = FLT_MAX;
        s.k = 0;
        s.rn = s.inLeft;
        s.ln = s.inRight;
      }),
      dc::sequential_loop([] (sar& s, par&) { s.lo = 0; s.hi = s.n; },
                          [] (sar& s, par&) { return std::make_pair(&s.lo, &s.hi); },
                          [] (sar& s, par&, int lo, int hi) {
        auto E = s.E;
        auto inRight = s.inRight;
        auto r = s.r;
        auto area = s.area;
        auto diameter = s.diameter;
        auto inLeft = s.inLeft;
        auto k = s.k;
        auto minCost = s.minCost;
        auto rn = s.rn;
        auto ln = s.ln;
        for (int i = lo; i < hi; i++) {
          float cost;
          if (IS_END(E[i])) inRight--;
          float leftLength = E[i].v - r.min;
          float leftArea = area + diameter * leftLength;
          float rightLength = r.max - E[i].v;
          float rightArea = area + diameter * rightLength;
          cost = (leftArea * inLeft + rightArea * inRight);
          if (cost < minCost) {
            rn = inRight;
            ln = inLeft;
            minCost = cost;
            k = i;
          }
          if (IS_START(E[i])) inLeft++;
        }
        s.k = k;
        s.inLeft = inLeft;
        s.inRight = inRight;
        s.minCost = minCost;
        s.rn = rn;
        s.ln = ln;
      }),
      dc::stmt([] (sar& s, par&) {
        *s.dest = cutInfo(s.minCost, s.E[s.k].v, s.ln, s.rn);
      })
    });
  }

};

heartbeat_pcfg_allocate(bestCutSerial, get_cfg)

class bestCut : public heartbeat::edsl::pcfg::shared_activation_record {
public:

  event* E; range r; range r1; range r2; intT n; cutInfo* dest;
  float orthogArea; float diameter; intT* upperC; intT u;
  float* cost; intT k; int lg_lt;
  
  bestCut(event* E, range r, range r1, range r2, intT n, cutInfo* dest)
    : E(E), r(r), r1(r1), r2(r2), n(n), dest(dest) { }
    
  heartbeat_private_activation_record_begin(heartbeat::edsl, bestCut, 2)
    int lo; int hi;
  heartbeat_private_activation_record_end(heartbeat::edsl, bestCut, sar, par, dc, get_dc)

  static
  dc get_dc() {
    return dc::stmts({
      dc::mk_if([] (sar& s, par& p) { return s.n < 1000; }, dc::stmts({
        dc::stmt([] (sar& s, par& p) {
          *s.dest = bestCutSerial0(s.E, s.r, s.r1, s.r2, s.n);
        }),
        dc::exit_function()
      })),
      dc::mk_if([] (sar& s, par& p) { return s.n < minParallelSize; }, dc::stmts({
        dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
          return heartbeat_call<bestCutSerial>(st, pt, s.E, s.r, s.r1, s.r2, s.n, s.dest);
        }),
        dc::exit_function()
      })),
      dc::mk_if([] (sar& s, par& p) { return s.r.max - s.r.min == 0.0; }, dc::stmts({
        dc::stmt([] (sar& s, par& p) {
          *s.dest = cutInfo(FLT_MAX, s.r.min, s.n, s.n);
        }),
        dc::exit_function()
      })),
      dc::stmt([] (sar& s, par& p) {
        // area of two orthogonal faces
        s.orthogArea = 2 * ((s.r1.max-s.r1.min) * (s.r2.max-s.r2.min));

        // length of diameter of orthogonal face
        s.diameter = 2 * ((s.r1.max-s.r1.min) + (s.r2.max-s.r2.min));

        // count number that end before i
        s.upperC = malloc_array<intT>(s.n);
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.n; },
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
          auto upperC = s.upperC;
          auto E = s.E;
          for (auto i = lo; i != hi; i++) {
            upperC[i] = IS_END(E[i]);
          }
      }),   
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        return sequence::plusScan(st, pt, s.upperC, s.upperC, s.n, &s.u);
      }),
      dc::stmt([] (sar& s, par& p) {
        s.cost = malloc_array<float>(s.n);
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.n; },
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
          auto upperC = s.upperC;
          auto E = s.E;
          auto n = s.n;
          auto r = s.r;
          auto diameter = s.diameter;
          auto orthogArea = s.orthogArea;
          auto cost = s.cost;
          for (auto i = lo; i != hi; i++) {
            intT inLeft = i - upperC[i];
            intT inRight = n/2 - (upperC[i] + IS_END(E[i]));
            float leftLength = E[i].v - r.min;
            float leftArea = orthogArea + diameter * leftLength;
            float rightLength = r.max - E[i].v;
            float rightArea = orthogArea + diameter * rightLength;
            cost[i] = (leftArea * inLeft + rightArea * inRight);
          }
      }),
      // find minimum across all (maxIndex with less is minimum)
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        return sequence::maxIndex4(st, pt, s.cost, s.n,less<float>(), &s.k);
      }),
      dc::stmt([] (sar& s, par& p) {
        float c = s.cost[s.k];
        intT ln = s.k - s.upperC[s.k];
        intT rn = s.n/2 - (s.upperC[s.k] + IS_END(s.E[s.k]));
        free(s.upperC); free(s.cost);
        *s.dest = cutInfo(c, s.E[s.k].v, ln, rn);
      })          
    });
  }
  
};

heartbeat_pcfg_allocate(bestCut, get_cfg)
  
typedef pair<_seq<event>, _seq<event> > eventsPair;

eventsPair splitEventsSerial0(range* boxes, event* events, 
                             float cutOff, intT n) {
  intT l = 0;
  intT r = 0;
  event* eventsLeft = newA(event,n);
  event* eventsRight = newA(event,n);
  for (intT i=0; i < n; i++) {
    intT b = GET_INDEX(events[i]);
    if (boxes[b].min < cutOff) {
      eventsLeft[l++] = events[i];
      if (boxes[b].max > cutOff) 
        eventsRight[r++] = events[i]; 
    } else eventsRight[r++] = events[i]; 
  }
  return eventsPair(_seq<event>(eventsLeft,l), 
		    _seq<event>(eventsRight,r));
}

class splitEventsSerial : public heartbeat::edsl::pcfg::shared_activation_record {
public:

  range* boxes; event* events; float cutOff; intT n; eventsPair* dest;
  intT l; intT r; event* eventsLeft; event* eventsRight; int lo; int hi;
  
  splitEventsSerial(range* boxes, event* events, 
                    float cutOff, intT n, eventsPair* dest)
    : boxes(boxes), events(events), cutOff(cutOff), n(n), dest(dest) { }
  
  heartbeat_dc_declare(heartbeat::edsl, splitEventsSerial, sar, par, dc, get_dc)

  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par& p) {
        s.l = 0;
        s.r = 0;
        s.eventsLeft = malloc_array<event>(s.n);
        s.eventsRight = malloc_array<event>(s.n);
      }),
      dc::sequential_loop([] (sar& s, par&) { s.lo = 0; s.hi = s.n; },
                          [] (sar& s, par&) { return std::make_pair(&s.lo, &s.hi); },
                          [] (sar& s, par&, int lo, int hi) {
        auto events = s.events;
        auto boxes = s.boxes;
        auto r = s.r;
        auto l = s.l;
        auto cutOff = s.cutOff;
        auto eventsLeft = s.eventsLeft;
        auto eventsRight = s.eventsRight;
        for (int i = lo; i < hi; i++) {
          intT b = GET_INDEX(events[i]);
          if (boxes[b].min < cutOff) {
            eventsLeft[l++] = events[i];
            if (boxes[b].max > cutOff) 
              eventsRight[r++] = events[i]; 
          } else eventsRight[r++] = events[i]; 
        }
        s.r = r;
        s.l = l;
      }),
      dc::stmt([] (sar& s, par& p) {
        *s.dest = eventsPair(_seq<event>(s.eventsLeft,s.l), 
                             _seq<event>(s.eventsRight,s.r));          
      })
    });
  }

};

heartbeat_pcfg_allocate(splitEventsSerial, get_cfg)

class splitEvents : public heartbeat::edsl::pcfg::shared_activation_record {
public:

  range* boxes; event* events; float cutOff; intT n; eventsPair* dest;
  bool* lower; bool* upper; _seq<event> L; _seq<event> R;
    
  splitEvents(range* boxes, event* events, float cutOff, intT n, eventsPair* dest)
    : boxes(boxes), events(events), cutOff(cutOff), n(n), dest(dest) { }
    
  heartbeat_private_activation_record_begin(heartbeat::edsl, splitEvents, 1)
    int lo; int hi;
  heartbeat_private_activation_record_end(heartbeat::edsl, splitEvents, sar, par, dc, get_dc)

  static
  dc get_dc() {
    return dc::stmts({
      dc::mk_if([] (sar& s, par&) { return s.n < 1000; }, dc::stmts({
        dc::stmt([] (sar& s, par& p) {
          *s.dest = splitEventsSerial0(s.boxes, s.events, s.cutOff, s.n);
        }),
        dc::exit_function()
      })),
      dc::mk_if([] (sar& s, par&) { return s.n < minParallelSize; }, dc::stmts({
        dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
          return heartbeat_call<splitEventsSerial>(st, pt, s.boxes, s.events, s.cutOff, s.n, s.dest);
        }),
        dc::exit_function()
      })),
      dc::stmt([] (sar& s, par& p) {
        s.lower = malloc_array<bool>(s.n);
        s.upper = malloc_array<bool>(s.n);
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.n; },
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto events = s.events;
        auto boxes = s.boxes;
        auto lower = s.lower;
        auto upper = s.upper;
        auto cutOff = s.cutOff;
        for (auto i = lo; i != hi; i++) {
          intT b = GET_INDEX(events[i]);
          lower[i] = boxes[b].min < cutOff;
          upper[i] = boxes[b].max > cutOff || boxes[b].min >= cutOff;
        }
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        return sequence::pack4(st, pt, s.events, s.lower, s.n, &s.L);
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        return sequence::pack4(st, pt, s.events, s.upper, s.n, &s.R);
      }),
      dc::stmt([] (sar& s, par& p) {
        free(s.lower); free(s.upper);
        
        *s.dest = eventsPair(s.L,s.R);
      })        
    });
  }
  
};

heartbeat_pcfg_allocate(splitEvents, get_cfg)

// n is the number of events (i.e. twice the number of triangles)
class generateNode : public heartbeat::edsl::pcfg::shared_activation_record {
public:

  Boxes boxes; Events events; BoundingBox B; intT n; intT maxDepth; treeNode** dest;
  cutInfo cuts[3]; int cutDim;
  range* cutDimRanges; float cutOff; float area; float bestCost; float origCost;
  /*BoundingBox*/ range BBL[3]; /*BoundingBox*/ range BBR[3];
  event* rightEvents[3]; intT nr; event* leftEvents[3]; intT nl; eventsPair X[3];
  treeNode *L; treeNode *R;
  
  generateNode(Boxes boxes, Events events, BoundingBox B, 
               intT n, intT maxDepth, treeNode** dest)
    : boxes(boxes), events(events), B(B), n(n), maxDepth(maxDepth), dest(dest) { }
    
  heartbeat_private_activation_record_begin(heartbeat::edsl, generateNode, 2)
    int lo; int hi;
  heartbeat_private_activation_record_end(heartbeat::edsl, generateNode, sar, par, dc, get_dc)

  static
  dc get_dc() {
    return dc::stmts({
      dc::mk_if([] (sar& s, par&) { return s.n <= 2 || s.maxDepth == 0; }, dc::stmts({
        dc::stmt([] (sar& s, par& p) {
          *s.dest = new treeNode(s.events, s.n, s.B);
        }),
        dc::exit_function()
      })),
      dc::stmt([] (sar& s, par& p) {
        p.lo = 0;
        p.hi = 3;
      }),
      dc::parallel_for_loop([] (par& p) { return std::make_pair(&p.lo, &p.hi); }, dc::stmts({
        dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
          auto d = p.lo;
          return heartbeat_call<bestCut>(st, pt, s.events[d], s.B[d], s.B[(d+1)%3], s.B[(d+2)%3], s.n, &s.cuts[d]);
        }),
        dc::stmt([] (sar& s, par& p) {
          p.lo++;
        })
      })),
      dc::stmt([] (sar& s, par& p) {
        s.cutDim = 0;
        for (int d = 1; d < 3; d++) 
          if (s.cuts[d].cost < s.cuts[s.cutDim].cost) s.cutDim = d;
        s.cutDimRanges = s.boxes[s.cutDim];
        s.cutOff = s.cuts[s.cutDim].cutOff;
        s.area = boxSurfaceArea(s.B);
        s.bestCost = CT + CL * s.cuts[s.cutDim].cost/s.area;
        s.origCost = (float) (s.n/2);
      }),
      dc::mk_if([] (sar& s, par&) { return s.bestCost >= s.origCost || 
            s.cuts[s.cutDim].numLeft + s.cuts[s.cutDim].numRight > maxExpand * s.n/2; }, dc::stmts({
        dc::stmt([] (sar& s, par& p) {
          *s.dest = new treeNode(s.events, s.n, s.B);
        }),
        dc::exit_function()
      })),
      dc::stmt([] (sar& s, par& p) {
        // declare structures for recursive calls
        for (int i=0; i < 3; i++) s.BBL[i] = s.B[i];
        s.BBL[s.cutDim] = range(s.BBL[s.cutDim].min, s.cutOff);

        for (int i=0; i < 3; i++) s.BBR[i] = s.B[i];
        s.BBR[s.cutDim] = range(s.cutOff, s.BBR[s.cutDim].max);
        p.lo = 0;
        p.hi = 3;
      }),
      dc::parallel_for_loop([] (par& p) { return std::make_pair(&p.lo, &p.hi); }, dc::stmts({
        dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
          auto d = p.lo;
          return heartbeat_call<splitEvents>(st, pt, s.cutDimRanges, s.events[d], s.cutOff, s.n, &s.X[d]);
        }),
        dc::stmt([] (sar& s, par& p) {
          p.lo++;
        })
      })),
      dc::stmt([] (sar& s, par& p) {
        for (int d = 0; d < 3; d++) {
          s.leftEvents[d] = s.X[d].first.A;
          s.rightEvents[d] = s.X[d].second.A;
          if (d == 0) {
            s.nl = s.X[d].first.n;
            s.nr = s.X[d].second.n;
          } else if (s.X[d].first.n != s.nl || s.X[d].second.n != s.nr) {
            std::cout << "kdTree: mismatched lengths, something wrong" << std::endl;
            exit(0);
          }
        }
        // free old events and make recursive calls
        for (int i=0; i < 3; i++) free(s.events[i]);
      }),
      dc::spawn2_join(
        [] (sar& s, par&, plt pt, stt st) {
          return heartbeat_call<generateNode>(st, pt, s.boxes, s.leftEvents, s.BBL, s.nl, s.maxDepth-1, &s.L); },
        [] (sar& s, par&, plt pt, stt st) {
          return heartbeat_call<generateNode>(st, pt, s.boxes, s.rightEvents, s.BBR, s.nr, s.maxDepth-1, &s.R);
      }),
      dc::stmt([] (sar& s, par& p) {
        *s.dest = new treeNode(s.L, s.R, s.cutDim, s.cutOff, s.B);
      })
    });
  }
  
};

heartbeat_pcfg_allocate(generateNode, get_cfg)
  
intT tcount = 0;
intT ccount = 0;

// Given an a ray, a bounding box, and a sequence of triangles, returns the 
// index of the first triangle the ray intersects inside the box.
// The triangles are given by n indices I into the triangle array Tri.
// -1 is returned if there is no intersection
intT findRay(rayT r, intT* I, intT n, pbbs::triangles<pointT> Tri, BoundingBox B) {
  if (STATS) { tcount += n; ccount += 1;}
  pointT* P = Tri.P;
  floatT tMin = FLT_MAX;
  intT k = -1;
  for (intT i = 0; i < n; i++) {
    intT j = I[i];
    pbbs::triangle* tr = Tri.T + j;
    pointT m[3] = {P[tr->C[0]],  P[tr->C[1]],  P[tr->C[2]]};
    floatT t = rayTriangleIntersect(r, m);
    if (t > 0.0 && t < tMin && inBox(r.o + r.d*t, B)) {
      tMin = t;
      k = j;
    }
  }
  return k;
}

using vect2d = pbbs::vect2d;
using point2d = pbbs::point2d;

// Given a ray and a tree node find the index of the first triangle the 
// ray intersects inside the box represented by that node.
// -1 is returned if there is no intersection
intT findRay(rayT r, treeNode* TN, trianglesT Tri) {
  //cout << "TN->n=" << TN->n << endl;
  if (TN->isLeaf()) 
    return findRay(r, TN->triangleIndices, TN->n, Tri, TN->box);
  pointT o = r.o;
  vectT d = r.d;

  floatT oo[3] = {o.x,o.y,o.z};
  floatT dd[3] = {d.x,d.y,d.z};

  // intersect ray with splitting plane
  int k0 = TN->cutDim;
  int k1 = (k0 == 2) ? 0 : k0+1;
  int k2 = (k0 == 0) ? 2 : k0-1;
  point2d o_p(oo[k1], oo[k2]);
  vect2d d_p(dd[k1], dd[k2]);
  // does not yet deal with dd[k0] == 0
  floatT scale = (TN->cutOff - oo[k0])/dd[k0];
  point2d p_i = o_p + d_p * scale;

  range rx = TN->box[k1];
  range ry = TN->box[k2];
  floatT d_0 = dd[k0];

  // decide which of the two child boxes the ray intersects
  enum {LEFT, RIGHT, BOTH};
  int recurseTo = LEFT;
  if      (p_i.x < rx.min) { if (d_p.x*d_0 > 0) recurseTo = RIGHT;}
  else if (p_i.x > rx.max) { if (d_p.x*d_0 < 0) recurseTo = RIGHT;}
  else if (p_i.y < ry.min) { if (d_p.y*d_0 > 0) recurseTo = RIGHT;}
  else if (p_i.y > ry.max) { if (d_p.y*d_0 < 0) recurseTo = RIGHT;}
  else recurseTo = BOTH;

  if (recurseTo == RIGHT) return findRay(r, TN->right, Tri);
  else if (recurseTo == LEFT) return findRay(r, TN->left, Tri);
  else if (d_0 > 0) {
    intT t = findRay(r, TN->left, Tri);
    if (t >= 0) return t;
    else return findRay(r, TN->right, Tri);
  } else {
    intT t = findRay(r, TN->right, Tri);
    if (t >= 0) return t;
    else return findRay(r, TN->left, Tri);
  }
}

template <class T>
using triangles = pbbs::triangles<T>;

template <class T>
using ray = pbbs::ray<T>;

class rayCast : public heartbeat::edsl::pcfg::shared_activation_record {
public:

  triangles<pointT> Tri; ray<pointT>* rays; intT numRays; intT** dest;
  /*Boxes*/ range* boxes[3]; intT n; /*Events*/ event* events[3];
  /*BoundingBox*/ range boundingBox[3];
  pointT* P; intT recursionDepth; treeNode* R; intT* results; int d; 
  
  rayCast(triangles<pointT> Tri, ray<pointT>* rays, intT numRays, intT** dest)
    : Tri(Tri), rays(rays), numRays(numRays), dest(dest) { }
    
  heartbeat_private_activation_record_begin(heartbeat::edsl, rayCast, 3)
    int lo; int hi;
  heartbeat_private_activation_record_end(heartbeat::edsl, rayCast, sar, par, dc, get_dc)

  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par& p) {
        s.n = s.Tri.numTriangles;
        for (int d = 0; d < 3; d++) s.boxes[d] = malloc_array<range>(s.n);
        s.P = s.Tri.P;
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.n; },
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto events = s.events;
        auto P = s.P;
        auto Tri = s.Tri;
        auto boxes = s.boxes;
        for (auto i = lo; i != hi; i++) {
          pointT p0 = P[Tri.T[i].C[0]];
          pointT p1 = P[Tri.T[i].C[1]];
          pointT p2 = P[Tri.T[i].C[2]];
          boxes[0][i] = fixRange(std::min(p0.x,std::min(p1.x,p2.x)),std::max(p0.x,std::max(p1.x,p2.x)));
          boxes[1][i] = fixRange(std::min(p0.y,std::min(p1.y,p2.y)),std::max(p0.y,std::max(p1.y,p2.y)));
          boxes[2][i] = fixRange(std::min(p0.z,std::min(p1.z,p2.z)),std::max(p0.z,std::max(p1.z,p2.z)));
        }
      }),
      dc::stmt([] (sar& s, par& p) {
        s.d = 0;
      }),
      dc::sequential_loop([] (sar& s, par&) { return s.d < 3; }, dc::stmts({
        dc::stmt([] (sar& s, par& p) {
          s.events[s.d] = malloc_array<event>(2*s.n); // freed while generating tree
        }),
        dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.n; },
                              [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                              [] (sar& s, par& p, int lo, int hi) {
          auto boxes = s.boxes;
          auto events = s.events;
          auto d = s.d;
          for (auto i = lo; i != hi; i++) {
            events[d][2*i] = event(boxes[d][i].min, i, START);
            events[d][2*i+1] = event(boxes[d][i].max, i, END);
          }
        }),
        dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
          return sampleSort3(st, pt, s.events[s.d], s.n*2, cmpVal());
        }),
        dc::stmt([] (sar& s, par& p) {
          s.boundingBox[s.d] = range(s.events[s.d][0].v, s.events[s.d][2*s.n-1].v);
          s.d++;
        })
      })),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        intT recursionDepth = std::min(maxRecursionDepth, pbbs::utils::log2Up(s.n)-1);
        return heartbeat_call<generateNode>(st, pt, s.boxes, s.events, s.boundingBox, s.n*2,
                                         recursionDepth, &s.R);
      }),
      dc::stmt([] (sar& s, par& p) {
        for (int d = 0; d < 3; d++) free(s.boxes[d]);
        s.results = malloc_array<intT>(s.numRays);
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.numRays; },
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto results = s.results;
        auto rays = s.rays;
        auto R = s.R;
        auto Tri = s.Tri;
        for (auto i = lo; i != hi; i++) {
          results[i] = findRay(rays[i], R, Tri);
        }
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        return heartbeat_call<treeNode_del>(st, pt, s.R);          
      }),
      dc::stmt([] (sar& s, par& p) {
        *s.dest = s.results;
      })
    });
  }
  
};

heartbeat_pcfg_allocate(rayCast, get_cfg)

} //end namespace
