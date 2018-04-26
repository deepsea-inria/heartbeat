// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2010 Guy Blelloch and the PBBS team
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
#include "geometry.hpp"

#if defined(TEST)
#include "test.hpp"
#endif
//#include "prandgen.hpp"
#include "geometrydata.hpp"
#include "readinputbinary.hpp"

#include "hull.h"

namespace sched = heartbeat::sched;
namespace cmdline = deepsea::cmdline;
namespace dsl = heartbeat::edsl;

using namespace std;
using namespace sequence;

using point2d = sptl::point2d;

template <class ET, class F>
pair<intT,intT> split(ET* A, intT n, F lf, F rf) {
  intT ll = 0, lm = 0;
  intT rm = n-1, rr = n-1;
  while (1) {
    while ((lm <= rm) && !(rf(A[lm]) > 0)) {
      if (lf(A[lm]) > 0) A[ll++] = A[lm];
      lm++;
    }
    while ((rm >= lm) && !(lf(A[rm]) > 0)) {
      if (rf(A[rm]) > 0) A[rr--] = A[rm];
      rm--;
    }
    if (lm >= rm) break;
    ET tmp = A[lm++];
    A[ll++] = A[rm--];
    A[rr--] = tmp;
  }
  intT n1 = ll;
  intT n2 = n-rr-1;
  return pair<intT,intT>(n1,n2);
}

struct aboveLine {
  intT l, r;
  point2d* P;
  aboveLine(point2d* _P, intT _l, intT _r) : P(_P), l(_l), r(_r) {}
  bool operator() (intT i) {return sptl::triangle_area(P[l], P[r], P[i]) > 0.0;}
};

intT serialQuickHull(intT* I, point2d* P, intT n, intT l, intT r) {
  if (n < 2) return n;
  intT maxP = I[0];
  double maxArea = sptl::triangle_area(P[l],P[r],P[maxP]);
  for (intT i=1; i < n; i++) {
    intT j = I[i];
    double a = sptl::triangle_area(P[l],P[r],P[j]);
    if (a > maxArea) {
      maxArea = a;
      maxP = j;
    }
  }
  
  pair<intT,intT> nn = split(I, n, aboveLine(P,l,maxP), aboveLine(P,maxP,r));
  intT n1 = nn.first;
  intT n2 = nn.second;
  
  intT m1, m2;
  m1 = serialQuickHull(I,      P, n1, l,   maxP);
  m2 = serialQuickHull(I+n-n2, P, n2, maxP,r);
  for (intT i=0; i < m2; i++) I[i+m1+1] = I[i+n-n2];
  I[m1] = maxP;
  return m1+1+m2;
}

struct triangArea {
  intT l, r;
  point2d* P;
  intT* I;
  triangArea(intT* _I, point2d* _P, intT _l, intT _r) : I(_I), P(_P), l(_l), r(_r) {}
  double operator() (intT i) {return sptl::triangle_area(P[l], P[r], P[I[i]]);}
};

template <class ET, class intT, class PRED> 
intT filterP(ET* In, ET* Out, intT n, PRED p) {
  ET* Orig = Out;
  for (intT i=0; i < n; i++) {
    if (p(In[i])) {
      *Out++ = In[i];
    }
  }
  return Out - Orig;
}

intT quickHullP(intT* I, intT* Itmp, point2d* P, intT n, intT l, intT r) {
  if (n < 2)
    return serialQuickHull(I, P, n, l, r);
  else {
    intT idx = pbbs::maxIndexSerial<double>((intT)0,n,greater<double>(),triangArea(I,P,l,r));
    intT maxP = I[idx];
    intT n1 = filterP(I, Itmp,    n, aboveLine(P, l, maxP));
    intT n2 = filterP(I, Itmp+n1, n, aboveLine(P, maxP, r));

    intT m1, m2;
    m1 = quickHullP(Itmp, I ,P, n1, l, maxP);
    m2 = quickHullP(Itmp+n1, I+n1, P, n2, maxP, r);

    for (intT i=0; i < m1; i++) I[i] = Itmp[i];
    I[m1] = maxP;
    for (intT i=0; i < m2; i++) I[i+m1+1] = Itmp[i+n1];
    return m1+1+m2;
  }
}


class quickHull : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  intT* I; intT* Itmp; point2d* P; intT n; intT l; intT r;
  intT* dest;
  intT idx; intT maxP; intT n1; intT n2; intT m1; intT m2;
  int lg_lt;
    
  quickHull(intT* I, intT* Itmp, point2d* P, intT n, intT l, intT r, intT* dest)
  : I(I), Itmp(Itmp), P(P), n(n), l(l), r(r), dest(dest) { }
  
  heartbeat_dc_declare(heartbeat::edsl, quickHull, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    using controller_type = heartbeat::grain::controller<heartbeat::grain::automatic, quickHull>;
    controller_type::set_ppt(__LINE__, __FILE__);
    return dc::stmts({
      dc::mk_if([&] (sar& s, par&) {
        s.lg_lt = controller_type::predict_lg_nb_iterations();
        auto lt = std::max(2, controller_type::predict_nb_iterations(s.lg_lt));
        return s.n < lt;
	}, dc::stmts({
        dc::stmt([] (sar& s, par&) {
          *s.dest = quickHullP(s.I, s.Itmp, s.P, s.n, s.l, s.r);
          controller_type::register_callback(s.lg_lt, controller_type::predict_nb_iterations(s.lg_lt));
        }),
        dc::exit_function()
      })),
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        auto f = greater<double>();
        auto g = triangArea(s.I,s.P,s.l,s.r);
        return heartbeat_call<maxIndex<double,intT,typeof(f),typeof(g)>>(st, pt, 0, s.n, f, g, &s.idx);
      }),
      dc::stmt([] (sar& s, par&) {
        s.maxP = s.I[s.idx];
      }),
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        auto a = aboveLine(s.P, s.l, s.maxP);
        return heartbeat_call<filterDPS<intT, intT, typeof(a)>>(st, pt, s.I, s.Itmp, s.n, a, &s.n1);
      }),
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        auto a = aboveLine(s.P, s.maxP, s.r);
        return heartbeat_call<filterDPS<intT, intT, typeof(a)>>(st, pt, s.I, s.Itmp+s.n1, s.n, a, &s.n2);
      }),
      dc::spawn2_join(
        [] (sar& s, par&, plt pt, stt st) {
          return heartbeat_call<quickHull>(st, pt, s.Itmp, s.I, s.P, s.n1, s.l, s.maxP, &s.m1);
        },
        [] (sar& s, par&, plt pt, stt st) {
          return heartbeat_call<quickHull>(st, pt, s.Itmp+s.n1, s.I+s.n1, s.P, s.n2, s.maxP, s.r, &s.m2);
        }),
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        return heartbeat_call<sequence::copy<intT*, intT*>>(st, pt, s.Itmp, s.Itmp + s.m1, s.I);
      }),
      dc::stmt([] (sar& s, par& p) {
        s.I[s.m1] = s.maxP;
      }),
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        auto Itmp = s.Itmp + s.n1;
        return heartbeat_call<sequence::copy<intT*, intT*>>(st, pt, Itmp, Itmp + s.m2, s.I + s.m1 + 1);
      }),
      dc::stmt([] (sar& s, par& p) {
        *s.dest = s.m1 + 1 + s.m2;
      })
    });
  }
  
};

heartbeat_pcfg_allocate(quickHull, get_cfg)

struct makePair {
  pair<intT,intT> operator () (intT i) { return pair<intT,intT>(i,i);}
};

struct minMaxIndex {
  point2d* P;
  minMaxIndex (point2d* _P) : P(_P) {}
  pair<intT,intT> operator () (pair<intT,intT> l, pair<intT,intT> r) {
    intT minIndex =
    (P[l.first].x < P[r.first].x) ? l.first :
    (P[l.first].x > P[r.first].x) ? r.first :
    (P[l.first].y < P[r.first].y) ? l.first : r.first;
    intT maxIndex = (P[l.second].x > P[r.second].x) ? l.second : r.second;
    return pair<intT,intT>(minIndex, maxIndex);
  }
};

class hull : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  point2d* P; intT n;
  _seq<intT>* dest;
  intT l; intT r; bool* fTop; bool* fBot; intT* I; intT* Itmp; intT n1; intT n2;
  intT m1; intT m2;
  _seq<intT> tmp;
  pair<intT,intT> minMax;
  
  hull(point2d* P, intT n, _seq<intT>* dest)
  : P(P), n(n), dest(dest) { }
  
  heartbeat_private_activation_record_begin(heartbeat::edsl, hull, 1)
    intT s; intT e;
  heartbeat_private_activation_record_end(heartbeat::edsl, hull, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        auto f = minMaxIndex(s.P);
        auto g = makePair();
        return heartbeat_call<reduce<pair<intT,intT>, intT, typeof(f), typeof(g)>>(st, pt, 0, s.n, f, g, &s.minMax);
      }),
      dc::stmt([] (sar& s, par& p) {
        s.l = s.minMax.first;
        s.r = s.minMax.second;
        s.fTop = malloc_array<bool>(s.n);
        s.fBot = malloc_array<bool>(s.n);
        s.I = malloc_array<intT>(s.n);
        s.Itmp = malloc_array<intT>(s.n);
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.s = 0; p.e = s.n; },
                            [] (par& p) { return std::make_pair(&p.s, &p.e); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto Itmp = s.Itmp;
        auto P = s.P;
        auto fTop = s.fTop;
        auto fBot = s.fBot;
        auto r = s.r;
        auto l = s.l;
        for (auto i = lo; i != hi; i++) {
          Itmp[i] = i;
          double a = triangle_area(P[l], P[r], P[i]);
          fTop[i] = a > 0;
          fBot[i] = a < 0;
        }
      }),
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        return pack5(st, pt, s.Itmp, s.I, s.fTop, s.n, &s.tmp);
      }),
      dc::stmt([] (sar& s, par& p) {
        s.n1 = (intT)s.tmp.n;
      }),
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        return pack5(st, pt, s.Itmp, s.I + s.n1, s.fBot, s.n, &s.tmp);
      }),
      dc::stmt([] (sar& s, par& p) {
        s.n2 = (intT)s.tmp.n;
        free(s.fTop);
        free(s.fBot);
      }),
      dc::spawn2_join(
        [] (sar& s, par&, plt pt, stt st) {
          return heartbeat_call<quickHull>(st, pt, s.I, s.Itmp, s.P, s.n1, s.l, s.r, &s.m1);
        },
        [] (sar& s, par&, plt pt, stt st) {
          return heartbeat_call<quickHull>(st, pt, s.I + s.n1, s.Itmp + s.n1, s.P, s.n2, s.r, s.l, &s.m2);
        }),
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        return heartbeat_call<sequence::copy<intT*, intT*>>(st, pt, s.I, s.I + s.m1, s.Itmp + 1);
      }),
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        auto I = s.I + s.n1;
        return heartbeat_call<sequence::copy<intT*, intT*>>(st, pt, I, I + s.m2, s.Itmp + s.m1 + 2);
      }),
      dc::stmt([] (sar& s, par& p) {
        free(s.I);
        s.Itmp[0] = s.l;
        s.Itmp[s.m1 + 1] = s.r;
        *s.dest = _seq<intT>(s.Itmp, s.m1 + 2 + s.m2);
      })
    });
  }
  
};

heartbeat_pcfg_allocate(hull, get_cfg)

namespace sptl {

#if defined(TEST)

/*---------------------------------------------------------------------*/
/* Quickcheck IO */

template <class Container>
std::ostream& operator<<(std::ostream& out, const sptl::container_wrapper<Container>& c) {
  out << c.c;
  return out;
}

/*---------------------------------------------------------------------*/
/* Quickcheck generators */
  
intT m = 1;

void generate(size_t _nb, parray<point2d>& dst) {
  intT nb = m * (intT)_nb;
  if (quickcheck::generateInRange(0, 1) == 0) {
    dst = plummer2d(nb);
  } else {
    bool inSphere = quickcheck::generateInRange(0, 1) == 0;
    bool onSphere = quickcheck::generateInRange(0, 1) == 0;
    dst = uniform2d(inSphere, onSphere, nb);
  }
}

void generate(size_t nb, container_wrapper<parray<point2d>>& c) {
  generate(nb, c.c);
}

/*---------------------------------------------------------------------*/
/* Quickcheck properties */

struct getX {
  point2d* P;
  getX(point2d* _P) : P(_P) {}
  double operator() (intT i) {return P[i].x;}
};

struct lessX {bool operator() (point2d a, point2d b) {
  return (a.x < b.x) ? 1 : (a.x > b.x) ? 0 : (a.y < b.y);} };

bool equals(point2d a, point2d b) {
  return (a.x == b.x) && (a.y == b.y);
}

//    if (f(v,r)) { r = v; k = j;}
bool check_hull(parray<point2d>& points, parray<intT> indices) {
  point2d* p = points.begin();
  intT n = (intT)points.size();
  intT hull_size = (intT)indices.size();
  point2d* convex_hull = newA(point2d, hull_size);
  for (intT i = 0; i < hull_size; i++) {
    convex_hull[i] = p[indices[i]];
  }
  intT idx = 0;
  for (int i = 0; i < hull_size; i++) {
    if (convex_hull[i].x > convex_hull[idx].x ||
        (convex_hull[i].x == convex_hull[idx].x && convex_hull[i].y > convex_hull[idx].y)) {
      idx = i;
    }
  }
  std::sort(p, p + n, lessX());
  if (!equals(p[0], convex_hull[0])) {
    cout << "checkHull: bad leftmost point" << endl;
    p[0].print();  convex_hull[0].print(); cout << endl;
    return 1;
  }
  if (!equals(p[n - 1], convex_hull[idx])) {
    cout << "checkHull: bad rightmost point" << endl;
    return 1;
  }
  intT k = 1;
  for (intT i = 0; i < idx; i++) {
    if (i > 0 && counter_clockwise(convex_hull[i - 1], convex_hull[i], convex_hull[i + 1])) {
      cout << "checkHull: not convex sides" << endl;
      return 1;
    }
    if (convex_hull[i].x > convex_hull[i + 1].x) {
      cout << "checkHull: upper hull not sorted by x" << endl;
      cout << indices << endl;
      for (int i = 0; i < hull_size; i++) {
        cout << convex_hull[i] << " ";
      }
      cout << std::endl;
      return 1;
    }
    while (k < n && !equals(p[k], convex_hull[i + 1]))
      if (counter_clockwise(convex_hull[i], convex_hull[i + 1], p[k++])) {
        cout << "checkHull: not convex" << endl;
        return 1;
      }
    if (k == n) {
      cout << "checkHull: unexpected points in hull" << endl;
      return 1;
    }
    k++;
  }
  k = n - 2;
  for (intT i = idx; i < hull_size - 1; i++) {
    if (i > idx && counter_clockwise(convex_hull[i - 1], convex_hull[i], convex_hull[i + 1])) {
      cout << "checkHull: not convex sides" << endl;
      return 1;
    }
    if (convex_hull[i].x < convex_hull[i + 1].x) {
      cout << "checkHull: lower hull not sorted by x" << endl;
      return 1;
    }
    while (k >= 0 && !equals(p[k], convex_hull[i + 1])) {
      if (counter_clockwise(convex_hull[i], convex_hull[i + 1], p[k--])) {
        cout << "checkHull: not convex" << endl;
      }
    }
    k--;
  }
  free(convex_hull);
  return 0;
}


using parray_wrapper = container_wrapper<parray<point2d>>;

class consistent_hulls_property : public quickcheck::Property<parray_wrapper> {
public:
  
  bool holdsFor(const parray_wrapper& _in) {
    parray_wrapper in(_in);
    _seq<intT> idxs;
    heartbeat::launch_interpreter<hull>(in.c.begin(), in.c.size(), &idxs);
    parray<intT> idxs2(idxs.A, idxs.A + idxs.n);
    idxs.del();
    return ! check_hull(in.c, idxs2);
  }
  
};

#endif

/*---------------------------------------------------------------------*/
/* Benchmarking */

parray<pbbs::_point2d<double>> to_pbbs(parray<sptl::_point2d<double>>& points) {
  parray<pbbs::_point2d<double>> result(points.size());
  for (int i = 0; i < points.size(); i++) {
    result[i] = pbbs::_point2d<double>(points[i].x, points[i].y);
  }
  return result;
}


void benchmark(std::string infile) {
  parray<sptl::_point2d<double>> x = sptl::read_from_file<parray<sptl::_point2d<double>>>(infile);
  std::string algorithm = cmdline::parse<std::string>("algorithm");
  deepsea::cmdline::dispatcher d;
  d.add("heartbeat", [&] {
    _seq<intT> idxs;
    heartbeat::launch_interpreter<hull>(x.begin(), x.size(), &idxs);
    idxs.del();
    if (deepsea::cmdline::parse_or_default_bool("check", false)) {
      parray<intT> idxs2(idxs.A, idxs.A + idxs.n);
      idxs.del();
#if defined(TEST)
      if (! check_hull(x, idxs2)) {
        assert(false);
      }
#endif
    }
  });
  d.add("pbbs", [&] {
    pbbs::_seq<intT> idxs;
    parray<pbbs::_point2d<double>> y = to_pbbs(x);
    heartbeatbench::run_and_report_elapsed_time([&] {
      idxs = pbbs::hull(y.begin(), (int)y.size());
    });
    idxs.del();
  });
  d.dispatch("algorithm");
}

} // end namespace

int main(int argc, char** argv) {
  heartbeatbench::initialize(argc, argv);
  std::string infile = deepsea::cmdline::parse_or_default_string("infile", "");
  if (infile != "") {
    sptl::benchmark(infile);
    return 0;
  }
  return 0;
}
