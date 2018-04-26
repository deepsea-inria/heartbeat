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

#ifndef _HEARTBEAT__PBBS_QUICKSORT_H_
#define _HEARTBEAT__PBBS_QUICKSORT_H_

namespace heartbeatbench {

template <class E, class BinPred, class intT>
void insertionSort2(E* A, intT n, BinPred f) {
  for (intT i=0; i < n; i++) {
    E v = A[i];
    E* B = A + i;
    while (--B >= A && f(v,*B)) *(B+1) = *B;
    *(B+1) = v;
  }
}

#define ISORT 25

template <class E, class BinPred>
E median2(E a, E b, E c, BinPred f) {
  return  f(a,b) ? (f(b,c) ? b : (f(a,c) ? c : a))
           : (f(a,c) ? a : (f(b,c) ? c : b));
}

// Quicksort based on median of three elements as pivot
//  and uses insertionSort for small inputs
template<class E, class BinPred, class intT>
class quicksort : public heartbeat::edsl::pcfg::shared_activation_record {
public:

  E* A; intT n; BinPred f;
  E* M; E* L; E* R; E p; bool b;
  using trampoline = enum { loop1, loop2 };
  trampoline t;

  quicksort(E* A, intT n, BinPred f)
    : A(A), n(n), f(f) { }

  heartbeat_dc_declare(heartbeat::edsl, quicksort, sar, par, dc, get_dc)

  static
  dc get_dc() {
    return dc::mk_if([] (sar& s, par&) { return s.n < ISORT; }, dc::stmt([] (sar& s, par& p) {
      insertionSort2(s.A, s.n, s.f);
    }), dc::stmts({
    dc::stmt([] (sar& s, par&) {
      auto n = s.n;
      auto A = s.A;
      s.p = median2(A[n/4],A[n/2],A[(3*n)/4],s.f);
      s.L = A;   // below L are less than pivot
      s.M = A;   // between L and M are equal to pivot
      s.R = A+n-1; // above R are greater than pivot
      s.t = loop1;
      s.b = true;
    }),
    dc::sequential_loop([] (sar& s, par&) { return s.b; }, dc::stmt([] (sar& s, par&) {
      using controller_type = heartbeat::grain::controller<heartbeat::grain::automatic, quicksort>;
      auto lg_lt = controller_type::predict_lg_nb_iterations();
      auto lt = controller_type::predict_nb_iterations(lg_lt);
      int fuel0 = lt;
      int fuel = fuel0;
      E p = s.p;
      E* L = s.L;
      E* M = s.M;
      E* R = s.R;
      auto f = s.f;
      trampoline t = s.t;
      while (true) {
        switch (t) {
          case loop1: {
            while (!f(p,*M)) {
              if (f(*M,p)) std::swap(*M,*(L++));
              if (M >= R) break; 
              M++;
              if (--fuel == 0) {
                goto exit;
              }
            }
            t = loop2;
          }
          case loop2: {
            while (f(p,*R)) {
              R--;
              if (--fuel == 0) {
                goto exit;
              }
            }
            t = loop1;
            if (M >= R) {
              s.b = false;
              goto exit;
            }
            std::swap(*M,*R--); 
            if (f(*M,p)) std::swap(*M,*(L++));
            M++;
            --fuel;
          }
        }
      }
    exit:
      s.t = t;
      s.M = M;
      s.L = L;
      s.R = R;
      controller_type::register_callback(lg_lt, fuel0 - fuel);
      return;
    })),
    dc::spawn2_join(
       [] (sar& s, par&, plt pt, stt st) {
         return heartbeat_call<quicksort<E,BinPred,intT>>(st, pt, s.A, s.L-s.A, s.f);
    }, [] (sar& s, par&, plt pt, stt st) {
         return heartbeat_call<quicksort<E,BinPred,intT>>(st, pt, s.M, s.A+s.n-s.M, s.f);
    })}));
  }
  
};

template<class E, class BinPred, class intT>
typename quicksort<E,BinPred,intT>::cfg_type quicksort<E,BinPred,intT>::cfg = quicksort<E,BinPred,intT>::get_cfg();
  
} // end namespace

#endif
