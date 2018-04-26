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

#include <cmath>

#include "transpose.hpp"
#include "quicksort.hpp"

#ifndef _HEARTBEAT__PBBS_SAMPLESORT_H_
#define _HEARTBEAT__PBBS_SAMPLESORT_H_

namespace heartbeatbench {

template<class E, class BinPred, class intT>
void split_positions2(E* a, E* b, intT* c, intT length_a, intT length_b, BinPred compare) {
  if (length_a == 0 || length_b == 0) {
    return;
  }
  int pos_a = 0;
  int pos_b = 0;
  int pos_c = 0;
  for (intT i = 0; i < 2 * length_b; i++) {
    c[i] = 0;
  }
  while (pos_b < length_b) {
    while (pos_a < length_a && compare(a[pos_a], b[pos_b])) {
      c[pos_c]++;
      pos_a++;
    }
    pos_c++;
    while (pos_a < length_a && (compare(a[pos_a], b[pos_b]) ^ compare(b[pos_b], a[pos_a]) ^ true)) {
      c[pos_c]++;
      pos_a++;
    }
    
    pos_b++;
    pos_c++;
    // The pivots are equal
    while (pos_b < length_b && !compare(b[pos_b - 1], b[pos_b])) {
      pos_b++;
      pos_c += 2;
    }
  }
  c[pos_c] = length_a - pos_a;
}

//#define SSORT_THR 100000
int SSORT_THR = 100000;
#define AVG_SEG_SIZE 2
#define PIVOT_QUOT 2

template<class E, class BinPred, class intT>
class sampleSort : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  E* a; intT n; BinPred compare;
  intT sq; intT row_length; intT rows; intT segments;
  int over_sample; intT sample_set_size; int pivots_size; E* sample_set;
  E* pivots; E *b; intT *segments_sizes; intT *offset_a; intT *offset_b;
  
  sampleSort(E* a, intT n, BinPred compare)
  : a(a), n(n), compare(compare) { }
  
  heartbeat_private_activation_record_begin(heartbeat::edsl, sampleSort, 5)
  intT s; intT e;
  intT offset; intT size; intT tmp;
  heartbeat_private_activation_record_end(heartbeat::edsl, sampleSort, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::mk_if([] (sar& s, par&) { return s.n < SSORT_THR; }, dc::stmts({
        dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
          return heartbeat_call<quicksort<E,BinPred,intT>>(st, pt, s.a, s.n, s.compare);
        }),
        dc::exit_function()
      })),
      dc::stmt([] (sar& s, par& p) {
        s.sq = (intT)sqrt(s.n);
        s.row_length = s.sq * AVG_SEG_SIZE;
        s.rows = (intT) ceil(1. * s.n / s.row_length);
        // number of pivots + 1
        s.segments = (s.sq - 1) / PIVOT_QUOT;
      }),
      dc::mk_if([] (sar& s, par&) { return s.segments <= 1; }, dc::stmts({
        dc::stmt([] (sar& s, par&) {
          std::sort(s.a, s.a + s.n, s.compare);
        }),
        dc::exit_function()
      })),
      dc::stmt([] (sar& s, par& p) {
        s.over_sample = 4;
        s.sample_set_size = s.segments * s.over_sample;
        s.sample_set = malloc_array<E>(s.sample_set_size);
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.s = 0; p.e = s.sample_set_size; },
                            [] (par& p) { return std::make_pair(&p.s, &p.e); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto sample_set = s.sample_set;
        auto a = s.a;
        auto n = s.n;
        for (auto i = lo; i != hi; i++) {
          intT o = pbbs::utils::hash(i)% n;
          sample_set[i] = a[o];
        }
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        return heartbeat_call<quicksort<E,BinPred,intT>>(st, pt, s.sample_set, s.sample_set_size, s.compare);
      }),
      dc::stmt([] (sar& s, par& p) {
        s.pivots_size = s.segments - 1;
        s.pivots = malloc_array<E>(s.pivots_size);
      }),
      dc::parallel_for_loop([] (sar& s, par& p) { p.s = 0; p.e = s.pivots_size; },
                            [] (par& p) { return std::make_pair(&p.s, &p.e); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto over_sample = s.over_sample;
        auto sample_set = s.sample_set;
        auto pivots = s.pivots;
        for (auto i = lo; i != hi; i++) {
          intT o = over_sample * i;
          pivots[i] = sample_set[o];
        }
      }),
      dc::stmt([] (sar& s, par& p) {
        free(s.sample_set);
        s.segments = 2 * s.segments - 1;
        s.b = malloc_array<E>(s.rows * s.row_length);
        s.segments_sizes = malloc_array<intT>(s.rows * s.segments);
        s.offset_a = malloc_array<intT>(s.rows * s.segments);
        s.offset_b = malloc_array<intT>(s.rows * s.segments);
        p.s = 0;
        p.e = s.rows;
      }),
      dc::parallel_for_loop([] (sar&, par& p) { return p.s < p.e; },
                            [] (par& p) { return std::make_pair(&p.s, &p.e); },
                            dc::stmts({
        dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
          p.offset = p.s * s.row_length;
          p.size =  (p.s < s.rows - 1) ? s.row_length : s.n - p.offset;
          return heartbeat_call<sampleSort>(st, pt, s.a + p.offset, p.size, s.compare);
        }),
        dc::stmt([] (sar& s, par& p) {
          split_positions2(s.a + p.offset, s.pivots, s.segments_sizes + p.s * s.segments, p.size, s.pivots_size, s.compare);
          p.s++;
        }),
      })),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        auto plus = [] (intT x, intT y) {
          return x + y;
        };
        return sequence::scan6(st, pt, s.segments_sizes, s.offset_a, s.rows*s.segments, plus, 0, &p.tmp);
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        return transpose2(st, pt, s.segments_sizes, s.offset_b, s.rows, s.segments);
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        auto plus = [] (intT x, intT y) {
          return x + y;
        };
        return sequence::scan6(st, pt, s.offset_b, s.offset_b, s.rows*s.segments, plus, 0, &p.tmp);
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        return blockTrans2(st, pt, s.a, s.b, s.offset_a, s.offset_b, s.segments_sizes, s.rows, s.segments);
      }),
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        return heartbeat_call<sequence::copy<E*, E*>>(st, pt, s.b, s.b + s.n, s.a);
      }),
      dc::stmt([] (sar& s, par& p) {
        free(s.b);
        free(s.offset_a);
        free(s.segments_sizes);
        p.s = 0;
        p.e = s.pivots_size + 1;
      }),
      dc::parallel_for_loop([] (sar&, par& p) { return p.s < p.e; },
                            [] (par& p) { return std::make_pair(&p.s, &p.e); },
                            dc::stmts({
        dc::stmt([] (sar& s, par& p) {
          p.offset = s.offset_b[(2 * p.s) * s.rows];
        }),
        dc::cond({
          std::make_pair([] (sar& s, par& p) {
            return (p.s == 0);
          }, dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
            return heartbeat_call<sampleSort<E,BinPred,intT>>(st, pt, s.a, s.offset_b[s.rows], s.compare);
          })),
          std::make_pair([] (sar& s, par& p) {
            return (p.s < s.pivots_size);
          }, dc::mk_if([] (sar& s, par& p) {
            return (s.compare(s.pivots[p.s-1],s.pivots[p.s]));
          }, dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
            return heartbeat_call<sampleSort<E,BinPred,intT>>(st, pt, s.a + p.offset, s.offset_b[(2 * p.s + 1) * s.rows] - p.offset, s.compare);
          })))
        }, dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
          return heartbeat_call<sampleSort<E,BinPred,intT>>(st, pt, s.a + p.offset, s.n - p.offset, s.compare);
        })),
        dc::stmt([] (sar& s, par& p) {
          p.s++;
        })
      })),
      dc::stmt([] (sar& s, par& p) {
        free(s.pivots);
        free(s.offset_b);
      })
    });
  }
  
};

template<class E, class BinPred, class intT>
typename sampleSort<E,BinPred,intT>::cfg_type sampleSort<E,BinPred,intT>::cfg = sampleSort<E,BinPred,intT>::get_cfg();

template<class E, class BinPred, class intT>
stack_type sampleSort3(stack_type st, plt_type pt, E* a, intT n, BinPred compare) {
  return sequence::heartbeat_call<sampleSort<E,BinPred,intT>>(st, pt, a, n, compare);
}

} // end namespace

#endif
