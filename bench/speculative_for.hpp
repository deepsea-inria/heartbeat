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
#include "sequence.hpp"

#ifndef _HEARTBEAT_SPECULATIVE_FOR_H_
#define _HEARTBEAT_SPECULATIVE_FOR_H_

namespace heartbeatbench {

struct reservation {
  intT r;
  reservation() : r(INT_T_MAX) {}
  void reserve(intT i) { pbbs::utils::writeMin(&r, i); }
  bool reserved() { return (r < INT_T_MAX);}
  void reset() {r = INT_T_MAX;}
  bool check(intT i) { return (r == i);}
  bool checkReset(intT i) {
    if (r==i) { r = INT_T_MAX; return 1;}
    else return 0;
  }
};

inline void reserveLoc(intT& x, intT i) {pbbs::utils::writeMin(&x,i);}

template <class S>
class speculative_for : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  S step; intT s; intT e; int granularity;
  bool hasState; int maxTries; intT* dest;
  
  intT maxRoundSize; intT currentRoundSize; intT *I; intT *Ihold; bool *keep;
  S *state; intT size; _seq<intT> tmp;
  
  int round; intT numberDone; intT numberKeep; intT totalProcessed;
    
  speculative_for(S step, intT s, intT e, int granularity,
                  bool hasState, int maxTries, intT* dest)
  : step(step), s(s), e(e), granularity(granularity),
  hasState(hasState), maxTries(maxTries), dest(dest) { }
  
  heartbeat_private_activation_record_begin(heartbeat::edsl, speculative_for, 4)
  int s; int e;
  heartbeat_private_activation_record_end(heartbeat::edsl, speculative_for, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par& p) {
        if (s.maxTries < 0) s.maxTries = 100 + 200*s.granularity;
        s.maxRoundSize = (s.e-s.s)/s.granularity+1;
        s.currentRoundSize = s.maxRoundSize;
        s.I = malloc_array<intT>(s.maxRoundSize);
        s.Ihold = malloc_array<intT>(s.maxRoundSize);
        s.keep = malloc_array<bool>(s.maxRoundSize);
      }),
      dc::mk_if([] (sar& s, par&) { return s.hasState; },
        dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
          s.state = malloc_array<S>(s.maxRoundSize);
          return sequence::fill3(st, pt, s.state, s.state + s.maxRoundSize, &(s.step));
      })),
      dc::stmt([] (sar& s, par& p) {
        s.round = 0;
        s.numberDone = s.s; // number of iterations done
        s.numberKeep = 0; // number of iterations to carry to next round
        s.totalProcessed = 0; // number done including wasteds tries
      }),
      dc::sequential_loop([] (sar& s, par&) { return s.numberDone < s.e; }, dc::stmts({
        dc::stmt([] (sar& s, par& p) {
          if (s.round++ > s.maxTries) {
            std::cout << "speculativeLoop: too many iterations, increase maxTries parameter" << std::endl;
            abort();
          }
          s.size = std::min(s.currentRoundSize, s.e - s.numberDone);
          s.totalProcessed += s.size;
        }),
        dc::mk_if([] (sar& s, par&) { return s.hasState; },
          dc::parallel_for_loop([] (sar& s, par& p) { p.s = 0; p.e = s.size; },
                                [] (par& p) { return std::make_pair(&p.s, &p.e); },
                                [] (sar& s, par& p, int lo, int hi) {
            auto numberKeep = s.numberKeep;
            auto I = s.I;
            auto numberDone = s.numberDone;
            auto keep = s.keep;
            S* state = s.state;
            for (auto i = lo; i != hi; i++) {
              if (i >= numberKeep) {
                I[i] = numberDone + i;
              }
              keep[i] = state[i].reserve(I[i]);
            }
          }),
          dc::parallel_for_loop([] (sar& s, par& p) { p.s = 0; p.e = s.size; },
                                [] (par& p) { return std::make_pair(&p.s, &p.e); },
                                [] (sar& s, par& p, int lo, int hi) {
            auto numberKeep = s.numberKeep;
            auto I = s.I;
            auto numberDone = s.numberDone;
            auto keep = s.keep;
            S step = s.step;
            for (auto i = lo; i != hi; i++) {
              if (i >= numberKeep) {
                I[i] = numberDone + i;
              }
              keep[i] = step.reserve(I[i]);
            }
          })),
        dc::mk_if([] (sar& s, par&) { return s.hasState; },
          dc::parallel_for_loop([] (sar& s, par& p) { p.s = 0; p.e = s.size; },
                                [] (par& p) { return std::make_pair(&p.s, &p.e); },
                                [] (sar& s, par& p, int lo, int hi) {
            auto I = s.I;
            auto keep = s.keep;
            S* state = s.state;
            for (auto i = lo; i != hi; i++) {
              if (keep[i]) {
                keep[i] = !state[i].commit(I[i]);
              }
            }
          }),
          dc::parallel_for_loop([] (sar& s, par& p) { p.s = 0; p.e = s.size; },
                                [] (par& p) { return std::make_pair(&p.s, &p.e); },
                                [] (sar& s, par& p, int lo, int hi) {
            auto I = s.I;
            auto keep = s.keep;
            S step = s.step;
            for (auto i = lo; i != hi; i++) {
              if (keep[i]) {
                keep[i] = !step.commit(I[i]);
              }
            }
          })),
        dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
          return sequence::pack5(st, pt, s.I, s.Ihold, s.keep, s.size, &s.tmp);
        }),
        dc::stmt([] (sar& s, par& p) {
          // keep iterations that failed for next round
          s.numberKeep = (intT)s.tmp.n;
          std::swap(s.I, s.Ihold);
          s.numberDone += s.size - s.numberKeep;
          // adjust round size based on number of failed attempts
          if (float(s.numberKeep)/float(s.size) > .2) {
            s.currentRoundSize = std::max(s.currentRoundSize/2, 
                                          std::max(s.maxRoundSize/64 + 1, s.numberKeep));
          } else if (float(s.numberKeep)/float(s.size) < .1) {
            s.currentRoundSize = std::min(s.currentRoundSize * 2, s.maxRoundSize);
          } 
        })
      })),
      dc::stmt([] (sar& s, par& p) {
        free(s.I); free(s.Ihold); free(s.keep);
        if(s.hasState) free(s.state);
        *s.dest = s.totalProcessed;
      })
    });
  }
  
};

template <class S>
typename speculative_for<S>::cfg_type speculative_for<S>::cfg = speculative_for<S>::get_cfg();

using stack_type = heartbeat::edsl::pcfg::stack_type;

template <class S>
stack_type speculative_for5(stack_type st, plt_type pt, S step, intT s, intT e, int granularity,
                            bool hasState, int maxTries, intT* dest) {
  return sequence::heartbeat_call<speculative_for<S>>(st, pt, step, s, e, granularity, hasState, maxTries, dest);
}

template <class S>
stack_type speculative_for4(stack_type st, plt_type pt, S step, intT s, intT e, int granularity, intT* dest) {
  return sequence::heartbeat_call<speculative_for<S>>(st, pt, step, s, e, granularity, 1, -1, dest);
}

} // end namespace
  
#endif
