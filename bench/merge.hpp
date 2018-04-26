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

#include "heartbeat.hpp"
#include "logging.hpp"
#include "sequence.hpp"

#ifndef _HEARTBEAT__PBBS_MERGE_H_
#define _HEARTBEAT__PBBS_MERGE_H_

namespace heartbeatbench {

using stack_type = heartbeat::edsl::pcfg::stack_type;

#define _MERGE_BSIZEE (1 << 12)

template <class ET, class F, class intT> 
int binSearch(ET* S, intT n, ET v, F f) {
  ET* T = S;
  while (n > 0) {
    intT mid = n/2;
    if (f(v,T[mid])) n = mid;
    else {
      n = (n-mid)-1;
      T = T + mid + 1;
    }
  }
  return T-S;
}

template <class ET, class F, class intT> 
class merge : public heartbeat::edsl::pcfg::shared_activation_record {
public:

  ET* S1; intT l1; ET* S2; intT l2; ET* R; F f;
  intT m1, m2; bool not_done = true;
  ET* pR; ET* pS1; ET* pS2;

  using trampoline = enum { copy1, copy2, loop };
  trampoline t = loop;

  merge(ET* S1, intT l1, ET* S2, intT l2, ET* R, F f)
    : S1(S1), l1(l1), S2(S2), l2(l2), R(R), f(f) { }
  
  heartbeat_dc_declare(heartbeat::edsl, merge, sar, par, dc, get_dc)

  class copy0 { };

  static
  dc get_dc() {
    using controller_type = heartbeat::grain::controller<heartbeat::grain::automatic, merge>;
    using controller_copy_type = heartbeat::grain::controller<heartbeat::grain::automatic, copy0>;
    controller_type::set_ppt(__LINE__, __FILE__);
    controller_copy_type::set_ppt(__LINE__, __FILE__);
    return dc::mk_if([] (sar& s, par& p) {
        return (s.l1 + s.l2) > _MERGE_BSIZEE;
      }, dc::mk_if([] (sar& s, par& p) {
          return s.l2 > s.l1;
        }, dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
          return heartbeat_call<merge>(st, pt, s.S2, s.l2, s.S1, s.l1, s.R, s.f);
        }), dc::stmts({
             dc::stmt([] (sar& s, par& p) {
              s.m1 = s.l1 / 2;
              s.m2 = binSearch(s.S2, s.l2, s.S1[s.m1], s.f);
             }),
             dc::spawn2_join(
               [] (sar& s, par&, plt pt, stt st) {
                 return heartbeat_call<merge>(st, pt, s.S1, s.m1, s.S2, s.m2, s.R, s.f);
               },
               [] (sar& s, par&, plt pt, stt st) {
                 return heartbeat_call<merge>(st, pt,
                              s.S1 + s.m1, s.l1 - s.m1, s.S2 + s.m2, s.l2 - s.m2,
                                                      s.R + s.m1 + s.m2, s.f);
             })
      })),
      dc::stmts({
        dc::stmt([] (sar& s, par& p) {
          s.pR = s.R; 
          s.pS1 = s.S1; 
          s.pS2 = s.S2;
        }),
        dc::sequential_loop([] (sar& s, par& p) { return s.not_done; }, dc::stmt([] (sar& s, par& ep) {
          int nb;
          ET* pR = s.pR; 
          ET* pS1 = s.pS1; 
          ET* pS2 = s.pS2;
          ET* eS1 = s.S1 + s.l1; 
          ET* eS2 = s.S2 + s.l2;
          auto f = s.f;
          trampoline t = s.t;
          int lg_lt;
          while (true) {
            switch (t) {
              case copy1: {
                lg_lt = controller_copy_type::predict_lg_nb_iterations();
                auto lt = controller_copy_type::predict_nb_iterations(lg_lt);
                int lst2 = std::min((int)(lt + (pS2 - s.S2)), s.l2);
                ET* eeS2 = s.S2 + lst2;
                std::copy(pS2, eeS2, pR);
                auto nnb = eeS2 - pS2;
                pS2 += nnb;
                pR += nnb;
                if (pS2 == eS2) {
                  s.not_done = false;
                  return;
                } else {
                  nb = nnb;
                  goto exit_copy;
                }
                break;
              }
              case copy2: {
                lg_lt = controller_copy_type::predict_lg_nb_iterations();
                auto lt = controller_copy_type::predict_nb_iterations(lg_lt);
                int lst1 = std::min((int)(lt + (pS1 - s.S1)), s.l1);
                ET* eeS1 = s.S1 + lst1;
                std::copy(pS1, eeS1, pR);
                auto nnb = eeS1 - pS1;
                pS1 += nnb;
                pR += nnb;
                if (pS1 == eS1) {
                  s.not_done = false;
                  return;
                } else {
                  nb = nnb;
                  goto exit_copy;
                }
                break;
              }
              case loop: {
                lg_lt = controller_type::predict_lg_nb_iterations();
                auto lt = std::max(1, controller_type::predict_nb_iterations(lg_lt));
                int lst1 = std::min((int)(lt + (pS1 - s.S1)), s.l1);
                ET* eeS1 = s.S1 + lst1;
                int lst2 = std::min((int)(lt + (pS2 - s.S2)), s.l2);
                ET* eeS2 = s.S2 + lst2;
                while (true) {
                  if (pS1==eeS1) {break;}
                  if (pS2==eeS2) {break;}
                  *pR++ = f(*pS2,*pS1) ? *pS2++ : *pS1++;
                }
                if (pS1 == eS1) {
                  t = copy1;
                } else if (pS2 == eS2) {
                  t = copy2; 
                } else {
                  nb = (eeS1 - pS1) + (eeS2 - pS2);
                  goto exit;
                }
                break;
              }
            }
          }
         exit:
          s.pR = pR;
          s.pS1 = pS1;
          s.pS2 = pS2;
          s.t = t;
          controller_type::register_callback(lg_lt, nb);
          return;
        exit_copy:
          s.pR = pR;
          s.pS1 = pS1;
          s.pS2 = pS2;
          s.t = t;
          controller_copy_type::register_callback(lg_lt, nb);
          return;
        }))
      }));
  }
  
};

template <class ET, class F, class intT> 
typename merge<ET,F,intT>::cfg_type merge<ET,F,intT>::cfg = merge<ET,F,intT>::get_cfg();

} // end namespace

#endif
