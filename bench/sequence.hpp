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

#include <algorithm>

#include "edsl.hpp"
#include "utils.h"
#include "cmdline.hpp"

#include <limits.h>

#if defined(LONG)
typedef long intT;
typedef unsigned long uintT;
#define INT_T_MAX LONG_MAX
#else
typedef int intT;
typedef int uintT;
#define INT_T_MAX INT_MAX
#endif

#ifndef _HEARTBEAT_PBBS_SEQUENCE_H_
#define _HEARTBEAT_PBBS_SEQUENCE_H_

using plt_type = heartbeat::edsl::pcfg::cactus::parent_link_type;

template <class T>
T* malloc_array(size_t n) {
  return (T*)malloc(n * sizeof(T));
}

template <class T>
struct _seq {
  T* A;
  long n;
  _seq() {A = NULL; n=0;}
  _seq(T* _A, long _n) : A(_A), n(_n) {}
  void del() {free(A);}
};

namespace orig {
  namespace sequence {
    
    template <class intT>
    struct boolGetA {
      bool* A;
      boolGetA(bool* AA) : A(AA) {}
      intT operator() (intT i) {return (intT) A[i];}
    };
    
    template <class ET, class intT>
    struct getA {
      ET* A;
      getA(ET* AA) : A(AA) {}
      ET operator() (intT i) {return A[i];}
    };
    
    template <class IT, class OT, class intT, class F>
    struct getAF {
      IT* A;
      F f;
      getAF(IT* AA, F ff) : A(AA), f(ff) {}
      OT operator () (intT i) {return f(A[i]);}
    };
    
    template <class OT, class intT, class F, class G>
    OT reduceSerial(intT s, intT e, F f, G g) {
      OT r = g(s);
      for (intT j=s+1; j < e; j++) r = f(r,g(j));
      return r;
    }
    
    template <class ET, class intT, class F, class G>
    intT maxIndexSerial(intT s, intT e, F f, G g) {
      ET r = g(s);
      intT k = 0;
      for (intT j=s+1; j < e; j++) {
        ET v = g(j);
        if (f(v,r)) { r = v; k = j;}
      }
      return k;
    }
    
    template <class ET, class intT, class F, class G>
    ET scanSerial(ET* Out, intT s, intT e, F f, G g, ET zero, bool inclusive, bool back) {
      ET r = zero;
      
      if (inclusive) {
        if (back) for (intT i = e-1; i >= s; i--) Out[i] = r = f(r,g(i));
        else for (intT i = s; i < e; i++) Out[i] = r = f(r,g(i));
      } else {
        if (back)
          for (intT i = e-1; i >= s; i--) {
            ET t = g(i);
            Out[i] = r;
            r = f(r,t);
          }
        else
          for (intT i = s; i < e; i++) {
            ET t = g(i);
            Out[i] = r;
            r = f(r,t);
          }
      }
      return r;
    }
    
    template <class ET, class intT, class F>
    ET scanSerial(ET *In, ET* Out, intT n, F f, ET zero) {
      return scanSerial(Out, (intT) 0, n, f, getA<ET,intT>(In), zero, false, false);
    }
    
  }
}

namespace sequence {
 
int block_size = 2048;
  
template <class intT>
struct boolGetA {
  bool* A;
  boolGetA(bool* AA) : A(AA) {}
  intT operator() (intT i) {return (intT) A[i];}
};

template <class ET, class intT>
struct getA {
  ET* A;
  getA(ET* AA) : A(AA) {}
  ET operator() (intT i) {return A[i];}
};

template <class IT, class OT, class intT, class F>
struct getAF {
  IT* A;
  F f;
  getAF(IT* AA, F ff) : A(AA), f(ff) {}
  OT operator () (intT i) {return f(A[i]);}
};
  
template <class intT>
static inline
int nb_blocks(intT n, intT bsize) {
  return 1 + ((n - 1) / bsize);
}
  
template <class intT>
class range_type {
public:
  intT lo;
  intT hi;
  range_type(intT lo, intT hi)
  : lo(lo), hi(hi) { }
};

template <class intT>
range_type<intT> get_rng(intT block_size, intT item_lo, intT item_hi, intT block_lo, intT block_hi) {
  assert(block_lo <= block_hi);
  assert(item_lo <= item_hi);
  intT block_lo2 = block_lo + (item_lo * block_size);
  intT block_hi2 = std::min(block_lo2 + block_size, std::min(item_hi * block_size, block_hi));
  return range_type<intT>(block_lo2, block_hi2);
}
  
template <class OT, class intT, class F, class G>
class reduceSerial : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  intT s; intT e; F f; G g; OT* dest;
  
  reduceSerial(intT s, intT e, F f, G g, OT* dest)
  : s(s), e(e), f(f), g(g), dest(dest) { }
  
  heartbeat_dc_declare(heartbeat::edsl, reduceSerial, sar, par, dc, get_dc)

  static
  dc get_dc() {
    return dc::sequential_loop([] (sar& s, par&) {
        *s.dest = s.g(s.s);
        s.s++;
      }, [] (sar& s, par&) {
        return std::make_pair(&s.s, &s.e);
      }, [] (sar& s, par&, int lo, int hi) {
        *s.dest = s.f(*s.dest, orig::sequence::reduceSerial<OT>(lo, hi, s.f, s.g));
      });
  }
  
};
  
template <class OT, class intT, class F, class G>
typename reduceSerial<OT,intT,F,G>::cfg_type reduceSerial<OT,intT,F,G>::cfg = reduceSerial<OT,intT,F,G>::get_cfg();
  
#ifndef HEARTBEAT_SEQUENCE_USE_PBBS_VERSIONS
  
template <class OT, class intT, class F, class G>
class reduce : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  intT s; intT e; F f; G g; OT* dest;
  
  reduce(intT s, intT e, F f, G g, OT* dest)
  : s(s), e(e), f(f), g(g), dest(dest) { }
  
  class private_activation_record
  : public heartbeat::edsl::pcfg::parallel_loop_private_activation_record<reduce,
  private_activation_record> {
  public:
    
    private_activation_record() {
      private_activation_record::initialize_descriptors();
    }
    
    heartbeat::edsl::pcfg::parallel_combine_activation_record _ar;
    heartbeat::edsl::pcfg::parallel_loop_activation_record* _heartbeat_loop_activation_record_of(heartbeat::edsl::pcfg::parallel_loop_id_type id) {
      return &_ar;
    }
    
    int s; int e; OT acc;
    
  };
  
  heartbeat_dc_loop_declare(heartbeat::edsl, reduce, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::parallel_combine_loop([] (sar& s, par& p) { p.s = s.s; p.e = s.e; },
                                [] (par& p) { return std::make_pair(&p.s, &p.e); },
                                [] (sar&, par& p) { p.acc = OT(); },
                                [] (sar& s, par& p, par& dest) { dest.acc = s.f(p.acc, dest.acc); },
                                [] (sar& s, par& p, int lo, int hi) {
        p.acc = s.f(p.acc, orig::sequence::reduceSerial<OT>(lo, hi, s.f, s.g));
      }),
      dc::stmt([] (sar& s, par& p) {
        *s.dest = p.acc;
      })
    });
  }
  
};
  
#else

template <class OT, class intT, class F, class G>
class reduce : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  intT s; intT e; F f; G g; OT* dest;
  intT l; OT* Sums;
  
  reduce(intT s, intT e, F f, G g, OT* dest)
  : s(s), e(e), f(f), g(g), dest(dest) { }
  
  heartbeat_private_activation_record_begin(heartbeat::edsl, reduce, 1)
    int s; int e;
  heartbeat_private_activation_record_end(heartbeat::edsl, reduce, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par&) {
        s.l = nb_blocks(s.e - s.s, block_size);
      }),
      dc::mk_if([] (sar& s, par&) { return s.l <= 1; }, dc::stmts({
        dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
          return heartbeat_call<reduceSerial<OT,intT,F,G>>(st, pt, s.s, s.e, s.f, s.g, s.dest);
        }),
        dc::exit_function()
      })),
      dc::stmt([] (sar& s, par& p) {
        s.Sums = malloc_array<OT>(s.l);
        p.s = 0;
        p.e = s.l;
      }),
      dc::parallel_for_loop([] (sar&, par& p) { return p.s < p.e; },
                            [] (par& p) { return std::make_pair(&p.s, &p.e); },
                            dc::stmts({
        dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
          auto rng = get_rng(block_size, p.s, p.e, s.s, s.e);
          return heartbeat_call<reduceSerial<OT,intT,F,G>>(st, pt, rng.lo, rng.hi, s.f, s.g, &s.Sums[p.s]);
        }),
        dc::stmt([] (sar& s, par& p) {
          p.s++;
        })
      })),
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        auto g = getA<OT,intT>(s.Sums);
        return heartbeat_call<reduce<OT,intT,F,typeof(g)>>(st, pt, 0, s.l, s.f, g, s.dest);
      }),
      dc::stmt([] (sar& s, par&) {
        free(s.Sums);
      })
    });
  }
  
};
  
#endif
  
template <class OT, class intT, class F, class G>
typename reduce<OT,intT,F,G>::cfg_type reduce<OT,intT,F,G>::cfg = reduce<OT,intT,F,G>::get_cfg();
  
using stack_type = heartbeat::edsl::pcfg::stack_type;

template <class Activation_record, class ...Args>
static stack_type heartbeat_call(stack_type s, plt_type p, Args... args) {
  return heartbeat::edsl::pcfg::push_call<Activation_record>(s, p, args...);
}

template <class OT, class intT, class F>
stack_type reduce4(stack_type s, plt_type p, OT* A, intT n, F f, OT* dest) {
  auto g = getA<OT,intT>(A);
  return heartbeat_call<reduce<OT,intT,F,typeof(g)>>(s, p, (intT)0,n,f,g,dest);
}
  
template <class OT, class intT>
stack_type plusReduce(stack_type s, plt_type p, OT* A, intT n, OT* dest) {
  auto f = pbbs::utils::addF<OT>();
  auto g = getA<OT,intT>(A);
  return heartbeat_call<reduce<OT,intT,typeof(f),typeof(g)>>(s, p, (intT)0,n,f,g,dest);
}

// g is the map function (applied to each element)
// f is the reduce function
template <class OT, class IT, class intT, class F, class G>
stack_type mapReduce(stack_type s, plt_type p, IT* A, intT n, F f, G g, OT* dest) {
  auto g2 = getAF<IT,OT,intT,G>(A,g);
  return heartbeat_call<reduce<OT,intT,F,decltype(g2)>>(s, p, (intT)0,n,f,g2,dest);
}

template <class ET, class intT, class F, class G>
class maxIndexSerial : public heartbeat::edsl::pcfg::shared_activation_record {
public:

  intT s; intT e; F f; G g; intT* dest;
  
  maxIndexSerial(intT s, intT e, F f, G g, intT* dest)
  : s(s), e(e), f(f), g(g), dest(dest) { }
  
  heartbeat_dc_declare(heartbeat::edsl, maxIndexSerial, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::sequential_loop([] (sar& s, par&) {
        *s.dest = 0;
      }, [] (sar& s, par&) {
        return std::make_pair(&s.s, &s.e);
      }, [] (sar& s, par&, int lo, int hi) {
        intT m = orig::sequence::maxIndexSerial<ET>(lo, hi, s.f, s.g);
        if (s.f(s.g(m), s.g(*s.dest))) {
          *s.dest = m;
        }
      });
  }
  
};
  
template <class ET, class intT, class F, class G>
typename maxIndexSerial<ET,intT,F,G>::cfg_type maxIndexSerial<ET,intT,F,G>::cfg = maxIndexSerial<ET,intT,F,G>::get_cfg();

#ifndef HEARTBEAT_SEQUENCE_USE_PBBS_VERSIONS
  
template <class ET, class intT, class F, class G>
class maxIndex : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  intT s; intT e; F f; G g; intT* dest;
  
  maxIndex(intT s, intT e, F f, G g, intT* dest)
  : s(s), e(e), f(f), g(g), dest(dest) { }
  
  class private_activation_record
  : public heartbeat::edsl::pcfg::parallel_loop_private_activation_record<maxIndex,
  private_activation_record> {
  public:
    
    private_activation_record() {
      private_activation_record::initialize_descriptors();
    }
    
    heartbeat::edsl::pcfg::parallel_combine_activation_record _ar;
    heartbeat::edsl::pcfg::parallel_loop_activation_record* _heartbeat_loop_activation_record_of(heartbeat::edsl::pcfg::parallel_loop_id_type id) {
      return &_ar;
    }
    
    int s; int e; ET r; intT k;
    
  };
  
  heartbeat_dc_loop_declare(heartbeat::edsl, maxIndex, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::parallel_combine_loop([] (sar& s, par& p) { p.s = s.s + 1; p.e = s.e; },
                                [] (par& p) { return std::make_pair(&p.s, &p.e); },
                                [] (sar& s, par& p) { p.k = s.s; p.r = s.g(p.k); },
                                [] (sar& s, par& p, par& dest) {
                                  if (s.f(s.g(p.k), s.g(dest.k))) {
                                    dest.r = p.r;
                                    dest.k = p.k;
                                  }
                                },
                                [] (sar& s, par& p, int lo, int hi) {
        auto g = s.g;
        auto f = s.f;
        auto k = p.k;
        ET r = p.r;
        for (auto j = lo; j != hi; j++) {
          ET v = g(j);
          if (f(v, r)) {
            r = v;
            k = j;
          }
        }
        p.r = r;
        p.k = k;
      }),
      dc::stmt([] (sar& s, par& p) {
        *s.dest = p.k;
      })
    });
  }
  
};
  
#else
  
template <class ET, class intT, class F, class G>
class maxIndex : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  intT s; intT e; F f; G g; intT* dest;
  intT l; intT* Idx; intT j;
  
  maxIndex(intT s, intT e, F f, G g, intT* dest)
  : s(s), e(e), f(f), g(g), dest(dest) { }
  
  heartbeat_private_activation_record_begin(heartbeat::edsl, maxIndex, 1)
  int s; int e;
  heartbeat_private_activation_record_end(heartbeat::edsl, maxIndex, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par&) {
        s.l = nb_blocks(s.e - s.s, block_size);
      }),
      dc::mk_if([] (sar& s, par&) { return s.l <= 1; }, dc::stmts({
        dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
          return heartbeat_call<maxIndexSerial<ET, intT, F, G>>(st, pt, s.s, s.e, s.f, s.g, s.dest);
        }),
        dc::exit_function()
      })),
      dc::stmt([] (sar& s, par& p) {
        s.Idx = malloc_array<intT>(s.l);
        p.s = 0;
        p.e = s.l;
      }),
      dc::parallel_for_loop([] (sar&, par& p) { return p.s != p.e; },
                            [] (par& p) { return std::make_pair(&p.s, &p.e); },
                            dc::stmts({
        dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
          auto rng = get_rng(block_size, p.s, p.e, s.s, s.e);
          return heartbeat_call<maxIndexSerial<ET, intT, F, G>>(st, pt, rng.lo, rng.hi, s.f, s.g, &s.Idx[p.s]);
        }),
        dc::stmt([] (sar& s, par& p) {
          p.s++;
        })
      })),
      dc::stmt([] (sar& s, par&) {
        *s.dest = s.Idx[0];
        s.j = 1;
      }),
      dc::sequential_loop([] (sar& s, par&) { return s.j < s.l; }, dc::stmt([] (sar& s, par&) {
        if (s.f(s.g(s.Idx[s.j]), s.g(*s.dest))) {
          *s.dest = s.Idx[s.j];
        }
        s.j++;
      })),
      dc::stmt([] (sar& s, par&) {
        free(s.Idx);
      })
    });
  }
  
};
  
#endif
  
template <class ET, class intT, class F, class G>
typename maxIndex<ET,intT,F,G>::cfg_type maxIndex<ET,intT,F,G>::cfg = maxIndex<ET,intT,F,G>::get_cfg();

template <class ET, class intT, class F>
stack_type maxIndex4(stack_type st, plt_type pt, ET* A, intT n, F f, intT* dest) {
  auto g = getA<ET,intT>(A);
  return heartbeat_call<maxIndex<ET,intT,F,typeof(g)>>(st, pt, (intT)0, n, f, g, dest);
}
  
template <class ET, class intT, class F, class G>
class scanSerial : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  ET* Out; intT s; intT e; F f; G g; ET zero; bool inclusive; bool back; ET* dest;
  intT i;
  
  scanSerial(ET* Out, intT s, intT e, F f, G g, ET zero, bool inclusive, bool back, ET* dest)
  : Out(Out), s(s), e(e), f(f), g(g), zero(zero), inclusive(inclusive), back(back), dest(dest) { }
  
  heartbeat_dc_declare(heartbeat::edsl, scanSerial, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par&) {
        *s.dest = s.zero;
      }),
      dc::mk_if([] (sar& s, par&) { return s.inclusive; },
        dc::mk_if([] (sar& s, par&) { return s.back; },
          dc::sequential_loop([] (sar& s, par&) {
            // no initialization
          }, [] (sar& s, par&) {
            return std::make_pair(&s.s, &s.e);
          }, [] (sar& s, par&, int lo, int hi) {
            auto f = s.f;
            auto g = s.g;
            ET r = *s.dest;
            ET* Out = s.Out;
            for (auto _i = hi; _i != lo; _i--) {
              auto i = _i - 1;
              Out[i] = r = f(r, g(i));
            }
            *s.dest = r;
          }, dc::backward_loop),
          dc::sequential_loop([] (sar& s, par&) {
            // no initialization
          }, [] (sar& s, par&) {
            return std::make_pair(&s.s, &s.e);
          }, [] (sar& s, par&, int lo, int hi) {
            auto f = s.f;
            auto g = s.g;
            ET r = *s.dest;
            ET* Out = s.Out;
            for (auto i = lo; i != hi; i++) {
              Out[i] = r = f(r, g(i));
            }
            *s.dest = r;
          }, dc::forward_loop)),
        dc::mk_if([] (sar& s, par&) { return s.back; },
          dc::sequential_loop([] (sar& s, par&) {
            // no initialization
          }, [] (sar& s, par&) {
            return std::make_pair(&s.s, &s.e);
          }, [] (sar& s, par&, int lo, int hi) {
            auto f = s.f;
            auto g = s.g;
            ET r = *s.dest;
            ET* Out = s.Out;
            for (auto _i = hi; _i != lo; _i--) {
              auto i = _i - 1;
              ET t = g(i);
              Out[i] = r;
              r = f(r, t);
            }
            *s.dest = r;
          }, dc::backward_loop),
          dc::sequential_loop([] (sar& s, par&) {
            // no initialization
          }, [] (sar& s, par&) {
            return std::make_pair(&s.s, &s.e);
          }, [] (sar& s, par&, int lo, int hi) {
            auto f = s.f;
            auto g = s.g;
            ET r = *s.dest;
            ET* Out = s.Out;
            for (auto i = lo; i != hi; i++) {
              ET t = g(i);
              Out[i] = r;
              r = f(r, t);
            }
            *s.dest = r;
          }, dc::forward_loop)))
    });
  }
  
};
  
template <class ET, class intT, class F, class G>
typename scanSerial<ET,intT,F,G>::cfg_type scanSerial<ET,intT,F,G>::cfg = scanSerial<ET,intT,F,G>::get_cfg();
  
template <class ET, class intT, class F>
stack_type scanSerial6(stack_type st, plt_type pt, ET *In, ET* Out, intT n, F f, ET zero, ET* dest) {
  auto g = getA<ET,intT>(In);
  return heartbeat_call<scanSerial<ET,intT,F,typeof(g)>>(st, pt, Out, (intT) 0, n, f, g, zero, false, false, dest);
}

template <class ET, class intT, class F, class G>
class scan : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  ET* Out; intT s; intT e; F f; G g; ET zero; bool inclusive; bool back; ET* dest;
  intT i; intT l; ET* Sums;
  
  scan(ET* Out, intT s, intT e, F f, G g, ET zero, bool inclusive, bool back, ET* dest)
  : Out(Out), s(s), e(e), f(f), g(g), zero(zero), inclusive(inclusive), back(back), dest(dest) { }
  
  heartbeat_private_activation_record_begin(heartbeat::edsl, scan, 2)
    int s; int e;
    ET tmp;
  heartbeat_private_activation_record_end(heartbeat::edsl, scan, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par&) {
        s.l = nb_blocks(s.e - s.s, block_size);
      }),
      dc::mk_if([] (sar& s, par&) { return s.l <= 2; }, dc::stmts({
        dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
          return heartbeat_call<scanSerial<ET,intT,F,G>>(st, pt, s.Out, s.s, s.e, s.f, s.g, s.zero, s.inclusive, s.back, s.dest);
        }),
        dc::exit_function()
      })),
      dc::stmt([] (sar& s, par& p) {
        s.Sums = malloc_array<ET>(s.l);
        p.s = 0;
        p.e = s.l;
      }),
      dc::parallel_for_loop([] (sar&, par& p) { return p.s < p.e; },
                            [] (par& p) { return std::make_pair(&p.s, &p.e); },
                            dc::stmts({
        dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
          auto rng = get_rng(block_size, p.s, p.e, s.s, s.e);
          return heartbeat_call<reduceSerial<ET,intT,F,G>>(st, pt, rng.lo, rng.hi, s.f, s.g, &s.Sums[p.s]);
        }),
        dc::stmt([] (sar& s, par& p) {
          p.s++;
        })
      })),
      dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
        return heartbeat_call<scan>(st, pt, s.Sums, (intT)0, s.l, s.f, getA<ET,intT>(s.Sums), s.zero, false, s.back, s.dest);
      }),
      dc::stmt([] (sar& s, par& p) {
        p.s = 0;
        p.e = s.l;
      }),
      dc::parallel_for_loop([] (sar&, par& p) { return p.s < p.e; },
                            [] (par& p) { return std::make_pair(&p.s, &p.e); },
                            dc::stmts({
        dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
          auto rng = get_rng(block_size, p.s, p.e, s.s, s.e);
          return heartbeat_call<scanSerial<ET,intT,F,G>>(st, pt, s.Out, rng.lo, rng.hi, s.f, s.g, s.Sums[p.s], s.inclusive, s.back, &p.tmp);
        }),
        dc::stmt([] (sar& s, par& p) {
          p.s++;
        })
      })),
      dc::stmt([] (sar& s, par&) {
        free(s.Sums);
      })
    });
  }
  
};

template <class ET, class intT, class F, class G>
typename scan<ET,intT,F,G>::cfg_type scan<ET,intT,F,G>::cfg = scan<ET,intT,F,G>::get_cfg();
  
template <class ET, class intT, class F>
stack_type scan6(stack_type st, plt_type pt, ET *In, ET* Out, intT n, F f, ET zero, ET* dest) {
  auto g = getA<ET,intT>(In);
  return heartbeat_call<scan<ET,intT,F,typeof(g)>>(st, pt, Out, (intT) 0, n, f, g, zero, false, false, dest);
}

template <class ET, class intT, class F>
stack_type scanBack(stack_type st, plt_type pt, ET *In, ET* Out, intT n, F f, ET zero, ET* dest) {
  auto g = getA<ET,intT>(In);
  return heartbeat_call<scan<ET,intT,F,typeof(g)>>(st, pt, Out, (intT) 0, n, f, g, zero, false, true, dest);
}

template <class ET, class intT, class F>
stack_type scanI(stack_type st, plt_type pt, ET *In, ET* Out, intT n, F f, ET zero, ET* dest) {
  auto g = getA<ET,intT>(In);
  return heartbeat_call<scan<ET,intT,F,typeof(g)>>(st, pt, Out, (intT) 0, n, f, g, zero, true, false, dest);
}

template <class ET, class intT, class F>
stack_type scanIBack(stack_type st, plt_type pt, ET *In, ET* Out, intT n, F f, ET zero, ET* dest) {
  auto g = getA<ET,intT>(In);
  return heartbeat_call<scan<ET,intT,F,typeof(g)>>(st, pt, Out, (intT) 0, n, f, g, zero, true, true, dest);
}
  
template <class ET, class intT>
stack_type plusScan(stack_type st, plt_type pt, ET *In, ET* Out, intT n, ET* dest) {
  auto f = pbbs::utils::addF<ET>();
  auto g = getA<ET,intT>(In);
  return heartbeat_call<scan<ET, intT, typeof(f), typeof(g)>>(st, pt, Out, (intT)0, n, f, g, (ET)0, false, false, dest);
}
  
template <class intT>
class sumFlagsSerial : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  bool *Fl; intT n; intT* dest;
  int* IFl; int k; intT r; intT j; 
  int lo; int hi;
  
  sumFlagsSerial(bool *Fl, intT n, intT* dest)
  : Fl(Fl), n(n), dest(dest) { }
  
  heartbeat_dc_declare(heartbeat::edsl, sumFlagsSerial, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par&) {
        s.r = 0;
        s.k = 0;
        s.j = 0;
      }),
      dc::mk_if([] (sar& s, par&) { return (s.n >= 128 && (s.n & 511) == 0 && ((long) s.Fl & 3) == 0); }, dc::stmts({
        dc::stmt([] (sar& s, par&) {
          s.IFl = (int*)s.Fl;
        }),
        dc::sequential_loop([] (sar& s, par&) { s.lo = 0; s.hi = (s.n >> 9); },
                            [] (sar& s, par&) { return std::make_pair(&s.lo, &s.hi); },
                            [] (sar& s, par&, int lo, int hi) {
          auto r = s.r;
          auto IFl = s.IFl;
          for (auto k = lo; k != hi; k++) {
            int rr = 0;
            for (int j=0; j < 128; j++) rr += IFl[j];
            r += (rr&255) + ((rr>>8)&255) + ((rr>>16)&255) + ((rr>>24)&255);
            IFl += 128;
          }
          s.r = r;
          s.IFl = IFl;
        })
      }), // else
      dc::sequential_loop([] (sar& s, par&) { s.j = 0;  },
                          [] (sar& s, par&) { return std::make_pair(&s.j, &s.n); },
                          [] (sar& s, par&, int lo, int hi) {
        intT r = s.r;
        auto Fl = s.Fl;
        for (intT j = lo; j < hi; j++) {
          r += Fl[j];
        }
        s.r = r;
      })),
      dc::stmt([] (sar& s, par&) {
        *s.dest = s.r;
      })
    });
  }
  
};

template <class intT>
typename sumFlagsSerial<intT>::cfg_type sumFlagsSerial<intT>::cfg = sumFlagsSerial<intT>::get_cfg();
  
template <class ET, class intT, class F>
class packSerial : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  ET* Out; bool* Fl; intT s; intT e; F f; _seq<ET>* dest;
  intT m; intT k;
  
  packSerial(ET* Out, bool* Fl, intT s, intT e, F f, _seq<ET>* dest)
  : Out(Out), Fl(Fl), s(s), e(e), f(f), dest(dest) { }
  
  heartbeat_dc_declare(heartbeat::edsl, packSerial, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::mk_if([] (sar& s, par&) { return s.Out == nullptr; }, dc::stmts({
        dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
          return heartbeat_call<sumFlagsSerial<intT>>(st, pt, s.Fl + s.s, s.e - s.s, &s.m);
        }),
        dc::stmt([] (sar& s, par&) {
          s.Out = malloc_array<ET>(s.m);
        })
      })),
      dc::sequential_loop([] (sar& s, par&) {
        s.k = 0;
      }, [] (sar& s, par&) {
        return std::make_pair(&s.s, &s.e);
      }, [] (sar& s, par&, int lo, int hi) {
        auto Fl = s.Fl;
        auto Out = s.Out;
        auto f = s.f;
        auto kk = s.k;
        for (auto i = lo; i != hi; i++) {
          if (Fl[i]) {
            Out[kk++] = f(i);
          }
        }
        s.k = kk;
      }),
      dc::stmt([] (sar& s, par&) {
        *s.dest = _seq<ET>(s.Out,s.k);
      })
    });
  }
  
};
  
template <class ET, class intT, class F>
typename packSerial<ET,intT,F>::cfg_type packSerial<ET,intT,F>::cfg = packSerial<ET,intT,F>::get_cfg();
  
template <class ET, class intT, class F>
class pack : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  ET* Out; bool* Fl; intT s; intT e; F f; _seq<ET>* dest;
  intT l; intT m; intT* Sums;
  
  pack(ET* Out, bool* Fl, intT s, intT e, F f, _seq<ET>* dest)
  : Out(Out), Fl(Fl), s(s), e(e), f(f), dest(dest) { }
  
  heartbeat_private_activation_record_begin(heartbeat::edsl, pack, 2)
    int s; int e;
    _seq<ET> tmp;
  heartbeat_private_activation_record_end(heartbeat::edsl, pack, sar, par, dc, get_dc)
  
  static inline
  int get_pack_block_size() {
    return 2 * block_size;
  }
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par&) {
        s.l = nb_blocks(s.e - s.s, get_pack_block_size());
      }),
      dc::mk_if([] (sar& s, par&) { return s.l <= 1; }, dc::stmts({
        dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
          return heartbeat_call<packSerial<ET, intT, F>>(st, pt, s.Out, s.Fl, s.s, s.e, s.f, s.dest);
        }),
        dc::exit_function()
      })),
      dc::stmt([] (sar& s, par& p) {
        s.Sums = malloc_array<intT>(s.l);
        p.s = 0;
        p.e = s.l;
      }),
      dc::parallel_for_loop([] (sar&, par& p) { return p.s < p.e; },
                            [] (par& p) { return std::make_pair(&p.s, &p.e); },
                            dc::stmts({
        dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
          auto rng = get_rng(get_pack_block_size(), p.s, p.e, s.s, s.e);
          return heartbeat_call<sumFlagsSerial<intT>>(st, pt, s.Fl + rng.lo, rng.hi - rng.lo, &s.Sums[p.s]);
        }),
        dc::stmt([] (sar& s, par& p) {
          p.s++;
        })
      })),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        return plusScan(st, pt, s.Sums, s.Sums, s.l, &s.m);
      }),
      dc::stmt([] (sar& s, par& p) {
        if (s.Out == nullptr) {
          s.Out = malloc_array<ET>(s.m);
        }
        p.s = 0;
        p.e = s.l;
      }),
      dc::parallel_for_loop([] (sar&, par& p) { return p.s < p.e; },
                            [] (par& p) { return std::make_pair(&p.s, &p.e); },
                            dc::stmts({
        dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
          auto rng = get_rng(get_pack_block_size(), p.s, p.e, s.s, s.e);
          return heartbeat_call<packSerial<ET, intT, F>>(st, pt, s.Out + s.Sums[p.s], s.Fl, rng.lo, rng.hi, s.f, &p.tmp);
        }),
        dc::stmt([] (sar& s, par& p) {
          p.s++;
        })
      })),
      dc::stmt([] (sar& s, par& p) {
        free(s.Sums);
        *s.dest = _seq<ET>(s.Out,s.m);
      }) 
    });
  }
  
};
  
template <class ET, class intT, class F>
typename pack<ET,intT,F>::cfg_type pack<ET,intT,F>::cfg = pack<ET,intT,F>::get_cfg();
  
template <class ET, class intT>
stack_type pack5(stack_type st, plt_type pt, ET* In, ET* Out, bool* Fl, intT n, _seq<ET>* dest) {
  auto f = getA<ET,intT>(In);
  return heartbeat_call<pack<ET,intT,typeof(f)>>(st, pt, Out, Fl, (intT) 0, n, f, dest);
}

template <class ET, class intT>
stack_type pack4(stack_type st, plt_type pt, ET* In, bool* Fl, intT n, _seq<ET>* dest) {
  auto f = getA<ET,intT>(In);
  return heartbeat_call<pack<ET,intT,typeof(f)>>(st, pt, (ET*)nullptr, Fl, (intT)0, n, f, dest);
}

template <class intT>
stack_type packIndex(stack_type st, plt_type pt, intT* Out, bool* Fl, intT n, _seq<intT>* dest) {
  auto f = pbbs::utils::identityF<intT>();
  return heartbeat_call<pack<intT,intT,decltype(f)>>(st, pt, Out, Fl, (intT) 0, n, f, dest);
}

template <class intT>
stack_type packIndex(stack_type st, plt_type pt, bool* Fl, intT n, _seq<intT>* dest) {
  auto f = pbbs::utils::identityF<intT>();
  return heartbeat_call<pack<intT,intT,decltype(f)>>(st, pt, (intT *) NULL, Fl, (intT) 0, n, f, dest);
}

template <class ET, class intT, class PRED>
class filterDPS : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  ET* In; ET* Out; intT n; PRED p; intT* dest;
  _seq<ET> tmp; intT m; bool* Fl;
  
  filterDPS(ET* In, ET* Out, intT n, PRED p, intT* dest)
  : In(In), Out(Out), n(n), p(p), dest(dest) { }
  
  heartbeat_private_activation_record_begin(heartbeat::edsl, filterDPS, 1)
    int s; int e;
  heartbeat_private_activation_record_end(heartbeat::edsl, filterDPS, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::parallel_for_loop([] (sar& s, par& p) {
        s.Fl = malloc_array<bool>(s.n);
        p.s = 0;
        p.e = s.n;
      }, [] (par& p) { return std::make_pair(&p.s, &p.e);
      }, [] (sar& s, par& p, int lo, int hi) {
        auto Fl = s.Fl;
        auto pr = s.p;
        auto In = s.In;
        for (auto i = lo; i != hi; i++) {
          Fl[i] = (bool)pr(In[i]);
        }
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        return pack5(st, pt, s.In, s.Out, s.Fl, s.n, &s.tmp);
      }),
      dc::stmt([] (sar& s, par& p) {
        s.m = (intT)s.tmp.n;
        free(s.Fl);
        *s.dest = s.m;
      })
    });
  }
  
};

template <class ET, class intT, class PRED>
typename filterDPS<ET,intT,PRED>::cfg_type filterDPS<ET,intT,PRED>::cfg = filterDPS<ET,intT,PRED>::get_cfg();

template <class ET, class intT, class PRED>
stack_type filterDPS5(stack_type st, plt_type pt, ET* In, ET* Out, intT n, PRED p, intT* dest) {
  return heartbeat_call<filterDPS<ET,intT,PRED>>(st, pt, In, Out, n, p, dest);
}
  
template <class ET, class intT, class PRED>
class filter : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  ET* In; intT n; PRED p; _seq<ET>* dest;
  bool* Fl;
  
  filter(ET* In, intT n, PRED p, _seq<ET>* dest)
  : In(In), n(n), p(p), dest(dest) { }
  
  heartbeat_private_activation_record_begin(heartbeat::edsl, filter, 1)
  int s; int e;
  heartbeat_private_activation_record_end(heartbeat::edsl, filter, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::parallel_for_loop([] (sar& s, par& p) {
        s.Fl = malloc_array<bool>(s.n);
        p.s = 0;
        p.e = s.n;
      }, [] (par& p) { return std::make_pair(&p.s, &p.e);
      }, [] (sar& s, par& p, int lo, int hi) {
        auto Fl = s.Fl;
        auto pr = s.p;
        auto In = s.In;
        for (auto i = lo; i != hi; i++) {
          Fl[i] = (bool)pr(In[i]);
        }
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        return pack4<ET,intT>(st, pt, s.In, s.Fl, s.n, s.dest);
      }),
      dc::stmt([] (sar& s, par& p) {
        free(s.Fl);
      })
    });
  }
  
};

template <class ET, class intT, class PRED>
typename filter<ET,intT,PRED>::cfg_type filter<ET,intT,PRED>::cfg = filter<ET,intT,PRED>::get_cfg();

template <class ET, class intT, class PRED>
stack_type filter4(stack_type st, plt_type pt, ET* In, intT n, PRED p, _seq<ET>* dest) {
  return heartbeat_call<filter<ET, intT, PRED>>(st, pt, In, n, p, dest);
}
  
template <class Iter, class Output_iterator>
class copy : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  Iter lo; Iter hi; Output_iterator dst;
    
  copy(Iter lo, Iter hi, Output_iterator dst)
  : lo(lo), hi(hi), dst(dst) { }
  
  heartbeat_private_activation_record_begin(heartbeat::edsl, copy, 1)
    int s; int e;
  heartbeat_private_activation_record_end(heartbeat::edsl, copy, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return 
      dc::parallel_for_loop([] (sar& s, par& p) { p.s = 0; p.e = s.hi - s.lo; },
                            [] (par& p) { return std::make_pair(&p.s, &p.e); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto lo2 = s.lo + lo;
        auto hi2 = s.lo + hi;
        std::copy(lo2, hi2, s.dst + (lo2 - s.lo));
      });
  }
  
};
  
template <class Iter, class Output_iterator>
typename copy<Iter,Output_iterator>::cfg_type copy<Iter,Output_iterator>::cfg = copy<Iter,Output_iterator>::get_cfg();

template <class Iter, class Output_iterator>
stack_type copy3(stack_type st, plt_type pt, Iter lo, Iter hi, Output_iterator dst) {
  return heartbeat_call<copy<Iter, Output_iterator>>(st, pt, lo, hi, dst);
}
  
template <class Iter, class Item>
class fill : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  Iter lo; Iter hi; const Item* dst;
    
  fill(Iter lo, Iter hi, const Item* dst)
  : lo(lo), hi(hi), dst(dst) { }
  
  heartbeat_private_activation_record_begin(heartbeat::edsl, fill, 1)
    int s; int e;
  heartbeat_private_activation_record_end(heartbeat::edsl, fill, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return 
      dc::parallel_for_loop([] (sar& s, par& p) { p.s = 0; p.e = s.hi - s.lo; },
                            [] (par& p) { return std::make_pair(&p.s, &p.e); },
                            [] (sar& s, par& p, int lo, int hi) {
        auto lo2 = s.lo + lo;
        auto hi2 = s.lo + hi;
        auto v = *s.dst;
        std::fill(lo2, hi2, v);
      });
  }
  
};
  
template <class Iter, class Item>
typename fill<Iter,Item>::cfg_type fill<Iter,Item>::cfg = fill<Iter,Item>::get_cfg();

template <class Iter, class Item>
stack_type fill3(stack_type st, plt_type pt, Iter lo, Iter hi, const Item* dst) {
  return heartbeat_call<fill<Iter, Item>>(st, pt, lo, hi, dst);
}

void initialize() {
  block_size = deepsea::cmdline::parse_or_default("block_size", sequence::block_size);
}
  
} // end namespace

#endif /*! _HEARTBEAT_PBBS_SEQUENCE_H_ */
