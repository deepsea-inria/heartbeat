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
#include <math.h>

#include "sequence.hpp"
#include "transpose.hpp"
#include "utils.h"
#include "logging.hpp"

#ifndef _HEARTBEAT_A_RADIX_INCLUDED
#define _HEARTBEAT_A_RADIX_INCLUDED

namespace heartbeatbench {
  
namespace intSort {
  
  // Cannot be greater than 8 without changing definition of bIndexT
  //    from unsigned char to unsigned int (or unsigned short)
#define MAX_RADIX 8
#define BUCKETS 256    // 1 << MAX_RADIX
  
  //typedef int bucketsT[BUCKETS];
  
  
  // a type that must hold MAX_RADIX bits
  typedef unsigned char bIndexT;
    
  // A is the input and sorted output (length = n)
  // B is temporary space for copying data (length = n)
  // Tmp is temporary space for extracting the bytes (length = n)
  // BK is an array of bucket sets, each set has BUCKETS integers
  //    it is used for temporary space for bucket counts and offsets
  // numBK is the length of BK (number of sets of buckets)
  // the first entry of BK is also used to return the offset of each bucket
  // m is the number of buckets per set (m <= BUCKETS)
  // extract is a function that extract the appropriate bits from A
  //  it must return a non-negative integer less than m
  template <class E, class F, class intT>
  class radixStep : public heartbeat::edsl::pcfg::shared_activation_record {
  public:
    
    E* A; E* B; bIndexT *Tmp; intT (*BK)[BUCKETS];
    intT numBK; intT n; intT m; bool top; F extract;
    
    int expand; intT blocks; intT nn; intT* cnts; intT* oA; intT* oB; intT ss;
    intT jlo; intT jhi;
    
    radixStep(E* A, E* B, bIndexT *Tmp, intT (*BK)[BUCKETS],
              intT numBK, intT n, intT m, bool top, F extract)
    : A(A), B(B), Tmp(Tmp), BK(BK), numBK(numBK), n(n),
    m(m), top(top), extract(extract) { }

    using trampoline = enum { entry, loop1, loop2, loop3, loop4 };
    
    using rbty = struct radixBlock2 {
      E* A; E* B; bIndexT *Tmp; intT* counts; intT* offsets;
      intT Boffset; intT n; intT m; F extract;

      radixBlock2() { }
      radixBlock2(E* A, E* B, bIndexT *Tmp, intT* counts, intT* offsets,
                 intT Boffset, intT n, intT m, F extract) :
        A(A), B(B), Tmp(Tmp), counts(counts), offsets(offsets), Boffset(Boffset), n(n), m(m), extract(extract) { }

      intT i,j,s; trampoline t = entry;
    };
    
    heartbeat_private_activation_record_begin(heartbeat::edsl, radixStep, 1)
    int s; int e;
    intT od; intT nni; bool not_done;
    rbty rb;
    heartbeat_private_activation_record_end(heartbeat::edsl, radixStep, sar, par, dc, get_dc)

    class loop_grain0 { };
    class loop_grain1 { };
    class loop_grain2 { };
    class loop_grain3 { };
    
    static
    dc get_dc() {
      using controller_type0 = heartbeat::grain::controller<heartbeat::grain::automatic, loop_grain0>;
      using controller_type1 = heartbeat::grain::controller<heartbeat::grain::automatic, loop_grain1>;
      using controller_type2 = heartbeat::grain::controller<heartbeat::grain::automatic, loop_grain2>;
      using controller_type3 = heartbeat::grain::controller<heartbeat::grain::automatic, loop_grain3>;
      controller_type0::set_ppt(__LINE__, __FILE__);
      controller_type1::set_ppt(__LINE__, __FILE__);
      controller_type2::set_ppt(__LINE__, __FILE__);
      controller_type3::set_ppt(__LINE__, __FILE__);
      auto radixBlock = [] (sar& , par& p) {
        auto i = p.rb.i; auto j = p.rb.j; auto s = p.rb.s; auto t = p.rb.t;
        auto A = p.rb.A; auto B = p.rb.B; bIndexT *Tmp = p.rb.Tmp; intT* counts = p.rb.counts;
        intT* offsets = p.rb.offsets;
        intT Boffset = p.rb.Boffset; intT n = p.rb.n; intT m = p.rb.m; F extract = p.rb.extract;
        int nbiters = 0;
        switch (t) {
          case entry: {
            i = 0;
            t = loop1;
          }
          case loop1: {
            auto lg_lt = controller_type0::predict_lg_nb_iterations();
            auto lt = controller_type0::predict_nb_iterations(lg_lt);
            auto lst = std::min(m, i + lt);
            nbiters = lst - i;
            while (i < lst) {
              counts[i] = 0;
              i++;
            }
            if (i != m) {
              controller_type0::register_callback(lg_lt, nbiters);
              goto exit;
            }
            j = 0;
            t = loop2;
          }
          case loop2: {
            auto lg_lt = 10; //controller_type1::predict_lg_nb_iterations();
            auto lt = controller_type1::predict_nb_iterations(lg_lt);
            auto lst = std::min(n, j + lt);
            nbiters = lst - j;
            while (j < lst) {
              intT k = Tmp[j] = extract(A[j]);
              counts[k]++;
              j++;
            }
            if (j != n) {
              //              controller_type1::register_callback(lg_lt, nbiters);
              goto exit;
            }
            s = Boffset;
            i = 0;
            t = loop3;
          }
          case loop3: {
            auto lg_lt = controller_type2::predict_lg_nb_iterations();
            auto lt = controller_type2::predict_nb_iterations(lg_lt);
            auto lst = std::min(m, i + lt);
            nbiters = lst - i;
            while (i < lst) {
              s += counts[i];
              offsets[i] = s;
              i++;
            }
            if (i != m) {
              controller_type2::register_callback(lg_lt, nbiters);
              goto exit;
            }
            j = n-1;
            t = loop4;
          }
          case loop4: {
            auto lg_lt = 10; //controller_type3::predict_lg_nb_iterations();
            int lt = controller_type3::predict_nb_iterations(lg_lt);
            int lst = std::max(j - lt, 0);
            nbiters = j - lst;
            while (j >= lst) {
              intT x =  --offsets[Tmp[j]];
              B[x] = A[j];
              j--;
            }
            if (j + 1 != 0) {
              //              controller_type3::register_callback(lg_lt, nbiters);
              goto exit;
            }
          }
        }
        p.not_done = false;
        return;
      exit:
        p.rb.t = t; p.rb.i = i; p.rb.j = j; p.rb.s = s;
        p.not_done = true;
        return;
      };
      return dc::stmts({
        dc::stmt([] (sar& s, par& p) {
          // need 3 bucket sets per block
          s.expand = (sizeof(E)<=4) ? 64 : 32;
          s.blocks = std::min(s.numBK/3,(1+s.n/(BUCKETS*s.expand)));
        }),
        dc::mk_if([] (sar& s, par& p) {
          return s.blocks < 2;
        }, dc::stmts({
          dc::stmt([] (sar& s, par& p) {
            p.not_done = true;
            new (&p.rb) rbty(s.A, s.B, s.Tmp, s.BK[0], s.BK[0], (intT) 0, s.n, s.m, s.extract);
          }),
          dc::sequential_loop([] (sar& s, par& p) { return p.not_done; }, dc::stmt(radixBlock)),
          dc::stmt([] (sar& s, par& p) {
            s.ss = 0;
          }),
          dc::sequential_loop([] (sar& s, par&) { return s.ss < s.n; },
                              [] (sar& s, par&) { return std::make_pair(&s.ss, &s.n); },
                              [] (sar& s, par&, int lo, int hi) {
            std::copy(s.B + lo, s.B + hi, s.A + lo);
          }), 
          dc::exit_function()
        })),
        dc::stmt([] (sar& s, par& p) {
          s.nn = (s.n+s.blocks-1)/s.blocks;
          s.cnts = (intT*) s.BK;
          s.oA = (intT*) (s.BK+s.blocks);
          s.oB = (intT*) (s.BK+2*s.blocks);
          p.s = 0;
          p.e = s.blocks;
        }), 
        dc::parallel_for_loop([] (sar&, par& p) { return p.s < p.e; },
                              [] (par& p) { return std::make_pair(&p.s, &p.e); },
                              dc::stmts({
          dc::stmt([] (sar& s, par& p) {
            intT od = p.s*s.nn;
            intT nni = std::min(std::max<intT>(s.n-od,0),s.nn);
            //            radixBlock(s.A+od, s.B, s.Tmp+od, s.cnts + s.m*p.s, s.oB + s.m*p.s, od, nni, s.m, s.extract);
            //            p.s++;
            p.not_done = true;
            new (&p.rb) rbty(s.A+od, s.B, s.Tmp+od, s.cnts + s.m*p.s, s.oB + s.m*p.s, od, nni, s.m, s.extract);
          }),
          dc::sequential_loop([] (sar& s, par& p) { return p.not_done; }, dc::stmt(radixBlock)),  
          dc::stmt([] (sar& s, par& p) {
            p.s++;
          })
        })),
        dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
          return transpose2(st, pt, s.cnts, s.oA, s.blocks, s.m);
        }),
        dc::mk_if([] (sar& s, par& p) {
          return s.top;
        }, dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
          return sequence::scan6(st, pt, s.oA, s.oA, s.blocks*s.m, pbbs::utils::addF<intT>(), 0, &s.ss);
        }), dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
          return sequence::scanSerial6(st, pt, s.oA, s.oA, s.blocks*s.m, pbbs::utils::addF<intT>(), 0, &s.ss);
        })),
        dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
          return blockTrans2(st, pt, s.B, s.A, s.oB, s.oA, s.cnts, s.blocks, s.m);
        }),
        dc::sequential_loop([] (sar& s, par&) {
          s.jlo = 0;
          s.jhi = s.m;
        }, [] (sar& s, par&) {
          return std::make_pair(&s.jlo, &s.jhi);
        }, [] (sar& s, par&, int lo, int hi) {
          // put the offsets for each bucket in the first bucket set of BK
          auto BK = s.BK; auto oA = s.oA; auto blocks = s.blocks;
          for (auto j = lo; j != hi; j++) {
            BK[0][j] = oA[j*blocks];
          }
        })
      });
    }
    
  };
  
  template <class E, class F, class intT>
  typename radixStep<E,F,intT>::cfg_type radixStep<E,F,intT>::cfg = radixStep<E,F,intT>::get_cfg();
  
  // a function to extract "bits" bits starting at bit location "offset"
  template <class E, class F>
  struct eBits {
    F _f;  intT _mask;  intT _offset;
    eBits() { }
    eBits(int bits, intT offset, F f): _mask((1<<bits)-1),
    _offset(offset), _f(f) {}
    intT operator() (E p) {return _mask&(_f(p)>>_offset);}
  };

  // Radix sort with low order bits first
  template <class E, class F, class intT>
  class radixLoopBottomUp : public heartbeat::edsl::pcfg::shared_activation_record {
  public:
    
    E *A; E *B; bIndexT *Tmp; intT (*BK)[BUCKETS];
    intT numBK; intT n; int bits; bool top; F f;
    
    int rounds; int rbits; int bitOffset;
    
    radixLoopBottomUp(E *A, E *B, bIndexT *Tmp, intT (*BK)[BUCKETS],
                      intT numBK, intT n, int bits, bool top, F f)
    : A(A), B(B), Tmp(Tmp), BK(BK), numBK(numBK), n(n), bits(bits), top(top), f(f) { }
    
    heartbeat_dc_declare(heartbeat::edsl, radixLoopBottomUp, sar, par, dc, get_dc)

    static
    dc get_dc() {
      return dc::stmts({
        dc::stmt([] (sar& s, par& p) {
          s.rounds = 1+(s.bits-1)/MAX_RADIX;
          s.rbits = 1+(s.bits-1)/s.rounds;
          s.bitOffset = 0;
        }),
        dc::sequential_loop([] (sar& s, par&) { return s.bitOffset < s.bits; }, dc::stmts({
          dc::stmt([] (sar& s, par& p) {
            if (s.bitOffset+s.rbits > s.bits) s.rbits = s.bits-s.bitOffset;
          }),
          dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
            auto f = eBits<E,F>(s.rbits,s.bitOffset,s.f);
            return heartbeat_call<radixStep<E,typeof(f),intT>>(st, pt, s.A, s.B, s.Tmp, s.BK, s.numBK, s.n, (intT)1 << s.rbits, s.top, f);
            }),
          dc::stmt([] (sar& s, par& p) {
            s.bitOffset += s.rbits;
          })
        }))
      });
    }
    
  };
  
  template <class E, class F, class intT>
  typename radixLoopBottomUp<E,F,intT>::cfg_type radixLoopBottomUp<E,F,intT>::cfg = radixLoopBottomUp<E,F,intT>::get_cfg();
  
  // Radix sort with high order bits first
  template <class E, class F, class intT>
  class radixLoopTopDown : public heartbeat::edsl::pcfg::shared_activation_record {
  public:
    
    E *A; E *B; bIndexT *Tmp; intT (*BK)[BUCKETS];
    intT numBK; intT n; int bits;  F f;
    intT* offsets; intT remain; float y;
    
    radixLoopTopDown(E *A, E *B, bIndexT *Tmp, intT (*BK)[BUCKETS],
                      intT numBK, intT n, int bits, F f)
    : A(A), B(B), Tmp(Tmp), BK(BK), numBK(numBK), n(n), bits(bits), f(f) { }
    
    heartbeat_private_activation_record_begin(heartbeat::edsl, radixLoopTopDown, 1)
      int s; int e;
      intT segOffset; intT segNextOffset; intT segLen; intT blocksOffset;
      intT blocksNextOffset; intT blockLen;
    heartbeat_private_activation_record_end(heartbeat::edsl, radixLoopTopDown, sar, par, dc, get_dc)
    
    static
    dc get_dc() {
      return dc::stmts({
        dc::cond({
          std::make_pair([] (sar& s, par& p) {
            return s.n == 0;
          }, dc::exit_function()),
          std::make_pair([] (sar& s, par& p) {
            return s.bits <= MAX_RADIX;
          }, dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
            auto f = eBits<E,F>(s.bits,0,s.f);
            return heartbeat_call<radixStep<E,typeof(f),intT>>(st, pt, s.A, s.B, s.Tmp, s.BK, s.numBK, s.n, (intT)1 << s.bits, true, f);
          })),
          std::make_pair([] (sar& s, par& p) {
            return s.numBK >= BUCKETS+1;
          }, dc::stmts({
            dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
              auto f = eBits<E,F>(MAX_RADIX,s.bits-MAX_RADIX,s.f);
              return heartbeat_call<radixStep<E,typeof(f),intT>>(st, pt, s.A, s.B, s.Tmp, s.BK, s.numBK, s.n, (intT)BUCKETS, true, f);
            }),
            dc::stmt([] (sar& s, par& p) {
              s.offsets = s.BK[0];
              s.remain = s.numBK - BUCKETS - 1;
              s.y = s.remain / (float) s.n;
              p.s = 0;
              p.e = BUCKETS;
            }),
            dc::parallel_for_loop([] (sar&, par& p) { return p.s < p.e; },
                                  [] (par& p) { return std::make_pair(&p.s, &p.e); },
                                  dc::stmts({
              dc::stmt([] (sar& s, par& p) {
                p.segOffset = s.offsets[p.s];
                p.segNextOffset = (p.s == BUCKETS-1) ? s.n : s.offsets[p.s+1];
                p.segLen = p.segNextOffset - p.segOffset;
                p.blocksOffset = ((intT) floor(p.segOffset * s.y)) + p.s + 1;
                p.blocksNextOffset = ((intT) floor(p.segNextOffset * s.y)) + p.s + 2;
                p.blockLen = p.blocksNextOffset - p.blocksOffset;
              }),
              dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
                return heartbeat_call<radixLoopTopDown<E,F,intT>>(st, pt, s.A + p.segOffset, s.B + p.segOffset, s.Tmp + p.segOffset,
                                                               s.BK + p.blocksOffset, p.blockLen, p.segLen,
                                                               s.bits-MAX_RADIX, s.f);
              }),
              dc::stmt([] (sar& s, par& p) {
                p.s++;
              })
            }))
          }))
            }, dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
                return heartbeat_call<radixLoopBottomUp<E, F, intT>>(st, pt, s.A, s.B, s.Tmp, s.BK, s.numBK, s.n, s.bits, false, s.f);
              }))
          });
    }
    
  };
  
  template <class E, class F, class intT>
  typename radixLoopTopDown<E,F,intT>::cfg_type radixLoopTopDown<E,F,intT>::cfg = radixLoopTopDown<E,F,intT>::get_cfg();
  
  template <class E, class intT>
  long iSortSpace(intT n) {
    typedef intT bucketsT[BUCKETS];
    intT numBK = 1+n/(BUCKETS*8);
    return sizeof(E)*n + sizeof(bIndexT)*n + sizeof(bucketsT)*numBK;
  }
  
  template <class E, class F, class intT>
  class iSort : public heartbeat::edsl::pcfg::shared_activation_record {
  public:
    
    E *A; intT* bucketOffsets; intT n; intT m; bool bottomUp;
    char* tmpSpace; F f;
    typedef intT bucketsT[BUCKETS];
    
    int bits;
    intT numBK;
    
    // the temporary space is broken into 3 parts: B, Tmp and BK
    E *B; intT Bsize; bIndexT *Tmp; intT tmpSize; bucketsT *BK;
    intT tmp;
   
    iSort(E *A, intT* bucketOffsets, intT n, intT m, bool bottomUp,
          char* tmpSpace, F f)
    : A(A), bucketOffsets(bucketOffsets), n(n), m(m), bottomUp(bottomUp),
    tmpSpace(tmpSpace), f(f) { }
    
    heartbeat_private_activation_record_begin(heartbeat::edsl, iSort, 2)
      int s; int e;
    heartbeat_private_activation_record_end(heartbeat::edsl, iSort, sar, par, dc, get_dc)
    
    static
    dc get_dc() {
      return dc::stmts({
        dc::stmt([] (sar& s, par& p) {
          s.bits = pbbs::utils::log2Up(s.m);
          s.numBK = 1+s.n/(BUCKETS*8);
          // the temporary space is broken into 3 parts: B, Tmp and BK
          s.B = (E*) s.tmpSpace;
          s.Bsize =sizeof(E)*s.n;
          s.Tmp = (bIndexT*) (s.tmpSpace+s.Bsize); // one byte per item
          s.tmpSize = sizeof(bIndexT)*s.n;
          s.BK = (bucketsT*) (s.tmpSpace+s.Bsize+s.tmpSize);
        }),
        dc::cond({
          std::make_pair([] (sar& s, par& p) {
            return s.bits <= MAX_RADIX;
          }, dc::stmts({
            dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
              auto f = eBits<E,F>(s.bits,0,s.f);
              return heartbeat_call<radixStep<E,typeof(f),intT>>(st, pt, s.A, s.B, s.Tmp, s.BK, s.numBK, s.n, (intT) 1 << s.bits, true, f);
            }),
            dc::mk_if([] (sar& s, par& p) {
              return s.bucketOffsets != NULL;
            }, dc::stmts({
              dc::stmt([] (sar& s, par& p) {
                p.s = 0;
                p.e = s.m;
              }),
              dc::spawn_join([] (sar& s, par&, plt pt, stt st) {
                return heartbeat_call<sequence::copy<intT*, intT*>>(st, pt, s.BK[0], s.BK[0] + s.m, s.bucketOffsets);
                })
            })),
            dc::exit_function()
          })),
          std::make_pair([] (sar& s, par& p) {
            return s.bottomUp;
          },
            dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
            return heartbeat_call<radixLoopBottomUp<E, F, intT>>(st, pt, s.A, s.B, s.Tmp, s.BK, s.numBK, s.n, s.bits, true, s.f);
          }))
        },
          dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
            return heartbeat_call<radixLoopTopDown<E, F, intT>>(st, pt, s.A, s.B, s.Tmp, s.BK, s.numBK, s.n, s.bits, s.f);
        })),
        dc::mk_if([] (sar& s, par& p) {
          return s.bucketOffsets != NULL;
        }, dc::stmts({
          dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
            return sequence::fill3(st, pt, s.bucketOffsets, s.bucketOffsets + s.m, &s.n);
          }),
          dc::parallel_for_loop([] (sar& s, par& p) { p.s = 0; p.e = s.n-1; },
                                [] (par& p) { return std::make_pair(&p.s, &p.e); },
                                [] (sar& s, par& p, int lo, int hi) {
            auto f = s.f;
            auto A = s.A;
            auto bucketOffsets = s.bucketOffsets;
            for (auto i = lo; i != hi; i++) {
              intT v = f(A[i]);
              intT vn = f(A[i+1]);
              if (v != vn) {
                bucketOffsets[vn] = i+1;
              }
            }
          }),
          dc::stmt([] (sar& s, par& p) {
            s.bucketOffsets[s.f(s.A[0])] = 0;
          }),
          dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
            return sequence::scanIBack(st, pt, s.bucketOffsets, s.bucketOffsets, s.m, pbbs::utils::minF<intT>(), s.n, &s.tmp);
            })
        }))
      });
    }
    
  };
  
  template <class E, class F, class intT>
  typename iSort<E,F,intT>::cfg_type iSort<E,F,intT>::cfg = iSort<E,F,intT>::get_cfg();


  template <class E, class F, class intT>
  class iSort6 : public heartbeat::edsl::pcfg::shared_activation_record {
  public:
    
    E *A; intT* bucketOffsets; intT n; intT m; bool bottomUp;
    F f; typedef intT bucketsT[BUCKETS];

    char* ss;
    
    iSort6(E *A, intT* bucketOffsets, intT n, intT m, bool bottomUp, F f)
    : A(A), bucketOffsets(bucketOffsets), n(n), m(m), bottomUp(bottomUp),
    f(f) { }

    heartbeat_dc_declare(heartbeat::edsl, iSort6, sar, par, dc, get_dc)
    
    static
    dc get_dc() {
      return dc::stmts({
       	dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
          long x = iSortSpace<E,intT>(s.n);
          s.ss = (char*) malloc(x);
          return heartbeat_call<iSort<E,F,intT>>(st, pt, s.A, s.bucketOffsets, s.n, s.m, s.bottomUp, s.ss, s.f);
        }),
        dc::stmt([] (sar& s, par& p) {
          free(s.ss);
        })
      });
    }
    
  };
  
  template <class E, class F, class intT>
  typename iSort6<E,F,intT>::cfg_type iSort6<E,F,intT>::cfg = iSort6<E,F,intT>::get_cfg();

  template <class E, class F, class intT>
  stack_type iSort5(stack_type s, plt_type p, E *A, intT* bucketOffsets, intT n, intT m, F f) {
    return sequence::heartbeat_call<iSort6<E,F,intT>>(s, p, A, bucketOffsets, n, m, false, f);
  }

  template <class E, class F, class intT>
  stack_type iSort4(stack_type s, plt_type p, E *A, intT n, intT m, F f) {
    return sequence::heartbeat_call<iSort6<E,F,intT>>(s, p, A, (intT*)NULL, n, m, false, f);
  }

} // end namespace

template <class intT>
class integerSort : public heartbeat::edsl::pcfg::shared_activation_record {
public:

  intT* A; intT n;
  intT maxV;

  integerSort(intT* A, intT n) : A(A), n(n) { }

  heartbeat_dc_declare(heartbeat::edsl, integerSort, sar, par, dc, get_dc)

  static
  dc get_dc() {
    return dc::stmts({
      dc::spawn_join([] (sar& s, par&, plt p, stt st) {
        return sequence::reduce4(st, p, s.A, s.n, pbbs::utils::maxF<uintT>(), &(s.maxV));
      }),
      dc::spawn_join([] (sar& s, par&, plt p, stt st) {
        auto f = pbbs::utils::identityF<uintT>();
        return intSort::iSort5(st, p, s.A, (intT*)NULL, s.n, s.maxV+1, f);
      })
    });
  }

};

template <class intT>
typename integerSort<intT>::cfg_type integerSort<intT>::cfg = integerSort<intT>::get_cfg();


template <class T, class intT>
class integerSortPair : public heartbeat::edsl::pcfg::shared_activation_record {
public:

  std::pair<intT,T> *A; intT n;
  intT maxV;

  integerSortPair(std::pair<intT,T> *A, intT n) : A(A), n(n) { }

  heartbeat_dc_declare(heartbeat::edsl, integerSortPair, sar, par, dc, get_dc)

  static
  dc get_dc() {
    return dc::stmts({
      dc::spawn_join([] (sar& s, par&, plt p, stt st) {
        auto f = pbbs::utils::maxF<uintT>();
        auto g = pbbs::utils::firstF<uintT,T>();
        return sequence::mapReduce(st, p, s.A, s.n, f, g, &(s.maxV));
      }),
      dc::spawn_join([] (sar& s, par&, plt p, stt st) {
        return intSort::iSort5(st, p, s.A, (intT*)NULL, s.n, s.maxV+1, pbbs::utils::firstF<uintT,T>());
      })
    });
  }

};

template <class T, class intT>
typename integerSortPair<T, intT>::cfg_type integerSortPair<T,intT>::cfg = integerSortPair<T,intT>::get_cfg();

} // end namespace

#endif
