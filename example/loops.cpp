#include <iostream>
#include <chrono>

#include "heartbeatbench.hpp"

namespace sched = heartbeat::sched;
namespace cmdline = deepsea::cmdline;
namespace dsl = heartbeat::edsl;

int cutoff = 1;

class sequential_loop_0 : public heartbeat::edsl::pcfg::shared_activation_record {
public:
    
  int n;
  int* a; int lo; int hi;
  
  sequential_loop_0() { }
  
  sequential_loop_0(int n)
  : n(n) { }
  
  heartbeat_dc_declare(heartbeat::edsl, sequential_loop_0, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par&) {
        s.lo = 0;
        s.hi = s.n;
        s.a = new int[s.n];
      }),
      dc::sequential_loop([] (sar& s, par&) { return s.lo != s.hi; }, dc::stmt([] (sar& s, par&) {
        int lo = s.lo;
        int hi = std::min(s.hi, lo + cutoff);
        int* a = s.a;
        for (; lo != hi; lo++) {
          a[lo] = 0xdeadbeef;
        }
        s.lo = lo;
      })),
      dc::stmt([] (sar& s, par&) {
        for (int i = 0; i < s.n; i++) {
          assert(s.a[i] == 0xdeadbeef);
        }
        delete [] s.a;
      })
    });
  }
  
};

heartbeat_pcfg_allocate(sequential_loop_0, get_cfg)

class sequential_loop_1 : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  int n;
  int* a; int lo1; int hi1; int lo2; int hi2;
  
  sequential_loop_1() { }
  
  sequential_loop_1(int n)
  : n(n) { }
  
  heartbeat_dc_declare(heartbeat::edsl, sequential_loop_1, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par&) {
        s.lo1 = 0;
        s.hi1 = s.n;
        s.a = new int[s.n * s.n];
      }),
      dc::sequential_loop([] (sar& s, par&) { return s.lo1 != s.hi1; },
        dc::stmts({ dc::stmt([] (sar& s, par&) {
          s.lo2 = 0;
          s.hi2 = s.n;
        }),
        dc::sequential_loop([] (sar& s, par&) { return s.lo2 != s.hi2; }, dc::stmt([] (sar& s, par&) {
          int i = s.lo1;
          int j = s.lo2;
          s.a[i * s.n + j] = 0xdeadbeef;
          s.lo2++;
        })),
        dc::stmt([] (sar& s, par&) {
          s.lo1++;
        })
      })),
      dc::stmt([] (sar& s, par&) {
        for (int i = 0; i < s.n * s.n; i++) {
          assert(s.a[i] == 0xdeadbeef);
        }
        delete [] s.a;
      })
    });
  }
  
};

heartbeat_pcfg_allocate(sequential_loop_1, get_cfg)

class sequential_loop_2 : public heartbeat::edsl::pcfg::shared_activation_record {
public:
    
  int n;
  int* a; int lo; int hi;
  
  sequential_loop_2() { }
  
  sequential_loop_2(int n)
  : n(n) { }
  
  heartbeat_dc_declare(heartbeat::edsl, sequential_loop_2, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par&) {
        s.a = new int[s.n];
      }),
      dc::sequential_loop([] (sar& s, par&) {
        s.lo = 0;
        s.hi = s.n;
      }, [] (sar& s, par&) {
        return std::make_pair(&s.lo, &s.hi);
      }, [] (sar& s, par&, int lo, int hi) {
        auto a = s.a;          
        for (int i = lo; i < hi; i++) {
          a[i] = 0xdeadbeef;
        }
      }),
      dc::stmt([] (sar& s, par&) {
        for (int i = 0; i < s.n; i++) {
          assert(s.a[i] == 0xdeadbeef);
        }
        delete [] s.a;
      })
    });
  }
  
};

heartbeat_pcfg_allocate(sequential_loop_2, get_cfg)

class parallel_loop_0 : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  int n;
  int* a;
  
  parallel_loop_0() { }
  
  parallel_loop_0(int n)
  : n(n) { }
  
  heartbeat_private_activation_record_begin(heartbeat::edsl, parallel_loop_0, 1)
    int lo; int hi;
  heartbeat_private_activation_record_end(heartbeat::edsl, parallel_loop_0, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par& p) {
        p.lo = 0;
        p.hi = s.n;
        s.a = new int[s.n];
      }),
      dc::parallel_for_loop([] (sar&, par& p) { return p.lo != p.hi; },
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            dc::stmt([] (sar& s, par& p) {
        s.a[p.lo] = 0xdeadbeef;
        p.lo++;
      })),
      dc::stmt([] (sar& s, par&) {
        for (int i = 0; i < s.n; i++) {
          assert(s.a[i] == 0xdeadbeef);
        }
        delete [] s.a;
      })
    });
  }
  
};

heartbeat_pcfg_allocate(parallel_loop_0, get_cfg)

class parallel_loop_1 : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  int n;
  int* a;
  
  parallel_loop_1() { }
  
  parallel_loop_1(int n)
  : n(n) { }
  
  heartbeat_private_activation_record_begin(heartbeat::edsl, parallel_loop_1, 2)
    int lo1; int hi1;
    int lo2; int hi2;
  heartbeat_private_activation_record_end(heartbeat::edsl, parallel_loop_1, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par& p) {
        p.lo1 = 0;
        p.hi1 = s.n;
        s.a = new int[s.n * s.n];
      }),
      dc::parallel_for_loop([] (sar&, par& p) { return p.lo1 != p.hi1; },
                            [] (par& p) { return std::make_pair(&p.lo1, &p.hi1); },
                            dc::stmts({
        dc::parallel_for_loop([] (sar& s, par& p) { p.lo2 = 0; p.hi2 = s.n; },
                              [] (par& p) { return std::make_pair(&p.lo2, &p.hi2); },
                              [] (sar& s, par& p, int lo2, int hi2) {
          auto a = s.a;
          auto n = s.n;
          auto i = p.lo1;
          for (auto j = lo2; j != hi2; j++) {
            a[i * n + j] = 0xdeadbeef;
          }
        }),
        dc::stmt([] (sar& s, par& p) {
          p.lo1++;
        })
      })),
      dc::stmt([] (sar& s, par&) {
        for (int i = 0; i < s.n * s.n; i++) {
          assert(s.a[i] == 0xdeadbeef);
        }
        delete [] s.a;
      })
    });
  }
  
};

heartbeat_pcfg_allocate(parallel_loop_1, get_cfg)

class parallel_loop_2 : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  int n;
  int* a;
  
  parallel_loop_2() { }
  
  parallel_loop_2(int n)
  : n(n) { }
  
  heartbeat_private_activation_record_begin(heartbeat::edsl, parallel_loop_2, 1)
    int lo; int hi;
  heartbeat_private_activation_record_end(heartbeat::edsl, parallel_loop_2, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::parallel_for_loop([] (sar& s, par& p) {
        p.lo = 0;
        p.hi = s.n;
        s.a = new int[s.n];
      }, [] (par& p) {
        return std::make_pair(&p.lo, &p.hi);
      }, [] (sar& s, par& p, int lo, int hi) {
        auto a = s.a;
        for (auto i = lo; i != hi; i++) {
          a[i] = 0xdeadbeef;
        }
      }, __LINE__, __FILE__),                       
      dc::stmt([] (sar& s, par&) {
        for (int i = 0; i < s.n; i++) {
          assert(s.a[i] == 0xdeadbeef);
        }
        delete [] s.a;
      })
    });
  }
  
};

heartbeat_pcfg_allocate(parallel_loop_2, get_cfg)

class parallel_combine_0 : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  static constexpr int nb_loops = 1;
  
  int n;
  int* a;
  int* b;
  int result;
  
  parallel_combine_0() { }
  
  parallel_combine_0(int n)
  : n(n) { }
  
  class private_activation_record
  : public dsl::pcfg::parallel_loop_private_activation_record<parallel_combine_0,
  private_activation_record> {
  public:
    
    private_activation_record() {
      private_activation_record::initialize_descriptors();
    }
    
    dsl::pcfg::parallel_combine_activation_record _ar;
    dsl::pcfg::parallel_loop_activation_record* _heartbeat_loop_activation_record_of(dsl::pcfg::parallel_loop_id_type id) {
      assert(id == 0);
      return &_ar;
    }
    
    int lo; int hi;
    int acc = 0;
    
  };
  
  heartbeat_dc_loop_declare(heartbeat::edsl, parallel_combine_0, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par& p) {
        p.lo = 0;
        p.hi = s.n;
        p.acc = 0;
        s.a = new int[s.n];
        s.b = new int[s.n];
        for (int i = 0; i < s.n; i++) {
          s.a[i] = i % 1024;
          s.b[i] = 0;
        }
      }),
      dc::parallel_combine_loop([] (sar&, par& p) { return p.lo != p.hi; },
                                [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                                [] (sar&, par& p) { p.acc = 0; },
                                [] (sar&, par& p, par& destination) { destination.acc += p.acc; },
                                dc::stmt([] (sar& s, par& p) {
        p.acc += s.a[p.lo];
        s.b[p.lo] = 0xdeadbeef;
        p.lo++;
      })),
      dc::stmt([] (sar& s, par& p) {
        s.result = p.acc;
#ifndef NDEBUG
        int acc2 = 0;
        for (int i = 0; i < s.n; i++) {
          acc2 += s.a[i];
          assert(s.b[i] == 0xdeadbeef);
        }
        assert(s.result == acc2);
#endif
        delete [] s.a;
        delete [] s.b;
      })
    });
  }
  
};

heartbeat_pcfg_allocate(parallel_combine_0, get_cfg)

class parallel_combine_1 : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  static constexpr int nb_loops = 1;
  
  int n;
  int* a;
  int* b;
  int result;
  
  parallel_combine_1() { }
  
  parallel_combine_1(int n)
  : n(n) { }
  
  class private_activation_record
  : public dsl::pcfg::parallel_loop_private_activation_record<parallel_combine_1,
  private_activation_record> {
  public:
    
    private_activation_record() {
      private_activation_record::initialize_descriptors();
    }
    
    dsl::pcfg::parallel_combine_activation_record _ar;
    dsl::pcfg::parallel_loop_activation_record* _heartbeat_loop_activation_record_of(dsl::pcfg::parallel_loop_id_type id) {
      assert(id == 0);
      return &_ar;
    }
    
    int lo; int hi;
    int acc = 0;
    
  };
  
  heartbeat_dc_loop_declare(heartbeat::edsl, parallel_combine_1, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par& p) {
        s.a = new int[s.n];
        s.b = new int[s.n];
        for (int i = 0; i < s.n; i++) {
          s.a[i] = i % 1024;
          s.b[i] = 0;
        }
      }),
      dc::parallel_combine_loop([] (sar& s, par& p) { p.lo = 0; p.hi = s.n; p.acc = 0; },
                                [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                                [] (sar&, par& p) { p.acc = 0; },
                                [] (sar&, par& p, par& destination) { destination.acc += p.acc; },
                                [] (sar& s, par& p, int lo, int hi) {
        auto a = s.a;
        auto b = s.b;
        auto acc = p.acc;
        for (auto i = lo; i != hi; i++) {
          acc += a[i];
          b[i] = 0xdeadbeef;
        }
        p.acc = acc;
      }),
      dc::stmt([] (sar& s, par& p) {
        s.result = p.acc;
#ifndef NDEBUG
        int acc2 = 0;
        for (int i = 0; i < s.n; i++) {
          acc2 += s.a[i];
          assert(s.b[i] == 0xdeadbeef);
        }
        assert(s.result == acc2);
#endif
        delete [] s.a;
        delete [] s.b;
      })
    });
  }
  
};

heartbeat_pcfg_allocate(parallel_combine_1, get_cfg)

int main(int argc, char** argv) {
  heartbeatbench::initialize(argc, argv);
  int n = cmdline::parse<int>("n");
  cutoff = cmdline::parse_or_default("cutoff", cutoff);
  heartbeatbench::run_and_report_elapsed_time([&] {
    cmdline::dispatcher d;
    d.add("sequential_loop_0", [=] {
      heartbeat::launch_interpreter<sequential_loop_0>(n);
    });
    d.add("sequential_loop_1", [=] {
      heartbeat::launch_interpreter<sequential_loop_1>(n);
    });
    d.add("sequential_loop_2", [=] {
      heartbeat::launch_interpreter<sequential_loop_2>(n);
    });
    d.add("parallel_loop_0", [=] {
      heartbeat::launch_interpreter<parallel_loop_0>(n);
    });
    d.add("parallel_loop_1", [=] {
      heartbeat::launch_interpreter<parallel_loop_1>(n);
    });
    d.add("parallel_loop_2", [=] {
      heartbeat::launch_interpreter<parallel_loop_2 >(n);
    });
    d.add("parallel_combine_0", [=] {
      heartbeat::launch_interpreter<parallel_combine_0>(n);
    });
    d.add("parallel_combine_1", [=] {
      heartbeat::launch_interpreter<parallel_combine_1>(n);
    });
    d.dispatch_or_default("function", "sequential_loop_0");
  });
  return 0;
}
