
#include <iostream>
#include <chrono>

#include "heartbeatbench.hpp"

int fib(int n) {
  if (n <= 1) {
    return n;
  }
  return fib(n - 1) + fib(n - 2);
}

int cutoff = 2;

#ifdef USE_CILK_PLUS
int fib_cilk(int n) {
  if (n <= cutoff) {
    return fib(n);
  }
  int d1 = cilk_spawn fib_cilk(n - 1);
  int d2 = fib_cilk(n - 2);
  cilk_sync;
  return d1 + d2;
}
#endif

namespace sched = heartbeat::sched;

class fib_manual : public sched::vertex {
public:
  
  enum { entry, join, exit };
  int trampoline = entry;
  
  void yield_with(int target) {
    trampoline = target;
  }
  
  int n; int* dp;
  int d1, d2;
  
  fib_manual(int n, int* dp)
  : vertex(), n(n), dp(dp) { }
  
  int nb_strands() {
    if (trampoline == exit) {
      return 0;
    }
    return 1;
  }
  
  heartbeat::fuel::check_type run() {
    switch (trampoline) {
      case entry: {
        if (n <= cutoff) {
          *dp = fib(n);
          yield_with(exit);
          break;
        }
        fib_manual* b1 = new fib_manual(n - 1, &d1);
        fib_manual* b2 = new fib_manual(n - 2, &d2);
        yield_with(join);
        sched::new_edge(b2, this);
        sched::new_edge(b1, this);
        sched::release(b2);
        sched::release(b1);
        break;
      }
      case join: {
        *dp = d1 + d2;
        yield_with(exit);
        break;
      }
      default:
        assert(false);
    }
    return heartbeat::fuel::check_no_promote;
  }
  
  sched::vertex_split_type split(int nb) {
    assert(false); // impossible
    return sched::make_vertex_split(nullptr, nullptr);
  }
  
};

class fib_cfg : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  int n; int* dp;
  int d1; int d2;
  
  fib_cfg() { }
  
  fib_cfg(int n, int* dp)
  : n(n), dp(dp) { }
  
  heartbeat_pcfg_default_private_activation_record(heartbeat::edsl::pcfg)
  heartbeat_pcfg_declare(heartbeat::edsl, fib_cfg, sar, par, bb)
  
  static
  cfg_type get_cfg() {
    enum { entry, branch1, branch2, combine, nb_blocks };
    auto lab = branch1;
    cfg_type cfg(nb_blocks);
    cfg[entry] = bb::conditional_jump([=] (sar& s, par&) {
      if (s.n <= cutoff) {
        *s.dp = fib(s.n);
        return heartbeat::edsl::pcfg::exit_block_label; // going to exit block
      }
      return (int)lab;  // going to branch1
    });
    cfg[branch1] = bb::spawn2_join([] (sar& s, par&, plt p, stt st) {
      return heartbeat_call<fib_cfg>(st, p, s.n - 1, &s.d1);
    }, branch2);
    cfg[branch2] = bb::spawn_join([] (sar& s, par&, plt p, stt st) {
      return heartbeat_call<fib_cfg>(st, p, s.n - 2, &s.d2);
    }, combine);
    cfg[combine] = bb::unconditional_jump([] (sar& s, par&) {
      *s.dp = s.d1 + s.d2;
    }, heartbeat::edsl::pcfg::exit_block_label);
    return cfg;
  }
  
};

heartbeat_pcfg_allocate(fib_cfg, get_cfg)

class fib_dc : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  int n; int* dp;
  int d1; int d2;
  
  fib_dc() { }
  
  fib_dc(int n, int* dp)
  : n(n), dp(dp) { }

  heartbeat_dc_declare(heartbeat::edsl, fib_dc, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return
    dc::mk_if([] (sar& s, par&) { return s.n <= cutoff; },
      dc::stmt([] (sar& s, par&) { *s.dp = fib(s.n); }),
      dc::stmts({
        dc::spawn2_join(
            [] (sar& s, par&, plt p, stt st) {
              return heartbeat_call<fib_dc>(st, p, s.n - 1, &s.d1); },
            [] (sar& s, par&, plt p, stt st) {
              return heartbeat_call<fib_dc>(st, p, s.n - 2, &s.d2); }),
         dc::stmt([] (sar& s, par&) { *s.dp = s.d1 + s.d2; }),
      })
    );
  }

};

heartbeat_pcfg_allocate(fib_dc, get_cfg)

namespace cmdline = deepsea::cmdline;

int main(int argc, char** argv) {
  heartbeatbench::initialize(argc, argv);
  int n = cmdline::parse<int>("n");
  cutoff = cmdline::parse_or_default("cutoff", cutoff);
  int result = -1;
  cmdline::dispatcher d;
  d.add("sequential", [&] {
    result = fib(n);
  });
  d.add("dag", [&] {
    heartbeat::launch(new fib_manual(n, &result));
  });
#ifdef USE_CILK_PLUS
  d.add("cilk", [&] {
    result = fib_cilk(n);
  });
#endif
  d.add("pcfg", [&] {
    heartbeat::launch_interpreter<fib_cfg>(n, &result);
  });
  d.add("dc", [&] {
    heartbeat::launch_interpreter<fib_dc>(n, &result);
  });
  heartbeatbench::run_and_report_elapsed_time([&] {
    d.dispatch("algorithm");
  });
  auto fn = fib(n);
  assert(result == fn);
  return 0;
}
