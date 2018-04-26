#include <iostream>
#include <chrono>

#ifdef USE_CILK_PLUS
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#else
#define cilk_for for
#endif

#include "heartbeatbench.hpp"

namespace sched = heartbeat::sched;
namespace cmdline = deepsea::cmdline;
namespace dsl = heartbeat::edsl;

int branching_factor = 1024;

using value_type = int;

template <class T>
T* malloc_array(size_t n) {
  return (T*)malloc(n * sizeof(T));
}

value_type reduce_serial(int lo, int hi, value_type z, value_type* src) {
  value_type r = z;
  for (int i = lo; i < hi; i++) {
    r += src[i];
  }
  return r;
}

value_type scan_serial(int lo, int hi, value_type z, value_type* src, value_type* dst) {
  value_type r = z;
  for (int i = lo; i < hi; i++) {
    dst[i] = r;
    r += src[i];
  }
  return r;
}

static inline
int get_nb_blocks(int k, int n) {
  return 1 + ((n - 1) / k);
}

class range_type {
public:
  int lo;
  int hi;
  range_type(int lo, int hi)
  : lo(lo), hi(hi) { }
};

static inline
range_type get_rng(int k, int n, int i) {
  int lo = i * k;
  int hi = std::min(lo + k, n);
  return range_type(lo, hi);
}

static inline
range_type get_rng(int k, int n, int lo, int hi) {
  assert(lo < hi);
  int lo2 = lo * k;
  int hi2 = std::min(lo2 + k, std::min(hi * k, n));
  return range_type(lo2, hi2);
}

class scan_dc : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  int n;
  value_type z;
  value_type* src;
  value_type* dst;
  
  scan_dc() { }
  
  scan_dc(int n, value_type z, value_type* src, value_type* dst)
  : n(n), z(z), src(src), dst(dst) { }
  
  value_type* partials;
  value_type* scans;
  int m;
  
  heartbeat_private_activation_record_begin(heartbeat::edsl, scan_dc, 2)
    int lo; int hi;
  heartbeat_private_activation_record_end(heartbeat::edsl, scan_dc, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return
    dc::mk_if([] (sar& s, par& p) { return s.n <= branching_factor; }, dc::stmt([] (sar& s, par& p) {
      scan_serial(0, s.n, s.z, s.src, s.dst);
    }), dc::stmts({
      dc::stmt([] (sar& s, par& p) {
        s.m = get_nb_blocks(branching_factor, s.n);
        s.partials = malloc_array<value_type>(s.m);
        p.lo = 0;
        p.hi = s.m;
      }),
      dc::parallel_for_loop([] (sar&, par& p) { return p.lo < p.hi; },
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            dc::stmt([] (sar& s, par& p) {
        auto rng = get_rng(branching_factor, s.n, p.lo, p.hi);
        s.partials[p.lo] = reduce_serial(rng.lo, rng.hi, 0, s.src);
        p.lo++;
      })),
      dc::stmt([] (sar& s, par& p) {
        s.scans = malloc_array<value_type>(s.m);
      }),
      dc::spawn_join([] (sar& s, par& p, plt pt, stt st) {
        return heartbeat_call<scan_dc>(st, pt, s.m, s.z, s.partials, s.scans);
      }),
      dc::stmt([] (sar& s, par& p) {
        free(s.partials);
        p.lo = 0;
        p.hi = s.m;
      }),
      dc::parallel_for_loop([] (sar&, par& p) { return p.lo < p.hi; },
                            [] (par& p) { return std::make_pair(&p.lo, &p.hi); },
                            dc::stmt([] (sar& s, par& p) {
        auto rng = get_rng(branching_factor, s.n, p.lo, p.hi);
        scan_serial(rng.lo, rng.hi, s.scans[p.lo], s.src, s.dst);
        p.lo++;
      })),
      dc::stmt([] (sar& s, par& p) {
        free(s.scans);
      })
    }));
  }
  
};

heartbeat_pcfg_allocate(scan_dc, get_cfg)

value_type scan_cilk(int n, value_type z, value_type* src, value_type* dst) {
  int k = branching_factor;
  value_type r = z;
  if (n <= k) {
    r = scan_serial(0, n, r, src, dst);
  } else {
    int m = get_nb_blocks(k, n);
    value_type* partials = malloc_array<int>(m);
    cilk_for (int i = 0; i < m; i++) {
      auto rng = get_rng(k, n, i);
      partials[i] = reduce_serial(rng.lo, rng.hi, 0, src);
    }
    value_type* scans = malloc_array<int>(m);
    r = scan_cilk(m, z, partials, scans);
    free(partials);
    cilk_for (int i = 0; i < m; i++) {
      auto rng = get_rng(k, n, i);
      scan_serial(rng.lo, rng.hi, scans[i], src, dst);
    }
    free(scans);
  }
  return r;
}

int main(int argc, char** argv) {
  heartbeatbench::initialize(argc, argv);
  int n = cmdline::parse<int>("n");
  branching_factor = cmdline::parse_or_default("branching_factor", branching_factor);
  value_type* src = malloc_array<value_type>(n);
  value_type* dst = malloc_array<value_type>(n);
  for (int i = 0; i < n; i++) {
    src[i] = 1;
  }
  cmdline::dispatcher d;
  d.add("serial", [&] {
    scan_serial(0, n, 0, src, dst);
  });
  d.add("cilk", [&] {
    scan_cilk(n, 0, src, dst);
  });
  d.add("heartbeat", [&] {
    heartbeat::launch_interpreter<scan_dc>(n, 0, src, dst);
  });
  heartbeatbench::run_and_report_elapsed_time([&] {
    d.dispatch_or_default("algorithm", "serial");
  });
#ifndef NDEBUG
  value_type* src2 = malloc_array<value_type>(n);
  value_type* dst2 = malloc_array<value_type>(n);
  std::copy(src, src + n, src2);
  scan_serial(0, n, 0, src2, dst2);
  for (int i = 0; i < n; i++) {
    auto test = dst[i];
    auto ref = dst2[i];
    assert(test == ref);
  }
  free(src2);
  free(dst2);
#endif
  free(src);
  free(dst);
  return 0;
}
