
#include <iostream>

#include "heartbeatbench.hpp"

namespace sched = heartbeat::sched;
namespace cmdline = deepsea::cmdline;
namespace dsl = heartbeat::edsl;

template <class T>
T* malloc_array(size_t n) {
  return (T*)malloc(n * sizeof(T));
}

int cutoff = 1;

class bintree_rec : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  int lo; int hi; int* a;
  sched::incounter** join;
  
  int mid;
  
  bintree_rec() { }
  
  bintree_rec(int lo, int hi, int* a, sched::incounter** join)
  : lo(lo), hi(hi), a(a), join(join) { }
  
  heartbeat_dc_declare(heartbeat::edsl, bintree_rec, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return
    dc::mk_if([] (sar& s, par&) {
      return s.hi - s.lo <= cutoff; },
      dc::stmt([] (sar& s, par&) {
        for (int i = s.lo; i < s.hi; i++) {
          s.a[s.lo]++;
        }
      }),
      dc::stmts({
        dc::stmt([] (sar& s, par&) { s.mid = (s.lo + s.hi) / 2; }),
        dc::spawn_minus([] (sar& s, par&, plt pt, stt st) {
          return heartbeat_call<bintree_rec>(st, pt, s.lo, s.mid, s.a, s.join); },
          [] (sar& s, par&) { return s.join; }),
        dc::spawn_minus([] (sar& s, par&, plt pt, stt st) {
          return heartbeat_call<bintree_rec>(st, pt, s.mid, s.hi, s.a, s.join); },
          [] (sar& s, par&) { return s.join; })
      }));
  }
  
};

heartbeat_pcfg_allocate(bintree_rec, get_cfg)

class bintree : public heartbeat::edsl::pcfg::shared_activation_record {
public:
  
  int n; int* a;
  sched::incounter* join = nullptr;
    
  bintree() { }
  
  bintree(int n, int* a)
  : n(n), a(a) { }
  
  heartbeat_dc_declare(heartbeat::edsl, bintree, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return
    dc::join_plus([] (sar& s, par&, plt pt, stt st) {
      return heartbeat_call<bintree_rec>(st, pt, 0, s.n, s.a, &s.join); },
      [] (sar& s, par&) { return &s.join; });
  }
  
};

heartbeat_pcfg_allocate(bintree, get_cfg)

int main(int argc, char** argv) {
  heartbeatbench::initialize(argc, argv);
  int n = cmdline::parse<int>("n");
  cutoff = cmdline::parse_or_default("cutoff", cutoff);
  int* a = malloc_array<int>(n);
  for (int i = 0; i < n; i++) {
    a[i] = 0xdeadbeef;
  }
  heartbeat::launch_interpreter<bintree>(n, a);
#ifndef NDEBUG
  for (int i = 0; i < n; i++) {
    assert(a[i] == 0xdeadbeef + 1);
  }
#endif
  free(a);
  return 0;
}
