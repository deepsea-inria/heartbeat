#ifdef USE_CILK_PLUS
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#endif

#include "heartbeat.hpp"

#ifndef _HEARTBEAT_BENCH_H_
#define _HEARTBEAT_BENCH_H_

void cilk_set_nb_cores(int proc) {
#ifdef USE_CILK_PLUS
  __cilkrts_set_param("nworkers", std::to_string(proc).c_str());
#endif
}

namespace heartbeatbench {

namespace cmdline = deepsea::cmdline;
  
void initialize(int argc, char** argv) {
  heartbeat::initialize_runtime(argc, argv);
  cilk_set_nb_cores(cmdline::parse_or_default("proc", 1));
}

void trigger_cilk() {
  printf("");
}

template <class Function>
void run_and_report_elapsed_time(const Function& f) {
#ifdef CILK_RUNTIME_WITH_STATS
  cilk_spawn trigger_cilk();
  cilk_sync;
  __cilkg_take_snapshot_for_stats();
#endif
  auto start = std::chrono::system_clock::now();
  f();
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<float> diff = end - start;
  printf ("exectime %.3lf\n", diff.count());
  /* recall: if using the custom cilk runtime, need to set the
   * environment variable as such:
   *   export LD_LIBRARY_PATH=/home/rainey/cilk-plus-rts/lib:$LD_LIBRARY_PATH
   */
#ifdef CILK_RUNTIME_WITH_STATS
  __cilkg_dump_heartbeat_stats_to_stderr();
#endif
}

} // end namespace

#endif /*! _HEARTBEAT_BENCH_H_ */
