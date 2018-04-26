
#include "machine.hpp"
#include "scheduler.hpp"
#include "edsl.hpp"
#include "cmdline.hpp"
#include "grain.hpp"
#include "fuel.hpp"

#ifndef _HEARTBEAT_H_
#define _HEARTBEAT_H_

namespace heartbeat {

namespace cmdline = deepsea::cmdline;
  
template <class Function>
class call_and_report_elapsed : public edsl::pcfg::shared_activation_record {
public:
  
  Function f;
  std::chrono::time_point<std::chrono::system_clock> start;
  
  call_and_report_elapsed() { }
  
  call_and_report_elapsed(const Function& f)
  : f(f) { }
  
  heartbeat_dc_declare(heartbeat::edsl, call_and_report_elapsed, sar, par, dc, get_dc)
  
  static
  dc get_dc() {
    return dc::stmts({
      dc::stmt([] (sar& s, par&) {
        s.start = std::chrono::system_clock::now();
      }),
      dc::profile_statement([] (sar& s, par& p) {
          auto work = nullptr;
          auto span = nullptr;
          return std::make_pair(work, span);
          
        },
        dc::spawn_join([] (sar& s, par&, plt, stt st) {
          return s.f(st);
        }),
        [] (sar& s, par& p, uint64_t work, uint64_t span) {
          /*
          double ticks_per_second = cpu_frequency_ghz * 1000000000.0;
          double work_sec = ((double)work) / ticks_per_second;
          double span_sec = ((double)span) / ticks_per_second;
          printf("work %.5lf\n", work_sec);
          printf("span %.5lf\n", span_sec); */
      }),
      dc::stmt([] (sar& s, par&) {
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<float> diff = end - s.start;
        printf ("exectime %.3lf\n", diff.count());
        sched::should_exit = true;
      })
    });
  }
  
};

template <class Function>
typename call_and_report_elapsed<Function>::cfg_type call_and_report_elapsed<Function>::cfg = call_and_report_elapsed<Function>::get_cfg();
  
void initialize_runtime(int argc, char** argv) {
  cmdline::set(argc, argv);
  atomic::init_print_lock();
  machine::initialize_hwloc();
  machine::initialize_cpuinfo();
  auto scheduler = cmdline::parse_or_default_string("scheduler", "work_stealing");
  if (scheduler == "steal_half_work_stealing") {
    sched::scheduler = sched::steal_half_work_stealing_tag;
  } else if (scheduler == "heartbeat_work_stealing") {
    sched::scheduler = sched::heartbeat_work_stealing_tag;
  } else if (scheduler == "steal_one_work_stealing") {
    sched::scheduler = sched::steal_one_work_stealing_tag;
  } else if (scheduler == "work_stealing") {
    sched::scheduler = sched::work_stealing_tag;
  } else if (scheduler == "concurrent_deques_work_stealing") {
    sched::scheduler = sched::concurrent_deques_work_stealing_tag;
  } else {
    atomic::die("bogus scheduler\n");
  }
  edsl::pcfg::never_promote = cmdline::parse_or_default_bool("never_promote", edsl::pcfg::never_promote);
  double promotion_threshold_usec = 30.0;
  if (edsl::pcfg::never_promote) {
    promotion_threshold_usec = 10000000.0;
  }
  promotion_threshold_usec = cmdline::parse_or_default_double("promotion_threshold", promotion_threshold_usec);
  printf("promotion_threshold_usec %.3f\n", promotion_threshold_usec);
  fuel::initialize(machine::cpu_frequency_ghz, promotion_threshold_usec * 1000.0);
  double grain_usec = promotion_threshold_usec / 4.0;
  grain_usec = cmdline::parse_or_default_double("grain", grain_usec);
  grain::initialize(machine::cpu_frequency_ghz, grain_usec * 1000.0, promotion_threshold_usec * 1000.0);
  std::cout << "cpu_frequency_ghz " << machine::cpu_frequency_ghz << std::endl;
}
  
template <class Init>
void launch(int nb_workers, const Init& init) {
  logging::log_buffer::initialize();
  stats::initialize();
  sched::vertex* v = init();
  stats::on_enter_launch();
  sched::launch_scheduler(nb_workers, v);
  stats::on_exit_launch();
  stats::report();
  logging::log_buffer::output();
  data::perworker::reset();
}

void launch(sched::vertex* v, int nb_workers) {
  launch(nb_workers, [=] { return v; });
}
  
void launch(sched::vertex* v) {
  launch(v, cmdline::parse_or_default("proc", 1));
}

template <class F>
void launch_interpreter_via_lambda(const F& f) {
  /* thanks to buggy GCC, the code here will crash the compiler...
  launch(cmdline::parse_or_default("proc", 1), [&] {
    auto interp = new edsl::pcfg::interpreter<edsl::pcfg::stack_type>;
    auto f = [=] (edsl::pcfg::stack_type st) {
      return edsl::pcfg::push_call<Shared_activation_record>(st, args...);
    };
    using t = call_and_report_elapsed<typeof(f)>;
    interp->stack = edsl::pcfg::push_call<t>(interp->stack, f);
    return interp;
  });
  */
  int nb_workers = cmdline::parse_or_default("proc", 1);
  logging::log_buffer::initialize();
  stats::initialize();
  auto interp = new edsl::pcfg::interpreter;
  using t = call_and_report_elapsed<F>;
  interp->stack = edsl::pcfg::push_call<t>(interp->stack,
                                           edsl::pcfg::cactus::Parent_link_sync,
                                           f);
  stats::on_enter_launch();
  sched::launch_scheduler(nb_workers, interp);
  stats::on_exit_launch();
  stats::report();
  logging::log_buffer::output();
  data::perworker::reset();
}

template <class Shared_activation_record, class ...Args>
void launch_interpreter(Args... args) {
  using sar = Shared_activation_record;
  launch_interpreter_via_lambda([=] (edsl::pcfg::stack_type st) {
    return edsl::pcfg::push_call<sar>(st,
                                      edsl::pcfg::cactus::Parent_link_sync,
                                      args...);
  });
}
  
} // end namespace

#endif /*! _HEARTBEAT_H_ */
