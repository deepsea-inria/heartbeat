#include <chrono>
#include <map>
#include <iostream>

#include "perworker.hpp"

#ifndef _HEARTBEAT_STATS_H_
#define _HEARTBEAT_STATS_H_

namespace heartbeat {

template <bool enabled>
class stats_base {
public:
  
  using time_point_type = std::chrono::time_point<std::chrono::system_clock>;
  
private:
  
  using counter_id_type = enum {
    nb_promotions,
    nb_steals,
    nb_stacklet_allocations,
    nb_stacklet_deallocations,
    nb_counters
  };
  
  static
  const char* name_of_counter(counter_id_type id) {
    std::map<counter_id_type, const char*> names;
    names[nb_promotions] = "nb_promotions";
    names[nb_steals] = "nb_steals";
    names[nb_stacklet_allocations] = "nb_stacklet_allocations";
    names[nb_stacklet_deallocations] = "nb_stacklet_deallocations";
    return names[id];
  }

  using private_counters = struct {
    long counters[nb_counters];
  };
  
  static
  data::perworker::array<private_counters> all_counters;
  
  static inline
  void increment(counter_id_type id) {
    if (! enabled) {
      return;
    }
    all_counters.mine().counters[id]++;
  }
  
  static
  time_point_type enter_launch_time;
  
  static
  double launch_duration;
  
  static
  data::perworker::array<double> all_total_idle_time;
  
  static
  double since(time_point_type start) {
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    return elapsed.count();
  }
  
public:
  
  static inline
  void on_promotion() {
    increment(nb_promotions);
  }
  
  static inline
  void on_steal() {
    increment(nb_steals);
  }
  
  static inline
  void on_stacklet_allocation() {
    increment(nb_stacklet_allocations);
  }
  
  static inline
  void on_stacklet_deallocation() {
    increment(nb_stacklet_deallocations);
  }
  
  static
  void on_enter_launch() {
    enter_launch_time = std::chrono::system_clock::now();
  }
  
  static
  void on_exit_launch() {
    launch_duration = since(enter_launch_time);
  }
  
  static
  time_point_type on_enter_acquire() {
    if (! enabled) {
      return time_point_type();
    }
    return std::chrono::system_clock::now();
  }
  
  static
  void on_exit_acquire(time_point_type enter_acquire_time) {
    if (! enabled) {
      return;
    }
    all_total_idle_time.mine() += since(enter_acquire_time);
  }
  
  static
  void initialize() {
    for (int counter_id = 0; counter_id < nb_counters; counter_id++) {
      all_counters.for_each([&] (int, private_counters& b) {
        b.counters[counter_id] = 0;
      });
    }
    all_total_idle_time.for_each([&] (int, double& d) {
      d = 0.0;
    });
  }
  
  static
  void report() {
    if (! enabled) {
      return;
    }
    for (int counter_id = 0; counter_id < nb_counters; counter_id++) {
      long counter_value = 0;
      all_counters.for_each([&] (int, private_counters& b) {
        counter_value += b.counters[counter_id];
      });
      const char* counter_name = name_of_counter((counter_id_type)counter_id);
      std::cout << counter_name << " " << counter_value << std::endl;
    }
    std::cout << "launch_duration " << launch_duration << std::endl;
    double cumulated_time = launch_duration * data::perworker::get_nb_workers();
    double total_idle_time = 0.0;
    all_total_idle_time.for_each([&] (int, double& d) {
      total_idle_time += d;
    });
    double relative_idle = total_idle_time / cumulated_time;
    double utilization = 1.0 - relative_idle;
    std::cout << "total_idle_time " << total_idle_time << std::endl;
    std::cout << "utilization " << utilization << std::endl;
  }
  
};
  
template <bool enabled>
data::perworker::array<typename stats_base<enabled>::private_counters> stats_base<enabled>::all_counters;
  
template <bool enabled>
typename stats_base<enabled>::time_point_type stats_base<enabled>::enter_launch_time;
  
template <bool enabled>
double stats_base<enabled>::launch_duration;
  
template <bool enabled>
data::perworker::array<double> stats_base<enabled>::all_total_idle_time;
  
#ifdef HEARTBEAT_ENABLE_STATS
using stats = stats_base<true>;
#else
using stats = stats_base<false>;
#endif
  
} // end namespace

#endif /* _HEARTBEAT_STATS_H_ */
