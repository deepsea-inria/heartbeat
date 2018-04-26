
#include "cycles.hpp"
#include "perworker.hpp"

#ifndef _HEARTBEAT_FUEL_H_
#define _HEARTBEAT_FUEL_H_

namespace heartbeat {
namespace fuel {

using check_type = enum {
  check_yes_promote,
  check_no_promote,
  check_suspend
};

data::perworker::array<uint64_t> target_for_next_promotion;

uint64_t promotion_threshold = 0;

check_type check(uint64_t now) {
  uint64_t next = target_for_next_promotion.mine();
  if (now < next) {
    return check_no_promote;
  }
  target_for_next_promotion.mine() = now + promotion_threshold;
  return check_yes_promote;
}

void initialize(double cpu_freq_ghz, double kappa_nsec) {
  double cycles_per_nsec = cpu_freq_ghz;
  promotion_threshold = (uint64_t) (cycles_per_nsec * kappa_nsec);
}

void initialize_worker() {
  target_for_next_promotion.mine() = cycles::now();
}

} // end namespace
} // end namespace

#endif /*! _HEARTBEAT_FUEL_H_ */
