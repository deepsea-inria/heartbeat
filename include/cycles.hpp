
#ifndef _HEARTBEAT_CYCLES_H_
#define _HEARTBEAT_CYCLES_H_

namespace heartbeat {
namespace cycles {

namespace {
  
static inline
uint64_t rdtsc() {
  unsigned int hi, lo;
  __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
  return  ((uint64_t) lo) | (((uint64_t) hi) << 32);
}

static inline
void rdtsc_wait(uint64_t n) {
  const uint64_t start = rdtsc();
  while (rdtsc() < (start + n)) {
    __asm__("PAUSE");
  }
}
  
} // end namespace
  
static inline
uint64_t diff(uint64_t start, uint64_t finish) {
  return finish - start;
}

static inline
uint64_t now() {
  return rdtsc();
}

static inline
uint64_t since(uint64_t start) {
  return diff(start, now());
}

static inline
void spin_for(uint64_t nb_cycles) {
  rdtsc_wait(nb_cycles);
}
  
} // end namespace
} // end namespace

#endif /*! _HEARTBEAT_CYCLES_H_ */
