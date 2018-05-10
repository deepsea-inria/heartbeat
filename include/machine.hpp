
#include <assert.h>

#ifdef HAVE_HWLOC
#include <hwloc.h>
#endif

#ifdef TARGET_MAC_OS
#include <sys/param.h>
#include <sys/sysctl.h>
#endif

#include "cmdline.hpp"
#include "atomic.hpp"

#ifndef _HEARTBEAT_MACHINE_H_
#define _HEARTBEAT_MACHINE_H_

namespace heartbeat {
namespace machine {

namespace cmdline = deepsea::cmdline;
    
#ifdef HAVE_HWLOC
hwloc_topology_t topology;
#endif

void initialize_hwloc(int nb_workers) {
  bool numa_interleave = false;
#ifdef HAVE_HWLOC
  hwloc_topology_init(&topology);
  hwloc_topology_load(topology);
  bool numa_alloc_interleaved = (nb_workers == 0) ? false : true;
  numa_alloc_interleaved = cmdline::parse_or_default("numa_alloc_interleaved", numa_alloc_interleaved);
  if (numa_alloc_interleaved) {
    hwloc_cpuset_t all_cpus =
      hwloc_bitmap_dup(hwloc_topology_get_topology_cpuset(topology));
    int err = hwloc_set_membind(topology, all_cpus, HWLOC_MEMBIND_INTERLEAVE, 0);
    if (err < 0) {
      printf("Warning: failed to set NUMA round-robin allocation policy\n");
    }
  }
#endif
  printf("hwloc_interleave %d\n", numa_interleave);
}

void initialize_hwloc() {
  initialize_hwloc(cmdline::parse_or_default("proc", 1));
}

double cpu_frequency_ghz = 1.2;
  
void initialize_cpuinfo() {
  float cpu_frequency_mhz = 0.0;
#ifdef TARGET_LINUX
  /* Get information from /proc/cpuinfo.     *
   * cpu MHz         : <float>             # cpu frequency in MHz
   */
  FILE *cpuinfo_file = fopen("/proc/cpuinfo", "r");
  char buf[1024];
  int cache_line_szb;
  if (cpuinfo_file != NULL) {
    while (fgets(buf, sizeof(buf), cpuinfo_file) != 0) {
      sscanf(buf, "cpu MHz : %f", &(cpu_frequency_mhz));
    }
    fclose (cpuinfo_file);
  }
#endif
#ifdef TARGET_MAC_OS
  uint64_t freq = 0;
  size_t size;
  size = sizeof(freq);
  if (sysctlbyname("hw.cpufrequency", &freq, &size, NULL, 0) < 0) {
    perror("sysctl");
  }
  cpu_frequency_mhz = (float)freq / 1000000.;
#endif
  if (cpu_frequency_mhz == 0.) {
    atomic::die("Failed to read CPU frequency\n");
  }
  cpu_frequency_ghz = (double) (cpu_frequency_mhz / 1000.0);
}
  
} // end namespace
} // end namespace

#endif /*! _HEARTBEAT_MACHINE_H_ */
