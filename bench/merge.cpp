
#include <iostream>
#include <chrono>

#include "heartbeatbench.hpp"
#include "merge.hpp"
#include "merge.h"

namespace cmdline = deepsea::cmdline;

int main(int argc, char** argv) {
  heartbeatbench::initialize(argc, argv);
  int n = cmdline::parse<int>("n");
  int* a1 = malloc_array<int>(n);
  int* a2 = malloc_array<int>(n);
  int* a3 = malloc_array<int>(2*n);
  for (int i = 0; i < n; i++) {
    a1[i] = rand() % 1000;
    a2[i] = rand() % 1000;
  }
  auto f = [] (int x, int y) {
    return x < y;
  };
  std::string algorithm = cmdline::parse<std::string>("algorithm");
  if (algorithm == "heartbeat") {
    heartbeat::launch_interpreter<heartbeatbench::merge<int,decltype(f),int>>(a1, n, a2, n, a3, f);
  } else if (algorithm == "pbbs") {
    heartbeatbench::run_and_report_elapsed_time([&] {
      pbbs::merge(a1, n, a2, n, a3, f);
    });
  }
  free(a1);
  free(a2);
  return 0;
}
