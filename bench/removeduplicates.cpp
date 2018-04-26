/*!
 * \file deterministichash_bench.cpp
 * \brief Benchmarking script for parallel hash table
 * \date 2016
 * \copyright COPYRIGHT (c) 2015 Umut Acar, Arthur Chargueraud, and
 * Michael Rainey. All rights reserved.
 * \license This project is released under the GNU Public License.
 *
 */

#include <math.h>
#include <functional>
#include <stdlib.h>

#include "heartbeatbench.hpp"
#include "deterministichash.hpp"
#include "readinputbinary.hpp" 
#include "deterministicHash.h"
#undef parallel_for
/***********************************************************************/

/*---------------------------------------------------------------------*/

template <class Item>
using parray = sptl::parray<Item>;

template <class Item>
void heartbeat_bench(parray<Item>& x) {
  std::string lib_type = deepsea::cmdline::parse_or_default_string("algorithm", "heartbeat");
  if (lib_type == "pbbs") {
    heartbeatbench::run_and_report_elapsed_time([&] {
        pbbs::_seq<Item> res = pbbs::removeDuplicates(pbbs::_seq<Item>(&x[0], x.size()));
    });
  } else {
    auto ptr = x.begin();
    auto n = x.size();
    heartbeat::launch_interpreter_via_lambda([=] (heartbeat::edsl::pcfg::stack_type st) {
      auto ty = heartbeat::edsl::pcfg::cactus::Parent_link_sync;
      auto s = _seq<Item>(ptr, n);
      _seq<Item> res;
      return heartbeatbench::removeDuplicates(st, ty, s, &res);
    });
  }
}

int main(int argc, char** argv) {
  heartbeatbench::initialize(argc, argv);
  std::string infile = deepsea::cmdline::parse_or_default_string("infile", "");
  if (infile != "") {
    deepsea::cmdline::dispatcher d;
    d.add("array_int", [&] {
      parray<int> x = sptl::read_from_file<parray<int>>(infile);
      heartbeat_bench(x);
    });
    d.add("array_string", [&] {
      parray<char*> x = sptl::read_from_file<parray<char*>>(infile);
      heartbeat_bench(x);
    });
    d.add("array_pair_string_int", [&] {
      parray<std::pair<char*, int>*> x = sptl::read_from_file<parray<std::pair<char*, int>*>>(infile);
      heartbeat_bench(x);
    });
    d.dispatch("type");
  }
  return 0;
}

/***********************************************************************/
