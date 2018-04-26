#include <iostream>
#include <chrono>
#include <array>
#include <fstream>
#include <iostream>
#include <istream>
#include <ostream>
#include <sstream>
#include <algorithm>

#include "heartbeatbench.hpp"
#include "blockradixsort.hpp"
#include "readinputbinary.hpp" 
#include "blockRadixSort.h"

namespace pbbs {

template <class Item>
using parray = sptl::parray<Item>;

void benchmark(parray<int>& x) {
  bool should_check = deepsea::cmdline::parse_or_default_bool("check", false);
  parray<int> ref;
  if (should_check) {
    ref = x;
  }
  deepsea::cmdline::dispatcher d;
  d.add("heartbeat", [&] {
    heartbeat::launch_interpreter<heartbeatbench::integerSort<int>>(x.begin(), (int)x.size());
  });
  d.add("pbbs", [&] {
    heartbeatbench::run_and_report_elapsed_time([&] {
      integerSort<int>(&x[0], (int)x.size());
    });
  });
  d.dispatch("algorithm");
  if (should_check) {
    std::sort(ref.begin(), ref.end());
    auto it_ref = ref.begin();
    for (auto it = x.begin(); it != x.end(); it++) {
      if (*it != *it_ref) {
        std::cerr << "bogus result" << std::endl;
        exit(0);
      }
      it_ref++;
    }
  }
}

void benchmark(parray<std::pair<int, int>>& x) {
  deepsea::cmdline::dispatcher d;
  d.add("heartbeat", [&] {
    heartbeat::launch_interpreter<heartbeatbench::integerSortPair<int,intT>>(x.begin(), (int)x.size());
  });
  d.add("pbbs", [&] {
    heartbeatbench::run_and_report_elapsed_time([&] {
      integerSort<int>(&x[0], (int)x.size());
    });
  });
  d.dispatch("algorithm");
}

void benchmark(std::string infile) {
  deepsea::cmdline::dispatcher d;
  d.add("int", [&] {
    parray<int> x = sptl::read_from_file<parray<int>>(infile);      
    benchmark(x);
  });
  d.add("pair_int_int", [&]  {
    parray<std::pair<int, int>> x = sptl::read_from_file<parray<std::pair<int, int>>>(infile);            
    benchmark(x);
  });
  d.dispatch("type");
}

} // end namespace

#ifdef TEST
#include "test.hpp"
#include "prandgen.hpp"
#include "sequencedata.hpp"
#endif

namespace sched = heartbeat::sched;
namespace cmdline = deepsea::cmdline;
namespace dsl = heartbeat::edsl;

namespace sptl {

#ifdef TEST

/*---------------------------------------------------------------------*/
/* Quickcheck IO */

template <class Container>
std::ostream& operator<<(std::ostream& out, const container_wrapper<Container>& c) {
  out << c.c;
  return out;
}

/*---------------------------------------------------------------------*/
/* Quickcheck generators */


using value_type = intT;
  
int m = 1;

void generate(size_t _nb, parray<value_type>& dst) {
  long n = _nb * m;
  std::cerr << "Size: " << n << "\n";
      dst = sequencedata::exp_dist<value_type>(0L, (value_type)n);
      return;
  int r = quickcheck::generateInRange(0, 2); // currently something is wrong with exp_dist
  //  std::cerr << r << "\n";
  if (r == 0) {
    dst = sequencedata::rand_int_range((value_type)0, (value_type)n, (value_type)INT_MAX);
  } else if (r == 1) {
    int m = quickcheck::generateInRange(0, 1 << 10);
    dst = sequencedata::almost_sorted<value_type>(0L, n, m);
  } else if (r == 2) {
    value_type x = (value_type)quickcheck::generateInRange(0, INT_MAX);
    dst = sequencedata::all_same(n, x);
  } else {
    dst = sequencedata::exp_dist<value_type>(0L, (value_type)n);
  }
}

void generate(size_t nb, container_wrapper<parray<value_type>>& c) {
  generate(nb, c.c);
}


/*---------------------------------------------------------------------*/
/* Quickcheck properties */

using parray_wrapper = container_wrapper<parray<value_type>>;

class sorted_property : public quickcheck::Property<parray_wrapper> {
public:
  
  bool holdsFor(const parray_wrapper& _in) {
    parray<value_type> a = _in.c;
    parray<value_type> b = _in.c;
    heartbeat::launch_interpreter<heartbeatbench::integerSort<int>>(a.begin(), (int)a.size());
    std::sort(b.begin(), b.end(), std::less<value_type>());
    return same_sequence(a.cbegin(), a.cend(), b.cbegin(), b.cend());
  }
  
};

#endif

} // end namespace

int main(int argc, char** argv) {
  heartbeatbench::initialize(argc, argv);
  std::string infile = deepsea::cmdline::parse_or_default_string("infile", "");
  if (infile == "") {
    int nb_tests = cmdline::parse_or_default_int("nb_tests", 1000);
#ifdef TEST
    checkit<pasl::pctl::sorted_property>(nb_tests, "radixsort is correct");
#endif
    return 0;
  }
  pbbs::benchmark(infile);
  return 0;
}
