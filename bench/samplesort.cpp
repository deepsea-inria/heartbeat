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
#include "sequence.hpp"
#include "samplesort.hpp"

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

using value_type = unsigned int;
  
int m = 1;

void generate(size_t _nb, parray<value_type>& dst) {
  long n = _nb * m;
  //  std::cerr << "Size: " << n << "\n";
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
    auto p = std::less<value_type>();
    heartbeat::launch_interpreter<heartbeatbench::sampleSort<value_type, typeof(p), intT>>(a.begin(), a.size(), p);
    std::sort(b.begin(), b.end(), p);
    return same_sequence(a.cbegin(), a.cend(), b.cbegin(), b.cend());
  }
  
};

#endif
  
} // end namespace

/*---------------------------------------------------------------------*/
/* Benchmarking */

#include "readinputbinary.hpp" 
#include "sampleSort.h"

namespace pbbs {

template <class Item>
using parray = sptl::parray<Item>;

template <class Item, class Compare_fn>
void benchmark(parray<Item>& xs, const Compare_fn& compare_fn) {
  deepsea::cmdline::dispatcher d;
  d.add("heartbeat", [&] {
    heartbeat::launch_interpreter<heartbeatbench::sampleSort<Item, Compare_fn, intT>>(xs.begin(), (int)xs.size(), compare_fn);
  });
  d.add("pbbs", [&] {
    heartbeatbench::run_and_report_elapsed_time([&] {
      sampleSort(xs.begin(), (int)xs.size(), compare_fn);
    });
  });
  d.dispatch("algorithm");
  if (deepsea::cmdline::parse_or_default_bool("check", false)) {
    parray<Item> ys(xs);
    std::sort(ys.begin(), ys.end(), compare_fn);
    assert(sptl::same_sequence(xs.begin(), xs.end(), ys.begin(), ys.end()));
  }
}

void benchmark(std::string infile) {
  deepsea::cmdline::dispatcher d;
  d.add("double", [&] {
    auto xs = sptl::read_from_file<parray<double>>(infile);
    benchmark(xs, std::less<double>());
  });
  d.add("int", [&] {
    auto xs = sptl::read_from_file<parray<int>>(infile);
    benchmark(xs, std::less<int>());
  });
  d.add("string", [&] {
    parray<char*> xs = sptl::read_from_file<parray<char*>>(infile);
    benchmark(xs, [&] (char* a, char* b) {
      return std::strcmp(a, b) < 0;
    });
    for (int i = 0; i < xs.size(); i++) {
      delete [] xs[i];
    }
  }); 
  d.dispatch("type");
}

} // end namespace

int main(int argc, char** argv) {
  heartbeatbench::initialize(argc, argv);
  heartbeatbench::SSORT_THR = cmdline::parse_or_default_int("threshold", heartbeatbench::SSORT_THR);
  std::string infile = deepsea::cmdline::parse_or_default_string("infile", "");
  if (infile != "") {
    pbbs::benchmark(infile);
  } else { // run unit tests
#ifdef TEST
    sptl::m = cmdline::parse_or_default("m", sptl::m);
    int nb_tests = cmdline::parse_or_default_int("nb_tests", 1000);
    checkit<sptl::sorted_property>(nb_tests, "samplesort is correct");
#endif
  }
  return 0;
}
