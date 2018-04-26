
#include <iostream>
#include <chrono>
#include <array>
#include <fstream>
#include <iostream>
#include <istream>
#include <ostream>
#include <sstream>

#include "heartbeatbench.hpp"
#include "sequence.hpp"
#include "test.hpp"
#include "prandgen.hpp"

#include "sequence.h"

namespace sched = heartbeat::sched;
namespace cmdline = deepsea::cmdline;
namespace dsl = heartbeat::edsl;

int scale_factor = 1;

/*---------------------------------------------------------------------*/
/* Benchmark for reduce */

void bench_reduce() {
  using intT = int;
  using value_type = int;
  intT n = cmdline::parse<intT>("n");
  bool check = cmdline::parse_or_default("check", false);
  std::string algorithm = cmdline::parse<std::string>("algorithm");
  value_type* input = malloc_array<value_type>(n);
  value_type result = (value_type)0;
  for (intT i = 0; i < n; i++) {
    input[i] = (value_type)(i % 1024);
  }
  auto f = pbbs::utils::addF<value_type>();
  auto g = sequence::getA<value_type,intT>(input);
  if (algorithm == "heartbeat") {
    heartbeat::launch_interpreter<sequence::reduce<value_type,intT,typeof(f),typeof(g)>>(0, n, f, g, &result);
  } else if (algorithm == "heartbeat_sequential") {
    heartbeat::launch_interpreter<sequence::reduceSerial<value_type,intT,typeof(f),typeof(g)>>(0, n, f, g, &result);
  } else if (algorithm == "sequential") {
    heartbeatbench::run_and_report_elapsed_time([&] {
      result = pbbs::sequence::reduceSerial<value_type>(0, n, f, g);
    });
  } else if (algorithm == "pbbs") {
    heartbeatbench::run_and_report_elapsed_time([&] {
      result = pbbs::sequence::reduce<value_type>(0, n, f, g);
    });
  }
  if (check) {
#ifndef NDEBUG
    value_type result2 = pbbs::sequence::reduceSerial<value_type>(0, n, f, g);
#endif
    assert(result == result2);
  }
  free(input);
}

/*---------------------------------------------------------------------*/
/* Benchmark for max_index */

void bench_max_index() {
  using intT = int;
  using value_type = int;
  intT n = cmdline::parse<intT>("n");
  bool check = cmdline::parse_or_default("check", false);
  std::string algorithm = cmdline::parse<std::string>("algorithm");
  value_type* input = malloc_array<value_type>(n);
  value_type result = (value_type)0;
  for (intT i = 0; i < n; i++) {
    input[i] = (value_type)i;
  }
  auto f = pbbs::utils::addF<value_type>();
  auto g = sequence::getA<value_type,intT>(input);
  if (algorithm == "heartbeat") {
    heartbeat::launch_interpreter<sequence::maxIndex<value_type,intT,typeof(f),typeof(g)>>(0, n, f, g, &result);
  } else if (algorithm == "heartbeat_sequential") {
    heartbeat::launch_interpreter<sequence::maxIndexSerial<value_type,intT,typeof(f),typeof(g)>>(0, n, f, g, &result);
  } else if (algorithm == "sequential") {
    heartbeatbench::run_and_report_elapsed_time([&] {
      result = pbbs::sequence::maxIndexSerial<value_type>(0, n, f, g);
    });
  } else if (algorithm == "pbbs") {
    heartbeatbench::run_and_report_elapsed_time([&] {
      result = pbbs::sequence::maxIndex<value_type>(0, n, f, g);
    });
  }
  if (check) {
#ifndef NDEBUG
    value_type result2 = pbbs::sequence::maxIndexSerial<value_type>(0, n, f, g);
#endif
    assert(result == result2);
  }
  free(input);
}

/*---------------------------------------------------------------------*/
/* Benchmark for scan */

template <class Num>
Num abs(Num x) {
  return x >= 0 ? x : -1 * x;
}
template <class Num, class Item>
Num count_diffs(Num n, Num m, Item* trusted, Item* untrusted) {
  if (n != m) {
    return abs(n - m);
  }
  Num nb_diffs = 0;
  for (Num i = 0; i < n; i++) {
    Item& t = trusted[i];
    Item& u = untrusted[i];
    if (t != u) {
      nb_diffs++;
    }
  }
  return nb_diffs;
}

template <class Num, class Item>
Num count_diffs(Num n, Item* trusted, Item* untrusted) {
  return count_diffs(n, n, trusted, untrusted);
}

void bench_scan() {
  using intT = int;
  using value_type = int;
  intT n = cmdline::parse<intT>("n");
  bool check = cmdline::parse_or_default("check", false);
  std::string algorithm = cmdline::parse<std::string>("algorithm");
  bool inclusive = cmdline::parse_or_default("inclusive", false);
  bool back = cmdline::parse_or_default("back", false);
  value_type* input = malloc_array<value_type>(n);
  value_type* output = malloc_array<value_type>(n);
  value_type zero = 0;
  value_type result = 0;
  for (intT i = 0; i < n; i++) {
    input[i] = 1;
    output[i] = -1234;
  }
  auto f = pbbs::utils::addF<value_type>();
  auto g = sequence::getA<value_type,intT>(input);
  if (algorithm == "heartbeat") {
    heartbeat::launch_interpreter<sequence::scan<value_type,intT,typeof(f),typeof(g)>>(output, 0, n, f, g, zero, inclusive, back, &result);
  } else if (algorithm == "heartbeat_sequential") {
    heartbeat::launch_interpreter<sequence::scanSerial<value_type,intT,typeof(f),typeof(g)>>(output, 0, n, f, g, zero, inclusive, back, &result);
  } else if (algorithm == "sequential") {
    heartbeatbench::run_and_report_elapsed_time([&] {
      result = pbbs::sequence::scanSerial<value_type>(output, (intT) 0, n, f, g, zero, inclusive, back);
    });
  } else if (algorithm == "pbbs") {
    heartbeatbench::run_and_report_elapsed_time([&] {
      result = pbbs::sequence::scan<value_type>(output, (intT) 0, n, f, g, zero, inclusive, back);
    });
  }
  if (check) {
#ifndef NDEBUG
    value_type* output2 = malloc_array<value_type>(n);
    value_type result2 = pbbs::sequence::scanSerial(output2, (intT) 0, n, f, g, zero, inclusive, back);
    auto nb_diffs = count_diffs(n, output2, output);
    assert(nb_diffs == 0);
    free(output2);
#endif
    assert(result == result2);
  }
  free(input);
  free(output);
}

/*---------------------------------------------------------------------*/
/* Benchmark for pack */

template <class T>
_seq<T> from_pbbs(pbbs::_seq<T> s) {
  return _seq<T>(s.A, s.n);
}

void bench_pack() {
  using intT = int;
  using value_type = int;
  intT n = std::max(2, cmdline::parse<intT>("n"));
  intT m = cmdline::parse_or_default("m", 3);
  bool check = cmdline::parse_or_default("check", false);
  std::string algorithm = cmdline::parse<std::string>("algorithm");
  value_type* input = malloc_array<value_type>(n);
  bool* flags = malloc_array<bool>(n);
  _seq<value_type> output;
  for (intT i = 0; i < n; i++) {
    input[i] = (value_type)i;
    flags[i] = rand() % m == 0;
  }
  auto f = sequence::getA<value_type,intT>(input);
  if (algorithm == "heartbeat") {
    heartbeat::launch_interpreter<sequence::pack<value_type, intT, typeof(f)>>((value_type*)nullptr, flags, (intT)0, n, f, &output);
  } else if (algorithm == "heartbeat_sequential") {
    heartbeat::launch_interpreter<sequence::packSerial<value_type, intT, typeof(f)>>((value_type*)nullptr, flags, (intT)0, n, f, &output);
  } else if (algorithm == "sequential") {
    heartbeatbench::run_and_report_elapsed_time([&] {
      output = from_pbbs(pbbs::sequence::packSerial((value_type*)nullptr, flags, (intT)0, n, f));
    });
  } else if (algorithm == "pbbs") {
    heartbeatbench::run_and_report_elapsed_time([&] {
      output = from_pbbs(pbbs::sequence::pack((value_type*)nullptr, flags, (intT)0, n, f));
    });
  }
  if (check) {
#ifndef NDEBUG
    _seq<value_type> output2 = from_pbbs(pbbs::sequence::packSerial((value_type*)nullptr, flags, (intT)0, n, f));
    assert(output.n == output2.n);
    for (int i = 0; i < output.n; i++) {
      auto u = output.A[i];
      auto t = output2.A[i];
      assert(u == t);
    }
    output2.del();
#endif
  }
  free(input);
  output.del();
}

/*---------------------------------------------------------------------*/
/* Benchmark for filter */

template <class ET, class intT, class PRED>
_seq<ET> filterSerial(ET* In, intT n, PRED p) {
  intT m = 0;
  for (int i = 0; i < n; i++) {
    if (p(In[i])) {
      m++;
    }
  }
  ET* A = malloc_array<ET>(m);
  _seq<ET> result(A, m);
  m = 0;
  for (int i = 0; i < n; i++) {
    if (p(In[i])) {
      A[m++] = In[i];
    }
  }
  return result;
}

void bench_filter() {
  using intT = int;
  using value_type = int;
  intT n = std::max(2, cmdline::parse<intT>("n"));
  intT m = cmdline::parse_or_default("m", 3);
  bool check = cmdline::parse_or_default("check", false);
  std::string algorithm = cmdline::parse<std::string>("algorithm");
  value_type* input = malloc_array<value_type>(n);
  bool* flags = malloc_array<bool>(n);
  _seq<value_type> output;
  for (intT i = 0; i < n; i++) {
    input[i] = (value_type)i;
    flags[i] = rand() % m == 0;
  }
  auto p = [&] (value_type v) {
    return flags[v];
  };
  if (algorithm == "heartbeat") {
    heartbeat::launch_interpreter<sequence::filter<value_type,intT,typeof(p)>>(input, n, p, &output);
  } else if (algorithm == "heartbeat_sequential") {
    assert(false);
  } else if (algorithm == "sequential") {
    heartbeatbench::run_and_report_elapsed_time([&] {
      output = filterSerial(input, n, p);
    });
  } else if (algorithm == "pbbs") {
    heartbeatbench::run_and_report_elapsed_time([&] {
      output = from_pbbs(pbbs::sequence::filter(input, n, p));
    });
  }
  if (check) {
#ifndef NDEBUG
    _seq<value_type> output2;
    output2 = filterSerial(input, n, p);
    assert(output.n == output2.n);
    for (int i = 0; i < output.n; i++) {
      auto u = output.A[i];
      auto t = output2.A[i];
      assert(t == u);
    }
    output2.del();
#endif
  }
  free(input);
  output.del();
}

/*---------------------------------------------------------------------*/
/* Quickcheck unit testing */

namespace pasl {
namespace pctl {
  
template <class Item1, class Item2>
std::ostream& operator<<(std::ostream& out, const std::pair<parray<Item1>, parray<Item2>>& p) {
  out << p.first;
  out << "\n";
  out << p.second;
  return out;
}

template <class Container>
std::ostream& operator<<(std::ostream& out, const container_wrapper<Container>& c) {
  out << c.c;
  return out;
}

void generate(size_t _nb, parray<intT>& dst) {
  dst = prandgen::gen_integ_parray(_nb, 0, 32);
}

void generate(size_t nb, container_wrapper<parray<intT>>& c) {
  nb *= scale_factor;
  generate(nb, c.c);
}
  
void generate(size_t _nb, std::pair<parray<intT>,parray<bool>>& dst) {
  generate(_nb, dst.first);
  dst.second = prandgen::gen_parray<bool>(_nb, [&] (long, unsigned int x) {
    return x % 2 == 0;
  });
}

void generate(size_t nb, container_wrapper<std::pair<parray<intT>,parray<bool>>>& c) {
  nb *= scale_factor; 
  generate(nb, c.c);
}
  
using int_array_wrapper = container_wrapper<parray<intT>>;
  
class consistent_copies_property : public quickcheck::Property<int_array_wrapper> {
public:
  
  using value_type = intT;
  
  bool holdsFor(const int_array_wrapper& _in) {
    int_array_wrapper in(_in);
    parray<intT> output(in.c.size());
    heartbeat::launch_interpreter<sequence::copy<intT*, intT*>>(in.c.begin(), in.c.end(), output.begin());
    return count_diffs(in.c.size(), in.c.begin(), output.begin()) == 0;
  }
  
};
  
class consistent_reductions_property : public quickcheck::Property<int_array_wrapper> {
public:
  
  using value_type = intT;
  
  bool holdsFor(const int_array_wrapper& _in) {
    int_array_wrapper in(_in);
    value_type* input = in.c.begin();
    intT n = (intT)in.c.size();
    value_type result = 0;
    auto f = pbbs::utils::addF<value_type>();
    auto g = sequence::getA<value_type,intT>(input);
    heartbeat::launch_interpreter<sequence::reduce<value_type,intT,typeof(f),typeof(g)>>(0, n, f, g, &result);
    value_type result2 = pbbs::sequence::reduceSerial<value_type>(0, n, f, g);
    return result == result2;
  }
  
};
  
class consistent_max_indexes_property : public quickcheck::Property<int_array_wrapper> {
public:
  
  using value_type = intT;
  
  bool holdsFor(const int_array_wrapper& _in) {
    int_array_wrapper in(_in);
    intT n = (intT)in.c.size();
    value_type result = 0;
    auto f = std::greater<intT>();
    auto g = sequence::getA<value_type,intT>(in.c.begin());
    heartbeat::launch_interpreter<sequence::maxIndex<value_type,intT,typeof(f),typeof(g)>>(0, n, f, g, &result);
    value_type result2 = pbbs::sequence::maxIndexSerial<value_type>(0, n, f, g);
    return result == result2;
  }
  
};

class consistent_scans_property : public quickcheck::Property<int_array_wrapper> {
public:
  
  using value_type = intT;
  
  bool holdsFor(const int_array_wrapper& _in) {
    int_array_wrapper in(_in);
    value_type* input = in.c.begin();
    intT n = (intT)in.c.size();
    value_type zero = 0;
    value_type result = 0;
    auto f = pbbs::utils::addF<value_type>();
    auto g = sequence::getA<value_type,intT>(input);
    parray<intT> output(n);
    bool inclusive = cmdline::parse_or_default("inclusive", false);
    bool back = cmdline::parse_or_default("back", false);
    heartbeat::launch_interpreter<sequence::scan<value_type,intT,typeof(f),typeof(g)>>(output.begin(), 0, n, f, g, zero, inclusive, back, &result);
    parray<intT> output2(n);
    value_type result2 = pbbs::sequence::scanSerial(output2.begin(), (intT) 0, n, f, g, zero, inclusive, back);
    return result == result2 && count_diffs(n, output2.begin(), output2.begin()) == 0;
  }
  
};
  
using pack_wrapper = container_wrapper<std::pair<parray<intT>, parray<bool>>>;
  
class consistent_packs_property : public quickcheck::Property<pack_wrapper> {
public:
  
  using value_type = intT;
  
  bool holdsFor(const pack_wrapper& _in) {
    pack_wrapper in(_in);
    _seq<value_type> output;
    intT n = (intT)in.c.first.size();
    assert(in.c.first.size() == in.c.second.size());
    auto f = sequence::getA<value_type,intT>(in.c.first.begin());
    auto flags = in.c.second.begin();
    heartbeat::launch_interpreter<sequence::pack<value_type, intT, typeof(f)>>((value_type*)nullptr, flags, (intT)0, n, f, &output);
    auto output2 = from_pbbs(pbbs::sequence::packSerial((value_type*)nullptr, flags, (intT)0, n, f));
    auto nb_diffs = count_diffs(output2.n, output.n, output2.A, output.A);
    output.del();
    output2.del();
    return nb_diffs == 0;
  }
  
};
  
class consistent_filters_property : public quickcheck::Property<int_array_wrapper> {
public:
  
  using value_type = intT;
  
  bool holdsFor(const int_array_wrapper& _in) {
    int_array_wrapper in(_in);
    _seq<value_type> output;
    intT n = (intT)in.c.size();
    intT k = cmdline::parse_or_default("keep_probability", 33);
    auto p = [&] (value_type v) {
      return v % k == 0;
    };
    heartbeat::launch_interpreter<sequence::filter<value_type,intT,typeof(p)>>(in.c.begin(), n, p, &output);
    auto output2 = from_pbbs(pbbs::sequence::filter(in.c.begin(), n, p));
    auto nb_diffs = count_diffs(output2.n, output.n, output2.A, output.A);
    output.del();
    output2.del();
    return nb_diffs == 0;
  }
  
};
  
class consistent_filterDPSs_property : public quickcheck::Property<int_array_wrapper> {
public:
  
  using value_type = intT;
  
  bool holdsFor(const int_array_wrapper& _in) {
    int_array_wrapper in(_in);
    intT n = (intT)in.c.size();
    intT k = cmdline::parse_or_default("keep_probability", 33);
    auto p = [&] (value_type v) {
      return v % k == 0;
    };
    auto output2 = from_pbbs(pbbs::sequence::filter(in.c.begin(), n, p));
    parray<value_type> output(n);
    intT dest;
    heartbeat::launch_interpreter<sequence::filterDPS<value_type,intT,typeof(p)>>(in.c.begin(), output.begin(), n, p, &dest);
    auto nb_diffs = count_diffs((intT)output2.n, dest, output2.A, output.begin());
    output2.del();
    return nb_diffs == 0;
  }
  
};

} // end namespace
} // end namespace

/*---------------------------------------------------------------------*/
/* Main routine */

void benchmark() {
  cmdline::dispatcher d;
  d.add("reduce", [&] {
    bench_reduce();
  });
  d.add("max_index", [&] {
    bench_max_index();
  });
  d.add("scan", [&] {
    bench_scan();
  });
  d.add("pack", [&] {
    bench_pack();
  });
  d.add("filter", [&] {
    bench_filter();
  });
  d.dispatch("operation");
}

void test() {
  scale_factor = deepsea::cmdline::parse_or_default("scale_factor", scale_factor);
  int nb_tests = deepsea::cmdline::parse_or_default_int("nb_tests", 1000);
  cmdline::dispatcher d;
  d.add("copy", [&] {
    checkit<pasl::pctl::consistent_copies_property>(nb_tests, "copy is correct");
  });
  d.add("reduce", [&] {
    checkit<pasl::pctl::consistent_reductions_property>(nb_tests, "reduce is correct");
  });
  d.add("max_index", [&] {
    checkit<pasl::pctl::consistent_max_indexes_property>(nb_tests, "max_index is correct");
  });
  d.add("scan", [&] {
    checkit<pasl::pctl::consistent_scans_property>(nb_tests, "scan is correct");
  });
  d.add("pack", [&] {
    checkit<pasl::pctl::consistent_packs_property>(nb_tests, "pack is correct");
  });
  d.add("filter", [&] {
    checkit<pasl::pctl::consistent_filters_property>(nb_tests, "filter is correct");
  });
  d.add("filterDPS", [&] {
    checkit<pasl::pctl::consistent_filterDPSs_property>(nb_tests, "filter is correct");
  });
  d.dispatch("operation");
}

int main(int argc, char** argv) {
  heartbeatbench::initialize(argc, argv);
  sequence::initialize();
  cmdline::dispatcher d;
  d.add("benchmark", [&] {
    benchmark();
  });
  d.add("test", [&] {
    test();
  });
  d.dispatch_or_default("mode", "benchmark");
  return 0;
}
