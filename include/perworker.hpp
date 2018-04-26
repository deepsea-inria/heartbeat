
#include <atomic>

#include "aligned.hpp"

#ifndef _HEARTBEAT_PERWORKER_H_
#define _HEARTBEAT_PERWORKER_H_

namespace heartbeat {
namespace data {
namespace perworker {
  
namespace {
  
static constexpr int default_max_nb_workers = 128;
  
std::atomic<int> fresh_id(0);
  
__thread int my_id = -1;
  
} // end namespace
  
void reset() {
  fresh_id.store(0);
  my_id = -1;
}
  
int get_my_id() {
  if (my_id == -1) {
    my_id = fresh_id++;
  }
  return my_id;
}
  
int get_nb_workers() {
  return fresh_id.load();
}
  
class my_fresh_id {
public:
  
  int operator()() {
    return get_my_id();
  }
  
};
  
// later: make default for My_id the heartbeat built-in one, when it exists

template <class Item, class My_id=my_fresh_id, int max_nb_workers=default_max_nb_workers>
class array {
private:
  
  cache_aligned_fixed_capacity_array<Item, max_nb_workers> items;
  
  int get_my_id() {
    My_id my_id;
    int id = my_id();
    assert(id >= 0);
    assert(id < max_nb_workers);
    return id;
  }
  
public:
  
  array() {
    for_each([&] (int, Item& x) {
      new (&x) Item();
    });
  }
  
  array(const Item& x) {
    for_each([&] (int, Item& y) {
      new (&y) Item(x);
    });
  }
  
  array(std::initializer_list<Item> l) {
    assert(l.size() == 1);
    items.init(*(l.begin()));
  }
  
  ~array() {
    for_each([&] (int, Item& x) {
      x.~Item();
    });
  }
  
  Item& mine() {
    return items[get_my_id()];
  }
  
  Item& operator[](std::size_t i) {
    assert(i >= 0);
    assert(i < max_nb_workers);
    return items[i];
  }
  
  template <class Body>
  void for_each(const Body& f) {
    items.for_each(f);
  }
  
};
  
  
} // end namespace
} // end namespace
} // end namespace

#endif /*! _PASL_PERWORKER_H_ */
