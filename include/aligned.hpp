#include <array>
#include <type_traits>
#include <assert.h>

#ifndef _HEARTBEAT_ALIGNED_H_
#define _HEARTBEAT_ALIGNED_H_

namespace heartbeat {
namespace data {
    
template <class Item, int capacity>
class cache_aligned_fixed_capacity_array {
private:
  
  static constexpr int cache_align_szb = 128;
  static constexpr int item_szb = sizeof(Item);
  
  using aligned_item_type = typename std::aligned_storage<item_szb, cache_align_szb>::type;
  
  aligned_item_type items[capacity];
  
  Item& at(std::size_t i) {
    assert(i >= 0);
    assert(i < capacity);
    return *reinterpret_cast<Item*>(items + i);
  }
  
public:
  
  Item& operator[](std::size_t i) {
    return at(i);
  }
  
  std::size_t size() const {
    return capacity;
  }
  
  template <class Body>
  void for_each(const Body& f) {
    for (int i = 0; i < size(); i++) {
      f(i, at(i));
    }
  }
  
};
  
} // end namespace
} // end namespace

#endif /*! _HEARTBEAT_ALIGNED_H_ */
