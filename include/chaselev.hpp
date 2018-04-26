#include <atomic>

#ifndef _HEARTBEAT_CHASELEV_H_
#define _HEARTBEAT_CHASELEV_H_

namespace heartbeat {
namespace sched {

static
vertex* const STEAL_RES_EMPTY = (vertex*) 0;
  
static
vertex* const STEAL_RES_ABORT = (vertex*) 1;

class chase_lev_deque {
public:

  using buffer_type = std::atomic<vertex*>*;

  volatile buffer_type  buf;         // deque contents
  std::atomic<int64_t> capacity;  // maximum number of elements
  std::atomic<int64_t> bottom;    // index of the first unused cell
  std::atomic<int64_t> top;       // index of the last used cell

  static
  vertex* cb_get(buffer_type buf, int64_t capacity, int64_t i)  {
    return buf[i % capacity].load();
  }
  
  static
  void cb_put (buffer_type buf, int64_t capacity, int64_t i, vertex* x) {
    buf[i % capacity].store(x);
  }
  
  buffer_type new_buffer(size_t capacity) {
    return new std::atomic<vertex*>[capacity];
  }
  
  void delete_buffer(buffer_type buf) {
    delete [] buf;
  }
  
  buffer_type grow (buffer_type old_buf,
                    int64_t old_capacity,
                    int64_t new_capacity,
                    int64_t b,
                    int64_t t) {
    chase_lev_deque::buffer_type new_buf = new_buffer(new_capacity);
    for (int64_t i = t; i < b; i++) {
      cb_put(new_buf, new_capacity, i, cb_get(old_buf, old_capacity, i));
    }
    return new_buf;
  }
  
  bool cas_top (int64_t old_val, int64_t new_val) {
    int64_t ov = old_val;
    return top.compare_exchange_strong(ov, new_val);
  }
  
  chase_lev_deque() : buf(nullptr), bottom(0l), top(0l) {
    capacity.store(0l);
  }
  
  void init(int64_t init_capacity) {
    capacity.store(init_capacity);
    buf = new_buffer(capacity.load());
    bottom.store(0l);
    top.store(0l);
    for (int64_t i = 0; i < capacity.load(); i++) {
      buf[i].store(nullptr); // optional
    }
  }
  
  void destroy() {
    assert (bottom.load() - top.load() == 0); // maybe wrong
    delete_buffer(buf);
  }
  
  void push_back(vertex* item) {
    int64_t b = bottom.load();
    int64_t t = top.load();
    if (b-t >= capacity.load() - 1) {
      chase_lev_deque::buffer_type old_buf = buf;
      int64_t old_capacity = capacity.load();
      int64_t new_capacity = capacity.load() * 2;
      buf = grow (old_buf, old_capacity, new_capacity, b, t);
      capacity.store(new_capacity);
      // UNSAFE! delete old_buf;
    }
    cb_put (buf, capacity.load(), b, item);
    // requires fence store-store
    bottom.store(b + 1);
  }
  
  vertex* pop_front() {
    int64_t t = top.load();
    // requires fence load-load
    int64_t b = bottom.load();
    // would need a fence load-load if were to read the pointer on deque->buf
    if (t >= b) {
      return STEAL_RES_EMPTY;
    }
    vertex* item = cb_get(buf, capacity.load(), t);
    // requires fence load-store
    if (! cas_top(t, t + 1)) {
      return STEAL_RES_ABORT;
    }
    return item;
  }

  vertex* pop_back() {
    int64_t b = bottom.load() - 1;
    bottom.store(b);
    // requires fence store-load
    int64_t t = top.load();
    if (b < t) {
      bottom.store(t);
      return nullptr;
    }
    vertex* item = cb_get(buf, capacity.load(), b);
    if (b > t) {
      return item;
    }
    if (! cas_top(t, t + 1)) {
      item = nullptr;
    }
    bottom.store(t + 1);
    return item;
  }

  size_t nb_threads() {
    return (size_t)bottom.load() - top.load();
  }
  
  bool empty() {
    return nb_threads() < 1;
  }
  
};
  
  
} // end namespace
} // end namespace

#endif /*! _HEARTBEAT_CHASELEV_H_ */
