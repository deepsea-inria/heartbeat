
#include <memory>

#include "forward.hpp"
#include "incounter.hpp"
#include "outset.hpp"
#include "fuel.hpp"

#ifndef _HEARTBEAT_SCHED_VERTEX_H_
#define _HEARTBEAT_SCHED_VERTEX_H_

namespace heartbeat {
namespace sched {
  
class vertex {
public:
  
  incounter_handle release_handle;
  
  bool is_suspended = false;
  
private:
  
  incounter in;

  outset out;
  
public:
  
  vertex() : in(this) {
    release_handle = in.increment(this);
  }
  
  virtual
  ~vertex() {
    assert(is_ready());
  }
  
  incounter* get_incounter() {
    return &in;
  }

#ifndef NDEBUG
  bool is_ready() {
    return ! in.is_nonzero();
  }
#endif

  void make_ready() {
    in.make_ready(this);
  }

  outset* get_outset() {
    return &out;
  }  
    
  virtual
  int nb_strands() = 0;
  
  virtual
  fuel::check_type run() = 0;
  
  virtual
  vertex_split_type split(int nb) = 0;
  
};
  
// invariant 1: each vertex is ready (incounter is zero)
// invariant 2: vertices ordered by priority
// invariant 3: vertex v0 may be the null pointer, but v1 and v2 must not
struct vertex_split_struct {
  vertex* v0;
  vertex* v1;
  vertex* v2;
};

vertex_split_type make_vertex_split(vertex* v0, vertex* v1, vertex* v2) {
  return {.v0 = v0, .v1 = v1, .v2 = v2};
}

vertex_split_type make_vertex_split(vertex* v1, vertex* v2) {
  return make_vertex_split(nullptr, v1, v2);
}
  
} // end namespace
} // end namespace

#endif /*! _HEARTBEAT_SCHED_VERTEX_H_ */
