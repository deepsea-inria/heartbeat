#include <assert.h>

#ifndef _HEARTBEAT_FORWARD_H_
#define _HEARTBEAT_FORWARD_H_

namespace heartbeat {
namespace sched {
    
class vertex;
class incounter;
class outset;
struct vertex_split_struct;

using vertex_split_type = struct vertex_split_struct;
  
void schedule(vertex* v);
void parallel_notify(outset*);
void release(vertex* v);
void suspend(vertex* v);
  
} // end namespace
} // end namespace

#endif /*! _HEARTBEAT_FORWARD_H_ */
