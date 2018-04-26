/*
 * A tagged value is a word which consists of a small integer
 * tag and a large integer value.
 *
 * The size in bits of the value is just large enough to contain a
 * (properly aligned) pointer.
 *
 * The operations typically use three types:
 *   - The type of the tag, value pair is named `Tagged`.
 *   - The type of the value is named `Value`.
 *   - The tag has type `long`.
 *
 * Currently this module supports only 64-bit words. As such,
 * the tag is a 3 bit integer and the value a 61 bit integer.
 *
 */

#include <assert.h>

#ifndef _HEARTBEAT_TAGGED_H_
#define _HEARTBEAT_TAGGED_H_

namespace heartbeat {
namespace tagged {

namespace {
// Size in bits of the tag
const long NUM_TAG_BITS = 3;
// Bit mask for the tag
const long TAG_MASK = (1<<NUM_TAG_BITS)-1;
} // end namespace

// Returns a tagged value consisting of the pair v and bits
template <typename Value, typename Tagged>
static inline Tagged create(Value v, long bits) {
  assert(sizeof(Value) == sizeof(int64_t));
  assert(sizeof(Tagged) == sizeof(int64_t));
  assert(bits <= TAG_MASK);
  union enc {
    Value         v;
    Tagged        t;
    uint64_t      bits;
  };
  enc e;
  e.v = v;
  assert((e.bits & TAG_MASK) == 0l);
  e.bits |= (uint64_t)bits;
  return e.t;
}
  
template <class Value>
Value* tag_with(Value* v, int t) {
  return create<Value*, Value*>(v, (long) t);
}

// Returns the value component of the pair
template <typename Value, typename Tagged>
static inline Value value_of(Tagged t) {
  assert(sizeof(Value) == sizeof(int64_t));
  assert(sizeof(Tagged) == sizeof(int64_t));
  union enc {
    Value         v;
    Tagged        t;
    uint64_t      bits;
  };
  enc e;
  e.t = t;
  e.bits &= ~ TAG_MASK;
  return e.v;
}
  
template <class Value>
Value* pointer_of(Value* v) {
  return value_of<Value*>(v);
}

// Returns the tag component of the pair
template <typename Value, typename Tagged>
static inline long tag_of(Tagged t) {
  assert(sizeof(Value) == sizeof(int64_t));
  assert(sizeof(Tagged) == sizeof(int64_t));
  union enc {
    Value         v;
    Tagged        t;
    uint64_t      bits;
  };
  enc e;
  e.t = t;
  e.bits &= TAG_MASK;
  return (long)e.bits;
}
  
template <typename Value>
static inline int tag_of(Value v) {
  return (int)tag_of<long>(v);
}

} // end namespace
} // end namespace

#endif /*! _HEARTBEAT_TAGGED_H_ */
