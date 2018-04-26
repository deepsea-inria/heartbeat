#include <atomic>
#include <array>
#include <type_traits>
#include <assert.h>
#include <cstdlib>
#include <memory>

#include "forward.hpp"
#include "tagged.hpp"
#include "atomic.hpp"

#ifndef _HEARTBEAT_SCHED_INCOUNTER_H_
#define _HEARTBEAT_SCHED_INCOUNTER_H_

namespace heartbeat {
namespace sched {

namespace gsnzi {
  
/*---------------------------------------------------------------------*/
/* Growable, scalable non-zero indicator */
    
static constexpr
int cache_align_szb = 128;
  
static constexpr
int nb_children = 2;

// pre: sizeof(ty) <= align_szb
#define DECLARE_PADDED_FIELD(ty, name, align_szb)   \
  ty name;                                          \
  char _padding_for_ ## name[align_szb - sizeof(ty)]
    
template <int saturation_upper_bound>
class node {
private:
  
  static constexpr
  int one_half = -1;
  
  static constexpr
  int root_node_tag = 1;
  
  using contents_type = struct {
    int c; // counter value
    int v; // version number
  };
  
  DECLARE_PADDED_FIELD(std::atomic<contents_type>, X, cache_align_szb);
  
  DECLARE_PADDED_FIELD(node*, parent, cache_align_szb);
  
  static
  bool is_root_node(const node* n) {
    return tagged::tag_of(n) == root_node_tag;
  }
  
  template <class Item>
  static
  node* create_root_node(Item x) {
    return (node*)tagged::tag_with(x, root_node_tag);
  }
  
public:
  
  // post: constructor zeros out all fields, except parent
  node(node* _parent = nullptr) {
    parent = (_parent == nullptr) ? create_root_node(_parent) : _parent;
    contents_type init;
    init.c = 0;
    init.v = 0;
    X.store(init);
  }
  
  bool is_saturated() const {
    return X.load().v >= saturation_upper_bound;
  }
  
  bool is_nonzero() const {
    return X.load().c > 0;
  }
  
  void increment() {
    bool succ = false;
    int undo_arr = 0;
    while (! succ) {
      contents_type x = X.load();
      if (x.c >= 1) {
        contents_type orig = x;
        contents_type next = x;
        next.c++;
        next.v++;
        succ = atomic::compare_exchange(X, orig, next);
      }
      if (x.c == 0) {
        contents_type orig = x;
        contents_type next = x;
        next.c = one_half;
        next.v++;
        if (atomic::compare_exchange(X, orig, next)) {
          succ = true;
          x.c = one_half;
          x.v++;
        }
      }
      if (x.c == one_half) {
        if (! is_root_node(parent)) {
          parent->increment();
        }
        contents_type orig = x;
        contents_type next = x;
        next.c = 1;
        if (! atomic::compare_exchange(X, orig, next)) {
          undo_arr++;
        }
      }
    }
    if (is_root_node(parent)) {
      return;
    }
    while (undo_arr > 0) {
      parent->decrement();
      undo_arr--;
    }
  }
  
  bool decrement() {
    while (true) {
      contents_type x = X.load();
      assert(x.c >= 1);
      contents_type orig = x;
      contents_type next = x;
      next.c--;
      if (atomic::compare_exchange(X, orig, next)) {
        bool s = (x.c == 1);
        if (is_root_node(parent)) {
          return s;
        } else if (s) {
          return parent->decrement();
        } else {
          return false;
        }
      }
    }
  }
  
  template <class Item>
  static
  void set_root_annotation(node* n, Item x) {
    node* m = n;
    assert(! is_root_node(m));
    while (! is_root_node(m->parent)) {
      m = m->parent;
    }
    assert(is_root_node(m->parent));
    m->parent = create_root_node(x);
  }
  
  template <class Item>
  static
  Item get_root_annotation(node* n) {
    node* m = n;
    while (! is_root_node(m)) {
      m = m->parent;
    }
    assert(is_root_node(m));
    return (Item)tagged::pointer_of(m);
  }
  
};

#undef DECLARE_PADDED_FIELD

template <
  int max_height = 6,  // must be >= 0
  int saturation_upper_bound = (1<<(max_height-1)) // any constant fraction of 2^max_height
>
class tree {
public:
  
  using node_type = node<saturation_upper_bound>;
  
private:
  
  static constexpr
  int nb_leaves = 1 << max_height;
  
  static constexpr
  int heap_size = 2 * nb_leaves;
  
  static constexpr
  int loading_heap_tag = 1;
  
  static
  unsigned int hashu(unsigned int a) {
    a = (a+0x7ed55d16) + (a<<12);
    a = (a^0xc761c23c) ^ (a>>19);
    a = (a+0x165667b1) + (a<<5);
    a = (a+0xd3a2646c) ^ (a<<9);
    a = (a+0xfd7046c5) + (a<<3);
    a = (a^0xb55a4f09) ^ (a>>16);
    return a;
  }
  
  template <class Item>
  static
  unsigned int random_path_for(Item x) {
    union {
      Item x;
      long b;
    } bits;
    bits.x = x;
    return std::abs((int)hashu((unsigned int)bits.b));
  }
  
  node_type root;
  
  // if (heap.load() != nullptr), then we have a representation of a tree of height
  // max_height, using the array-based binary tree representation
  std::atomic<node_type*> heap;
  
  // only called once
  // pre: assume tagged::tag_of(heap.load()) == loading_heap_tag
  // post: heap.load() points to an array of pre-initialized SNZI nodes
  void create_heap() {
    assert(tagged::tag_of(heap.load()) == loading_heap_tag);
    size_t szb = heap_size * sizeof(node_type);
    node_type* h = (node_type*)malloc(szb);
    // cells at indices 0 and 1 are not used
    for (int i = 2; i < 4; i++) {
      new (&h[i]) node_type(&root);
    }
    for (int i = 4; i < heap_size; i++) {
      new (&h[i]) node_type(&h[i / 2]);
    }
    heap.store(h);
  }
  
public:
  
  tree() {
    assert(max_height > 1);
    heap.store(nullptr);
  }
  
  ~tree() {
    node_type* h = heap.load();
    assert(tagged::tag_of(h) != loading_heap_tag);
    if (h != nullptr) {
      free(heap);
    }
  }
  
  bool is_nonzero() const {
    return root.is_nonzero();
  }
  
  node_type* get_target_of_path(unsigned int path) {
    node_type* h = heap.load();
    if ((h != nullptr) && (tagged::tag_of(h) != loading_heap_tag)) {
      int i = nb_leaves + (path & (nb_leaves - 1));
      assert(i >= 2 && i < heap_size);
      return &h[i];
    } else if ((h == nullptr) && (root.is_saturated())) {
      node_type* orig = nullptr;
      node_type* next = tagged::tag_with<node_type>(nullptr, loading_heap_tag);
      if (atomic::compare_exchange(heap, orig, next)) {
        create_heap();
      }
    }
    return &root;
  }
  
  template <class Item>
  node_type* get_target_of_value(Item x) {
    return get_target_of_path(random_path_for(x));
  }
  
  template <class Item>
  void set_root_annotation(Item x) {
    node_type::set_root_annotation(&root, x);
  }
  
};
  
} // end namespace
  
/*---------------------------------------------------------------------*/
/* Incounter */

using incounter_tag_type = enum {
  incounter_tag_fetch_add,
  incounter_tag_tree
};

constexpr
int snzi_tree_height = 9;

using gsnzi_tree_type = gsnzi::tree<snzi_tree_height>;

using fetch_add_cell_type = struct {
  vertex* v;
  std::atomic<int> counter;
};
  
class incounter_handle {
public:

  using gsnzi_tree_node_type = typename gsnzi_tree_type::node_type;

  void* h = nullptr;
  
  incounter_handle() noexcept = default;

  void decrement() {
    assert(h != nullptr);
    vertex* v = nullptr;
    incounter_tag_type tag = static_cast<incounter_tag_type>(tagged::tag_of(h));
    switch (tag) {
      case incounter_tag_fetch_add: {
        auto fetch_add = tagged::value_of<fetch_add_cell_type*, void*>(h);
        auto n = --fetch_add->counter;
        if (n == 0) {
          v = fetch_add->v;
        }
        break;
      }
      case incounter_tag_tree: {
        auto n = tagged::value_of<gsnzi_tree_node_type*, void*>(h);
        if (n->decrement()) {
          v = gsnzi_tree_node_type::get_root_annotation<vertex*>(n);
        }
        break;
      }
    }
    if (v != nullptr) {
      schedule(v);
    }
  }

};

class incounter {
public:

  incounter_tag_type tag = incounter_tag_fetch_add;

  union variants_union {
    fetch_add_cell_type fetch_add;
    std::unique_ptr<gsnzi_tree_type> tree;
    variants_union() { }
    ~variants_union() { }
  } u;

  void construct_fetch_add(vertex* v) {
    new (&u.fetch_add) fetch_add_cell_type;
    u.fetch_add.v = v;
    u.fetch_add.counter.store(0);
  }

  void construct_tree(vertex* v) {
    u.tree.reset(new gsnzi_tree_type);
    u.tree->set_root_annotation(v);
  }

  incounter() { }

  incounter(vertex* v) {
    assert(tag == incounter_tag_fetch_add);
    construct_fetch_add(v);    
  }

  incounter(incounter_tag_type tag, vertex* v)
    : tag(tag) {
    switch (tag) {
      case incounter_tag_fetch_add: {
        construct_fetch_add(v);
        break;
      }
      case incounter_tag_tree: {
        new (&u.tree) std::unique_ptr<gsnzi_tree_type>;
        construct_tree(v);
        break;
      }
    }
  }

  ~incounter() {
    if (tag == incounter_tag_tree) {
      u.tree.~unique_ptr<gsnzi_tree_type>();
    }
  }

  template <class Item>
  incounter_handle increment(Item* x) {
    incounter_handle h;
    switch (tag) {
      case incounter_tag_fetch_add: {
        u.fetch_add.counter++;
        h.h = tagged::tag_with(&u.fetch_add, incounter_tag_fetch_add);
        break;
      }
      case incounter_tag_tree: {
        auto n = u.tree->get_target_of_value(x);
        n->increment();
        h.h = tagged::tag_with(n, incounter_tag_tree);
        break;
      }
    }
    return h;
  }

  static
  void decrement(incounter_handle h) {
    h.decrement();
  }

  void make_ready(vertex* v) {
    switch (tag) {
      case incounter_tag_fetch_add: {
        u.fetch_add.counter.store(0);
        break;
      }
      case incounter_tag_tree: {
        construct_tree(v);
        break;
      }
    }
  }

#ifndef NDEBUG
  bool is_nonzero() {
    switch (tag) {
      case incounter_tag_fetch_add: {
        return u.fetch_add.counter.load() != 0;
      }
      case incounter_tag_tree: {
        return u.tree->is_nonzero();
      }
    }
    assert(false);
    return true;
  }
#endif

};

} // end namespace
} // end namespace

#endif /*! _HEARTBEAT_SCHED_INCOUNTER_H_ */
