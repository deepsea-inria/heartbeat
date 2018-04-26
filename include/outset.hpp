#include <deque>
#include <random>
#include <memory>

#include "incounter.hpp"
#include "cycles.hpp"
#include "aligned.hpp"
#include "perworker.hpp"

#ifndef _HEARTBEAT_SCHED_OUTSET_H_
#define _HEARTBEAT_SCHED_OUTSET_H_

namespace heartbeat {
namespace sched {
  
/*---------------------------------------------------------------------*/
/* Growable, scalable, concurrent bag, for storing the outgoing
 * edges on a vertex
 */

#if 0
  
static constexpr
int finished_tag = 1;
  
static constexpr
int creating_shortcuts_tag = 3;

template <class Item, int capacity, bool concurrent_inserts>
class block {
public:
  
  using cell_type = std::atomic<Item>;
  
private:
  
  cell_type start[capacity];
  
  std::atomic<cell_type*> head;
  
  static
  bool try_insert_item_at(cell_type& cell, Item x) {
    assert(concurrent_inserts);
    while (true) {
      Item y = cell.load();
      if (tagged::tag_of(y) == finished_tag) {
        return false;
      }
      assert(y.h == nullptr);
      Item orig = y;
      Item next = x;
      if (atomic::compare_exchange(cell, orig, next)) {
        return true;
      }
    }
  }
  
  static
  Item try_notify_item_at(cell_type& cell) {
    assert(concurrent_inserts);
    while (true) {
      Item y = cell.load();
      assert(tagged::tag_of(y) != finished_tag);
      if (y.h != nullptr) {
        return y;
      }
      Item orig;
      Item next = tagged::tag_with(Item(), finished_tag);
      if (atomic::compare_exchange(cell, orig, next)) {
        return Item();
      }
    }
  }
  
public:
  
  block() {
    assert(sizeof(Item) == sizeof(void*));
    head.store(start);
    if (concurrent_inserts) {
      for (int i = 0; i < capacity; i++) {
        start[i].store(Item());
      }
    }
  }
  
  bool is_full() const {
    return head.load() >= (start + capacity);
  }
  
  using try_insert_result_type = enum {
    succeeded_and_not_filled, succeeded_and_filled,
    failed_because_finish, failed_because_full
  };
  
  try_insert_result_type try_insert(Item x) {
    assert(x.h != nullptr);
    while (true) {
      cell_type* h = head.load();
      if (tagged::tag_of(h) == finished_tag) {
        return failed_because_finish;
      }
      if (h >= (start + capacity)) {
        return failed_because_full;
      }
      cell_type* orig = h;
      if (concurrent_inserts) {
        if (atomic::compare_exchange(head, orig, h + 1)) {
          if (! try_insert_item_at(*h, x)) {
            return failed_because_finish;
          }
          if (start + capacity == h + 1) {
            return succeeded_and_filled;
          } else {
            return succeeded_and_not_filled;
          }
        }
      } else {
        h->store(x);
        if (atomic::compare_exchange(head, orig, h + 1)) {
          if (start + capacity == h + 1) {
            return succeeded_and_filled;
          } else {
            return succeeded_and_not_filled;
          }
        }
      }
    }
  }
  
  std::pair<cell_type*, cell_type*> notify_init() {
    while (true) {
      cell_type* h = head.load();
      cell_type* orig = h;
      cell_type* next = tagged::tag_with<cell_type>(nullptr, finished_tag);
      if (atomic::compare_exchange(head, orig, next)) {
        assert(h <= start + capacity);
        return std::make_pair(start, h);
      }
    }
  }
  
  template <class Visit>
  static
  void notify_rng(cell_type* lo, cell_type* hi, const Visit& visit) {
    for (cell_type* it = lo; it != hi; it++) {
      if (concurrent_inserts) {
        Item x = try_notify_item_at(*it);
        if (x.h != nullptr) {
          visit(x);
        }
      } else {
        Item x = it->load();
        assert(x.h != nullptr);
        visit(x);
      }
    }
  }
  
  template <class Visit>
  void notify(const Visit& visit) {
    auto rng = notify_init();
    notify_rng(rng.first, rng.second, visit);
  }
  
};
  
template <class Item, int branching_factor, int block_capacity>
class node {
public:
  
  using block_type = block<Item, block_capacity, false>;
  
  block_type items;
  
  data::cache_aligned_fixed_capacity_array<std::atomic<node*>, branching_factor> children;
  
  node() {
    for (int i = 0; i < branching_factor; i++) {
      children[i].store(nullptr);
    }
  }
  
};

template <class Item, int branching_factor, int block_capacity>
class tree {
public:
  
  using node_type = node<Item, branching_factor, block_capacity>;
  
  std::atomic<node_type*> root;
  
  tree() {
    root.store(nullptr);
  }
  
  template <class Random_int>
  node_type* try_insert(const Random_int& random_int) {
    node_type* new_node = nullptr;
    std::atomic<node_type*>* current = &root;
    while (true) {
      node_type* target = current->load();
      if (target == nullptr) {
        if (new_node == nullptr) {
          new_node = new node_type;
        }
        assert(new_node != nullptr);
        node_type* orig = nullptr;
        if (atomic::compare_exchange(*current, orig, new_node)) {
          return new_node;
        }
        target = current->load();
      }
      if (tagged::tag_of(target) == finished_tag) {
        if (new_node != nullptr) {
          delete new_node;
          new_node = nullptr;
        }
        assert(new_node == nullptr);
        return nullptr;
      }
      int i = random_int(0, branching_factor);
      current = &(target->children[i]);
    }
  }
  
};
  
data::perworker::array<std::mt19937> outset_rngs;  // random-number generators
  
/*---------------------------------------------------------------------*/
/* Tree-based outset for scalable, high outdegree vertices */
  
class tree_outset {
private:
  
  using value_type = incounter_handle;
  
  static constexpr
  int branching_factor = 4;
  
  static constexpr
  int block_capacity = 4096;
  
public:
  
  static constexpr
  int max_nb_procs = 64;
  
  using tree_type = tree<value_type, branching_factor, block_capacity>;
  using node_type = typename tree_type::node_type;
  using item_cell_type = typename node_type::block_type::cell_type;
  using item_iterator = item_cell_type*;
  
private:
  
  using block_type = typename node_type::block_type;
  using shortcuts_type = data::cache_aligned_fixed_capacity_array<block_type*, max_nb_procs>;
  
  static constexpr
  int small_block_capacity = 16;
  
  using small_block_type = block<value_type, small_block_capacity, true>;
  
  small_block_type items;
  
  tree_type blocks;
  
  std::atomic<shortcuts_type*> shortcuts;
  
  template <class Random_int>
  void create_shortcuts(const Random_int& random_int) {
    while (true) {
      shortcuts_type* s = shortcuts.load();
      if (tagged::tag_of(s) == finished_tag) {
        return;
      }
      shortcuts_type* orig = s;
      shortcuts_type* next = tagged::tag_with(s, creating_shortcuts_tag);
      if (atomic::compare_exchange(shortcuts, orig, next)) {
        break;
      }
    }
    node_type* root = blocks.try_insert(random_int);
    assert(root != nullptr);
    shortcuts_type* new_shortcuts = new shortcuts_type;
    (*new_shortcuts)[0] = &(root->items);
    for (int i = 1; i < new_shortcuts->size(); i++) {
      node_type* n = blocks.try_insert(random_int);
      assert(n != nullptr);
      (*new_shortcuts)[i] = &(n->items);
    }
    shortcuts.store(new_shortcuts);
  }
  
  template <class Random_int>
  bool insert(value_type x, int my_id, const Random_int& random_int) {
    if (! items.is_full()) {
      auto status = items.try_insert(x);
      if (status == items.failed_because_finish) {
        return false;
      } else if (status == items.succeeded_and_not_filled) {
        return true;
      } else if (status == items.succeeded_and_filled) {
        create_shortcuts(random_int);
        return true;
      }
      assert(status == items.failed_because_full);
    }
    assert(items.is_full());
    while (true) {
      shortcuts_type* s = shortcuts.load();
      if (tagged::tag_of(s) == finished_tag) {
        return false;
      }
      if ( (s == nullptr) || (tagged::tag_of(s) == creating_shortcuts_tag) ) {
        cycles::spin_for(128);
        continue;
      }
      assert(s != nullptr);
      block_type* b = (*s)[my_id];
      auto status = b->try_insert(x);
      if (status == b->failed_because_finish) {
        return false;
      } else if (   (status == b->succeeded_and_filled)
                 || (status == b->succeeded_and_not_filled) ) {
        return true;
      }
      assert(status == b->failed_because_full);
      node_type* n = blocks.try_insert(random_int);
      if (n == nullptr) {
        return false;
      }
      (*s)[my_id] = &(n->items);
    }
  }
  
public:
  
  tree_outset() {
    shortcuts.store(nullptr);
  }
  
  ~tree_outset() {
    shortcuts_type* s = tagged::pointer_of(shortcuts.load());
    shortcuts.store(nullptr);
    if (s != nullptr) {
      delete s;
    }
    assert(false); // todo: free memory held by tree nodes
  }
  
  bool insert(value_type x) {
    int my_id = data::perworker::get_my_id();
    return insert(x, my_id, [&] (int lo, int hi) {
      std::uniform_int_distribution<int> distribution(lo, hi-1);
      return distribution(outset_rngs.mine());
    });
  }
  
  node_type* get_root() {
    return tagged::pointer_of(blocks.root.load());
  }
  
  template <class Visit>
  node_type* notify_init(const Visit& visit) {
    while (true) {
      shortcuts_type* s = shortcuts.load();
      assert(tagged::tag_of(s) != finished_tag);
      if (tagged::tag_of(s) == creating_shortcuts_tag) {
        cycles::spin_for(128);
        continue;
      }
      shortcuts_type* orig = s;
      shortcuts_type* next = tagged::tag_with(s, finished_tag);
      if (atomic::compare_exchange(shortcuts, orig, next)) {
        break;
      }
    }
    items.notify(visit);
    while (true) {
      node_type* n = blocks.root.load();
      node_type* orig = n;
      node_type* next = tagged::tag_with(n, finished_tag);
      if (atomic::compare_exchange(blocks.root, orig, next)) {
        return n;
      }
    }
  }
  
  template <class Visit, class Deque>
  static
  void notify_nb(const std::size_t nb, item_iterator& lo,
                        item_iterator& hi, Deque& todo,
                        const Visit& visit) {
    int k = 0;
    while ( (k < nb) && ((! todo.empty()) || ((hi - lo) > 0)) ) {
      if ((hi - lo) > 0) {
        item_iterator lo_next = std::min(hi, lo + (nb - k));
        block_type::notify_rng(lo, lo_next, visit);
        k += lo_next - lo;
        lo = lo_next;
      } else {
        assert(! todo.empty());
        node_type* current = todo.back();
        todo.pop_back();
        for (int i = 0; i < branching_factor; i++) {
          while (true) {
            node_type* child = current->children[i].load();
            assert(tagged::tag_of(child) == 0);
            node_type* orig = child;
            node_type* next = tagged::tag_with(child, finished_tag);
            if (atomic::compare_exchange(current->children[i], orig, next)) {
              if (child != nullptr) {
                todo.push_back(child);
              }
              break;
            }
          }
        }
        auto rng = current->items.notify_init();
        lo = rng.first;
        hi = rng.second;
      }
    }
  }
  
  static
  void deallocate_nb(const int nb, std::deque<node_type*>& todo) {
    int k = 0;
    while ((k < nb) && (! todo.empty())) {
      node_type* current = todo.back();
      todo.pop_back();
      for (int i = 0; i < branching_factor; i++) {
        node_type* child = tagged::pointer_of(current->children[i].load());
        if (child != nullptr) {
          todo.push_back(child);
        }
      }
      delete current;
      k++;
    }
  }
  
};

#endif
  
/*---------------------------------------------------------------------*/
/* Chain-based outset for fast, low outdegree vertices */
  
class chain_outset {
public:
  
  using value_type = incounter_handle;
  using node_type = int;
  using item_iterator = int*;
  
private:
  
  using concurrent_list_type = struct concurrent_list_struct {
    incounter_handle h;
    struct concurrent_list_struct* next;
  };
  
  std::atomic<concurrent_list_type*> head;
  
  static constexpr
  int finished_tag = 1;
  
public:
  
  chain_outset()
  : head(nullptr) { }
  
  node_type* get_root() {
    return nullptr;
  }
  
  bool insert(value_type x) {
    bool result = true;
    concurrent_list_type* cell = new concurrent_list_type;
    cell->h = x;
    while (true) {
      concurrent_list_type* orig = head.load();
      if (tagged::tag_of(orig) == finished_tag) {
        result = false;
        delete cell;
        break;
      } else {
        cell->next = orig;
        if (atomic::compare_exchange(head, orig, cell)) {
          break;
        }
      }
    }
    return result;
  }
  
  template <class Visit>
  void notify(const Visit& visit) {
    concurrent_list_type* todo = nullptr;
    while (true) {
      concurrent_list_type* orig = head.load();
      concurrent_list_type* v = (concurrent_list_type*)nullptr;
      concurrent_list_type* next = tagged::tag_with(v, finished_tag);
      if (atomic::compare_exchange(head, orig, next)) {
        todo = orig;
        break;
      }
    }
    while (todo != nullptr) {
      visit(todo->h);
      concurrent_list_type* next = todo->next;
      delete todo;
      todo = next;
    }
  }
    
};
  
/*---------------------------------------------------------------------*/
/* Parallel futures */

using future_tag_type = enum {
  future_tag_chain,
  future_tag_tree
};

class future {
public:

  future_tag_type tag = future_tag_chain;

  struct {
    std::shared_ptr<chain_outset> chain;
    //    std::shared_ptr<tree_outset> tree;
  } u;

  future() { }

  future(chain_outset* p) {
    tag = future_tag_chain;
    u.chain.reset(p);
  }
  /*
  future(tree_outset* p) {
    tag = future_tag_tree;
    u.tree.reset(p);
    } */

  explicit operator bool() const noexcept {
    bool b = false;
    switch (tag) {
      case future_tag_chain: {
        if (u.chain) {
          b = true;
        }
        break;
      }
      case future_tag_tree: {
        /*        if (u.tree) {
          b = true;
          } */
        break;
      }
    }
    return b;
  }
  
  bool insert(incounter_handle h) {
    bool b;
    switch (tag) {
      case future_tag_chain: {
        b = u.chain->insert(h);
        break;
      }
      case future_tag_tree: {
        //        b = u.tree->insert(h);
        break;
      }
    }
    return b;
  }

  template <class Visit>
  void notify(const Visit& visit) {
    switch (tag) {
      case future_tag_chain: {
        u.chain->notify(visit);
        break;
      }
      case future_tag_tree: { /*
        using outset_tree_node_type = typename tree_outset::node_type;
        using item_iterator = typename tree_outset::item_iterator;
        auto n = u.tree->notify_init(visit);        
        std::deque<outset_tree_node_type*> todo;
        todo.push_back(n);
        item_iterator lo = nullptr;
        item_iterator hi = nullptr;
        auto is_finished = [&] {
          return hi == lo && todo.empty();
        };
        while (! is_finished()) {
          u.tree->notify_nb(128, lo, hi, todo, [&] (incounter_handle h) {
            incounter::decrement(h);
          });
          } */
        break;
      }
    }
  }
  
  void* get_pointer() {
    void* p;
    switch (tag) {
      case future_tag_chain: {
        p = u.chain.get();
        break;
      }
      case future_tag_tree: {
        //        p = u.tree.get();
        break;
      }
    }
    return p;
  }

};


/*---------------------------------------------------------------------*/
/* Outset */

using outset_tag_type = enum {
  outset_tag_unary,
  outset_tag_chain,
  outset_tag_future
};

class outset {
public:

  using value_type = incounter_handle;

  outset_tag_type tag = outset_tag_chain;

  union variants_union {
    value_type unary;
    chain_outset chain;
    future fut;
    variants_union() { }
    ~variants_union() { }
  } u;

  outset() {
    new (&u.chain) chain_outset;
  }

  ~outset() {
    switch (tag) {
      case outset_tag_unary: {
        new (&u.unary) value_type;
        break;
      }
      case outset_tag_chain: {
        u.chain.~chain_outset();
        break;
      }
      case outset_tag_future: {
        u.fut.~future();
        break;
      }
    }
  }
  
  bool insert(value_type x) {
    bool b = true;
    switch (tag) {
      case outset_tag_unary: {
        assert(u.unary.h == nullptr);
        u.unary = x;
        break;
      }
      case outset_tag_chain: {
        b = u.chain.insert(x);
        break;
      }
      case outset_tag_future: {
        b = u.fut.insert(x);
        break;
      }
    }
    return b;
  }

  template <class Visit>
  void notify(const Visit& visit) {
    switch (tag) {
      case outset_tag_unary: {
        incounter::decrement(u.unary);
        u.unary = incounter_handle();
        break;
      }
      case outset_tag_chain: {
        u.chain.notify(visit);
        break;
      }
      case outset_tag_future: {
        u.fut.notify(visit);
        break;
      }
    }
  }

  void make_unary() {
    assert(tag == outset_tag_chain);
    tag = outset_tag_unary;
    u.unary = incounter_handle();
  }

  future make_chain_future() {
    assert(tag == outset_tag_chain);
    tag = outset_tag_future;
    new (&u.fut) future(new chain_outset);
    return u.fut;
  }

  future make_tree_future() {
    assert(false);
    /*
    assert(tag == outset_tag_chain);
    tag = outset_tag_future;
    new (&u.fut) future(new tree_outset); */
    return u.fut;
  }
  
};
  
} // end namespace
} // end namespace

#endif /*! _HEARTBEAT_SCHED_OUTSET_H_ */
