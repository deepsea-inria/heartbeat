
#include <functional>
#include <vector>
#include <map>
#include <algorithm>
#include <memory>
#include <array>

#include "vertex.hpp"
#include "cactus-plus.hpp"
#include "stats.hpp"
#include "scheduler.hpp"

#ifndef _HEARTBEAT_PCFG_H_
#define _HEARTBEAT_PCFG_H_

namespace heartbeat {
namespace edsl {
  
/*---------------------------------------------------------------------*/
/* Parallel Control-Flow Graph */
  
namespace pcfg {
  
namespace cactus = cactus_stack::plus;
  
using stack_type = cactus::stack_type;
  
class shared_activation_record;
  
class private_activation_record;
  
class interpreter;
    
using basic_block_label_type = int;

static constexpr
basic_block_label_type entry_block_label = 0;
  
static constexpr
basic_block_label_type exit_block_label = -1;

using basic_block_tag_type = enum {
  tag_unconditional_jump, tag_conditional_jump,
  tag_spawn_join, tag_spawn2_join,
  tag_tail,
  tag_join_plus, tag_spawn_minus,
  tag_spawn_plus, tag_join_minus,
  tag_none
};

template <class Shared_activation_record, class Private_activation_record>
class basic_block_type {
public:
  
  using sar_type = Shared_activation_record;
  using par_type = Private_activation_record;
  
  using unconditional_jump_code_type = std::function<void(sar_type&, par_type&)>;
  using conditional_jump_code_type = std::function<basic_block_label_type(sar_type&, par_type&)>;
  using procedure_call_code_type = std::function<stack_type(sar_type&, par_type&, cactus::parent_link_type, stack_type)>;
  using incounter_getter_code_type = std::function<sched::incounter**(sar_type&, par_type&)>;
  using future_getter_code_type = std::function<sched::future*(sar_type&, par_type&)>;
  
  basic_block_tag_type tag;
  
  union {
    struct {
      unconditional_jump_code_type code;
      basic_block_label_type next;
    } variant_unconditional_jump;
    struct {
      conditional_jump_code_type code;
    } variant_conditional_jump;
    struct {
      procedure_call_code_type code;
      basic_block_label_type next;
    } variant_spawn_join;
    struct {
      procedure_call_code_type code1;
      basic_block_label_type next; // must be the label of a spawn_join block
    } variant_spawn2_join;
    struct {
      procedure_call_code_type code;
      basic_block_label_type next;
    } variant_tail;
    struct {
      procedure_call_code_type code;
      incounter_getter_code_type getter;
      basic_block_label_type next;
    } variant_join_plus;
    struct {
      procedure_call_code_type code;
      incounter_getter_code_type getter;
      basic_block_label_type next;
    } variant_spawn_minus;
    struct {
      procedure_call_code_type code;
      future_getter_code_type getter;
      basic_block_label_type next;
    } variant_spawn_plus;
    struct {
      future_getter_code_type getter;
      basic_block_label_type next;
    } variant_join_minus;
  };
  
  basic_block_type() : tag(tag_none) { }
  
private:
  
  void copy_constructor(basic_block_type const& other) {
    tag = other.tag;
    switch (tag) {
      case tag_unconditional_jump: {
        new (&variant_unconditional_jump.code) unconditional_jump_code_type(other.variant_unconditional_jump.code);
        variant_unconditional_jump.next = other.variant_unconditional_jump.next;
        break;
      }
      case tag_conditional_jump: {
        new (&variant_conditional_jump.code) conditional_jump_code_type(other.variant_conditional_jump.code);
        break;
      }
      case tag_spawn_join: {
        new (&variant_spawn_join.code) procedure_call_code_type(other.variant_spawn_join.code);
        variant_spawn_join.next = other.variant_spawn_join.next;
        break;
      }
      case tag_spawn2_join: {
        new (&variant_spawn2_join.code1) procedure_call_code_type(other.variant_spawn2_join.code1);
        variant_spawn2_join.next = other.variant_spawn2_join.next;
        break;
      }
      case tag_tail: {
        new (&variant_tail.code) procedure_call_code_type(other.variant_tail.code);
        break;
      }
      case tag_join_plus: {
        new (&variant_join_plus.code) procedure_call_code_type(other.variant_join_plus.code);
        new (&variant_join_plus.getter) incounter_getter_code_type(other.variant_join_plus.getter);
        variant_join_plus.next = other.variant_join_plus.next;
        break;
      }
      case tag_spawn_minus: {
        new (&variant_spawn_minus.code) procedure_call_code_type(other.variant_spawn_minus.code);
        new (&variant_spawn_minus.getter) incounter_getter_code_type(other.variant_spawn_minus.getter);
        variant_spawn_minus.next = other.variant_spawn_minus.next;
        break;
      }
      case tag_spawn_plus: {
        new (&variant_spawn_plus.code) procedure_call_code_type(other.variant_spawn_plus.code);
        new (&variant_spawn_plus.getter) future_getter_code_type(other.variant_spawn_plus.getter);
        variant_spawn_plus.next = other.variant_spawn_plus.next;
        break;
      }
      case tag_join_minus: {
        new (&variant_join_minus.getter) future_getter_code_type(other.variant_join_minus.getter);
        variant_join_minus.next = other.variant_join_minus.next;
        break;
      }
      default:
        break;
    }
  }
  
public:
  
  basic_block_type(basic_block_type const& other) {
    copy_constructor(other);
  }
  
private:
  
  void move_constructor(basic_block_type&& other) {
    tag = other.tag;
    other.tag = tag_none;
    switch (tag) {
      case tag_unconditional_jump: {
        new (&variant_unconditional_jump.code) unconditional_jump_code_type;
        variant_unconditional_jump.code = std::move(other.variant_unconditional_jump.code);
        variant_unconditional_jump.next = std::move(other.variant_unconditional_jump.next);
        break;
      }
      case tag_conditional_jump: {
        new (&variant_conditional_jump.code) conditional_jump_code_type;
        variant_conditional_jump.code = std::move(other.variant_conditional_jump.code);
        break;
      }
      case tag_spawn_join: {
        new (&variant_spawn_join.code) procedure_call_code_type;
        variant_spawn_join.code = std::move(other.variant_spawn_join.code);
        variant_spawn_join.next = std::move(other.variant_spawn_join.next);
        break;
      }
      case tag_spawn2_join: {
        new (&variant_spawn2_join.code1) procedure_call_code_type;
        variant_spawn2_join.code1 = std::move(other.variant_spawn2_join.code1);
        variant_spawn2_join.next = std::move(other.variant_spawn2_join.next);
        break;
      }
      case tag_tail: {
        new (&variant_tail.code) procedure_call_code_type;
        variant_tail.code = std::move(other.variant_tail.code);
        break;
      }
      case tag_join_plus: {
        new (&variant_join_plus.code) procedure_call_code_type;
        variant_join_plus.code = std::move(other.variant_join_plus.code);
        new (&variant_join_plus.getter) incounter_getter_code_type;
        variant_join_plus.getter = std::move(other.variant_join_plus.getter);
        variant_join_plus.next = std::move(other.variant_join_plus.next);
        break;
      }
      case tag_spawn_minus: {
        new (&variant_spawn_minus.code) procedure_call_code_type;
        variant_spawn_minus.code = std::move(other.variant_spawn_minus.code);
        new (&variant_spawn_minus.getter) future_getter_code_type;
        variant_spawn_minus.getter = std::move(other.variant_spawn_minus.getter);
        variant_spawn_minus.next = std::move(other.variant_spawn_minus.next);
        break;
      }
      case tag_spawn_plus: {
        new (&variant_spawn_plus.code) procedure_call_code_type;
        variant_spawn_plus.code = std::move(other.variant_spawn_plus.code);
        new (&variant_spawn_plus.getter) future_getter_code_type;
        variant_spawn_plus.getter = std::move(other.variant_spawn_plus.getter);
        variant_spawn_plus.next = std::move(other.variant_spawn_plus.next);
        break;
      }
      case tag_join_minus: {
        new (&variant_join_minus.getter) future_getter_code_type;
        variant_join_minus.getter = std::move(other.variant_join_minus.getter);
        variant_join_minus.next = std::move(other.variant_join_minus.next);
        break;
      }
      default:
        break;
    }
  }
  
public:
  
  basic_block_type(basic_block_type&& other) {
    move_constructor(std::move(other));
  }
  
  ~basic_block_type() {
    switch (tag) {
      case tag_unconditional_jump: {
        variant_unconditional_jump.code.~unconditional_jump_code_type();
        break;
      }
      case tag_conditional_jump: {
        variant_conditional_jump.code.~conditional_jump_code_type();
        using blks = std::vector<basic_block_label_type>;
        break;
      }
      case tag_spawn_join: {
        variant_spawn_join.code.~procedure_call_code_type();
        break;
      }
      case tag_spawn2_join: {
        variant_spawn2_join.code1.~procedure_call_code_type();
        break;
      }
      case tag_tail: {
        variant_tail.code.~procedure_call_code_type();
        break;
      }
      case tag_join_plus: {
        variant_join_plus.getter.~incounter_getter_code_type();
        variant_join_plus.code.~procedure_call_code_type();
        break;
      }
      case tag_spawn_minus: {
        variant_spawn_minus.getter.~incounter_getter_code_type();
        variant_spawn_minus.code.~procedure_call_code_type();
        break;
      }
      case tag_spawn_plus: {
        variant_spawn_plus.getter.~future_getter_code_type();
        variant_spawn_plus.code.~procedure_call_code_type();
        break;
      }
      case tag_join_minus: {
        variant_join_minus.getter.~future_getter_code_type();
        break;
      }
      default:
        break;
    }
  }
  
  basic_block_type& operator=(basic_block_type const& other) {
    copy_constructor(other);
    return *this;
  }
  
  basic_block_type& operator=(basic_block_type&& other) {
    move_constructor(std::move(other));
    return *this;
  }
  
  static
  basic_block_type unconditional_jump(unconditional_jump_code_type code,
                                      basic_block_label_type next) {
    basic_block_type b;
    b.tag = tag_unconditional_jump;
    new (&b.variant_unconditional_jump.code) unconditional_jump_code_type(code);
    b.variant_unconditional_jump.next = next;
    return b;
  }
  
  static
  basic_block_type conditional_jump(conditional_jump_code_type code) {
    basic_block_type b;
    b.tag = tag_conditional_jump;
    new (&b.variant_conditional_jump.code) conditional_jump_code_type(code);
    return b;
  }
  
  static
  basic_block_type spawn_join(procedure_call_code_type code,
                              basic_block_label_type next) {
    basic_block_type b;
    b.tag = tag_spawn_join;
    new (&b.variant_spawn_join.code) procedure_call_code_type(code);
    b.variant_spawn_join.next = next;
    return b;
  }
  
  static
  basic_block_type spawn2_join(procedure_call_code_type code1,
                               basic_block_label_type next) {
    basic_block_type b;
    b.tag = tag_spawn2_join;
    new (&b.variant_spawn2_join.code1) procedure_call_code_type(code1);
    b.variant_spawn2_join.next = next;
    return b;
  }
  
  static
  basic_block_type tail(procedure_call_code_type code,
                        basic_block_label_type next) {
    basic_block_type b;
    b.t = tag_tail;
    new (&b.variant_tail.code) procedure_call_code_type(code);
    b.variant_tail.next = next;
    return b;
  }
  
  static
  basic_block_type join_plus(procedure_call_code_type code,
                             incounter_getter_code_type getter,
                             basic_block_label_type next) {
    basic_block_type b;
    b.tag = tag_join_plus;
    new (&b.variant_join_plus.code) procedure_call_code_type(code);
    new (&b.variant_join_plus.getter) incounter_getter_code_type(getter);
    b.variant_join_plus.next = next;
    return b;
  }
  
  static
  basic_block_type spawn_minus(procedure_call_code_type code,
                               incounter_getter_code_type getter,
                               basic_block_label_type next) {
    basic_block_type b;
    b.tag = tag_spawn_minus;
    new (&b.variant_spawn_minus.code) procedure_call_code_type(code);
    new (&b.variant_spawn_minus.getter) incounter_getter_code_type(getter);
    b.variant_spawn_minus.next = next;
    return b;
  }
  
  static
  basic_block_type spawn_plus(procedure_call_code_type code,
                              future_getter_code_type getter,
                              basic_block_label_type next) {
    basic_block_type b;
    b.tag = tag_spawn_plus;
    new (&b.variant_spawn_plus.code) procedure_call_code_type(code);
    new (&b.variant_spawn_plus.getter) future_getter_code_type(getter);
    b.variant_spawn_plus.next = next;
    return b;
  }
  
  static
  basic_block_type join_minus(future_getter_code_type getter,
                              basic_block_label_type next) {
    basic_block_type b;
    b.tag = tag_join_minus;
    new (&b.variant_join_minus.getter) future_getter_code_type(getter);
    b.variant_join_minus.next = next;
    return b;
  }
  
};
  
using trampoline_type = struct {
  basic_block_label_type pred;
  basic_block_label_type succ;
};
  
class children_record {
public:
  
  int nb = 0;
  
  using future_type = std::pair<sched::future, private_activation_record*>;
  
  std::vector<future_type> futures;
  
};
  
class parallel_loop_activation_record {
public:
  
  virtual
  int nb_strands() = 0;
  
  virtual
  void split(parallel_loop_activation_record*, int) = 0;
  
  virtual
  sched::vertex*& get_join() = 0;
  
  virtual
  std::unique_ptr<children_record>& get_children() = 0;
  
  virtual
  private_activation_record*& get_destination() = 0;

#ifndef NDEBUG
  virtual
  std::pair<int, int> loop_range() = 0;
#endif
  
};

using parallel_loop_id_type = int;

template <class Private_activation_record>
class parallel_loop_descriptor_type {
public:
  
  using par = Private_activation_record;
  
  using join_type = enum {
    join_trivial,
    join_binary_associative_combine,
    join_binary_associative_commutative_combine
  };
  
  join_type join = join_binary_associative_commutative_combine;
  
  trampoline_type entry;
  
  trampoline_type exit;
  
  // ids in parents vector are ordered from outermost to innermost
  std::vector<parallel_loop_id_type> parents;
  
  std::function<void(par&, parallel_loop_activation_record*)> initializer;
  
};
  
template <class Shared_activation_record>
using private_activation_record_of = typename Shared_activation_record::private_activation_record;
  
template <class Shared_activation_record>
class cfg_type {
public:
  
  using sar = Shared_activation_record;
  using par = private_activation_record_of<sar>;
  
  using basic_block_type = basic_block_type<sar, par>;
  
  cfg_type() { }
  
  cfg_type(size_t nb_blocks)
  : basic_blocks(nb_blocks) { }
  
  cfg_type(size_t nb_blocks, size_t nb_loops)
  : basic_blocks(nb_blocks), loop_descriptors(nb_loops), loop_of(nb_blocks) { }
  
  std::vector<basic_block_type> basic_blocks;
  
  std::vector<parallel_loop_descriptor_type<par>> loop_descriptors;
  
  // to map from basic_block_id to parallel_loop_id
  // loop_of[i] yields the id of the parallel loop that most immediately
  // encloses the basic block labeled i
  std::vector<parallel_loop_id_type> loop_of;
  
  size_t nb_basic_blocks() const {
    return basic_blocks.size();
  }
  
  size_t nb_loops() const {
    return loop_descriptors.size();
  }
  
  basic_block_type& operator[](size_t i) {
    return basic_blocks[i];
  }
  
};

} // end namespace
} // end namespace
} // end namespace

#endif /*! _HEARTBEAT_PCFG_H_ */
