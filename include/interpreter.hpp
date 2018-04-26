
#include <functional>
#include <vector>
#include <map>
#include <algorithm>
#include <memory>
#include <array>

#include "vertex.hpp"
#include "cactus-plus.hpp"
#include "stats.hpp"
#include "logging.hpp"
#include "scheduler.hpp"
#include "pcfg.hpp"
#include "grain.hpp"

#ifndef _HEARTBEAT_INTERPRETER_H_
#define _HEARTBEAT_INTERPRETER_H_

namespace heartbeat {
namespace edsl {
namespace pcfg {

/*---------------------------------------------------------------------*/
/* Parallel control-flow graph interpreter */
  
class shared_activation_record {
public:
  
  virtual
  std::pair<stack_type, fuel::check_type> run(stack_type) const = 0;
  
  virtual
  void promote_mark(interpreter*, private_activation_record*) = 0;
  
  virtual
  sched::future get_dependency_of_join_minus(stack_type stack) = 0;
  
  virtual
  const char* get_name() {
    return "no name";
  }
  
};

using vertex_split_type = sched::vertex_split_type;
  
sched::vertex* dummy_join = nullptr;
  
std::unique_ptr<children_record> dummy_children = nullptr;
  
private_activation_record* dummy_destination = nullptr;
  
class private_activation_record {
public:
  
  trampoline_type trampoline = { .pred = entry_block_label, .succ = entry_block_label };
  
  virtual
  int nb_strands() {
    return 1;
  }
  
  virtual
  vertex_split_type split(interpreter*, int) {
    assert(false); // impossible
    return sched::make_vertex_split(nullptr, nullptr);
  }
  
  virtual
  sched::vertex*& get_join(parallel_loop_id_type) {
    assert(false); // impossible
    return dummy_join;
  }
  
  virtual
  parallel_loop_activation_record* loop_activation_record_of(parallel_loop_id_type id) {
    assert(false); // impossible
    return nullptr;
  }

#ifndef NDEBUG
  virtual
  std::vector<std::pair<int,int>> loop_activation_records() {
    return std::vector<std::pair<int,int>>();
  }
#endif
  
};
  
template <class Shared_activation_record>
Shared_activation_record* get_shared_frame_pointer(char* _ar) {
  return (Shared_activation_record*)(_ar + sizeof(size_t));
}
  
template <class Private_activation_record>
Private_activation_record* get_private_frame_pointer(char* _ar) {
  size_t shared_szb = *((size_t*)_ar);
  return (Private_activation_record*)(_ar + sizeof(size_t) + shared_szb);
}
  
bool is_splittable(char* _ar) {
  return (get_private_frame_pointer<private_activation_record>(_ar)->nb_strands() >= 2);
}
  
template <class Shared_activation_record, class Private_activation_record>
stack_type create_stack(Shared_activation_record* sar, Private_activation_record* par) {
  static constexpr size_t shared_szb = sizeof(Shared_activation_record);
  static constexpr size_t frame_szb = sizeof(size_t) + shared_szb + sizeof(Private_activation_record);
  return cactus::create_stack<frame_szb>(cactus::Parent_link_sync, [&] (char* _ar) {
    new ((size_t*)_ar) size_t(sizeof(Shared_activation_record*));
    new (get_shared_frame_pointer<Shared_activation_record*>(_ar)) Shared_activation_record*(sar);
    new (get_private_frame_pointer<Private_activation_record>(_ar)) Private_activation_record(*par);
  }, [&] (char* _ar) {
    return is_splittable(_ar);
  });
}

bool empty_stack(stack_type s) {
  return cactus::empty(s);
}

bool empty_mark_stack(stack_type s) {
  return cactus::empty_mark(s);
}
  
template <class Shared_activation_record>
Shared_activation_record* resolve_shared_frame_pointer(char* _ar, cactus::shared_frame_type sft) {
  Shared_activation_record* sar;
  if (sft == cactus::Shared_frame_indirect) {
    Shared_activation_record** p = get_shared_frame_pointer<Shared_activation_record*>(_ar);
    assert(p != nullptr);
    sar = *p;
  } else {
    sar = get_shared_frame_pointer<Shared_activation_record>(_ar);
  }
  return sar;
}

template <class Shared_activation_record>
Shared_activation_record& peek_newest_shared_frame(stack_type s) {
  Shared_activation_record* sar;
  cactus::peek_back(s, [&] (cactus::shared_frame_type sft, cactus::call_link_type, char* _ar) {
    sar = resolve_shared_frame_pointer<Shared_activation_record>(_ar, sft);
  });
  return *sar;
}

template <class Private_activation_record>
Private_activation_record& peek_newest_private_frame(stack_type s) {
  Private_activation_record* par;
  cactus::peek_back(s, [&] (cactus::shared_frame_type, cactus::call_link_type, char* _ar) {
    par = get_private_frame_pointer<Private_activation_record>(_ar);
  });
  return *par;
}
  
using peek_mark_result_tag = enum {
  Peek_mark_none, Peek_mark_fork, Peek_mark_loop_split
};
  
using peek_mark_result_type = struct {
  peek_mark_result_tag tag;
  shared_activation_record* sar;
  private_activation_record* par;
};

peek_mark_result_type peek_mark(stack_type s) {
  peek_mark_result_type r;
  if (empty_mark_stack(s)) {
    r.tag = Peek_mark_none;
    return r;
  }
  cactus::peek_mark(s, [&] (cactus::shared_frame_type sft,
                            cactus::call_link_type clt,
                            char* _ar,
                            cactus::shared_frame_type pred_sft,
                            char* _pred_ar) {
    r.par = get_private_frame_pointer<private_activation_record>(_ar);
    if ((clt == cactus::Call_link_async) && (_pred_ar != nullptr)) {
      r.tag = Peek_mark_fork;
      r.sar = resolve_shared_frame_pointer<shared_activation_record>(_pred_ar, pred_sft);
      r.par = get_private_frame_pointer<private_activation_record>(_pred_ar);
    } else if (r.par->nb_strands() >= 2) {
      r.tag = Peek_mark_loop_split;
      r.sar = resolve_shared_frame_pointer<shared_activation_record>(_ar, sft);
    } else {
      r.tag = Peek_mark_none;
    }
  });
  return r;
}

template <class Shared_activation_record>
Shared_activation_record& peek_marked_shared_frame(stack_type s) {
  Shared_activation_record* sar;
  cactus::peek_mark(s, [&] (cactus::shared_frame_type sft,
                            cactus::call_link_type,
                            char* _ar,
                            cactus::shared_frame_type,
                            char*) {
    sar = resolve_shared_frame_pointer<Shared_activation_record>(_ar, sft);
  });
  return *sar;
}

template <class Private_activation_record>
Private_activation_record& peek_marked_private_frame(stack_type s) {
  Private_activation_record* par;
  cactus::peek_mark(s, [&] (cactus::shared_frame_type sft,
                            cactus::call_link_type,
                            char* _ar,
                            cactus::shared_frame_type,
                            char*) {
    par = get_private_frame_pointer<Private_activation_record>(_ar);
  });
  return *par;
}

template <class Shared_activation_record, class ...Args>
stack_type push_call(stack_type s, cactus::parent_link_type plt, Args... args) {
  // The line below is commented out to avoid a GCC bug whereby GCC segfaults.
  // using private_activation_record = private_activation_record_of<Shared_activation_record>;
  using private_activation_record = typename Shared_activation_record::private_activation_record;
  static constexpr size_t shared_szb = sizeof(Shared_activation_record);
  static constexpr size_t frame_szb = sizeof(size_t) + shared_szb + sizeof(private_activation_record);
  stack_type t = cactus::push_back<frame_szb>(s, plt, [&] (char* _ar)  {
    new ((size_t*)_ar) size_t(shared_szb);
    new (get_shared_frame_pointer<Shared_activation_record>(_ar)) Shared_activation_record(args...);
    new (get_private_frame_pointer<private_activation_record>(_ar)) private_activation_record;
  }, [&] (char* _ar) {
    return is_splittable(_ar);
  });
  return t;
}
  
template <class Shared_activation_record>
stack_type pop_call(stack_type s) {
  using private_activation_record = typename Shared_activation_record::private_activation_record;
  return cactus::pop_back(s, [&] (char* _ar, cactus::shared_frame_type sft) {
    if (sft == cactus::Shared_frame_direct) {
      get_shared_frame_pointer<Shared_activation_record>(_ar)->~Shared_activation_record();
    } else {
      *get_shared_frame_pointer<Shared_activation_record*>(_ar) = nullptr;
    }
    get_private_frame_pointer<private_activation_record>(_ar)->~private_activation_record();
  });
}
  
std::pair<stack_type, stack_type> split_stack(stack_type s) {
  return cactus::split_mark(s, [&] (char* _ar) {
    return is_splittable(_ar);
  });
}

std::pair<stack_type, stack_type> fork_stack(stack_type s) {
  return cactus::fork_mark(s, [&] (char* _ar) {
    return is_splittable(_ar);
  });
}

#ifdef DEBUG_HEARTBEAT_STACK
  
using frame_summary_type = struct {
  cactus_stack::plus::frame_header_type* fp;
  cactus_stack::plus::call_link_type link_type;
  cactus_stack::plus::loop_link_type loop_link_type;
  std::vector<std::pair<int,int>> loops;
};

bool operator==(std::vector<std::pair<int,int>> loops1,
                std::vector<std::pair<int,int>> loops2) {
  auto n = loops1.size();
  if (n != loops2.size()) {
    return false;
  }
  for (auto i = 0; i < n; i++) {
    if (loops1[i].first != loops2[i].first) {
      return false;
    }
    if (loops1[i].second != loops2[i].second) {
      return false;
    }
  }
  return true;
}

bool operator!=(std::vector<std::pair<int,int>> loops1,
                std::vector<std::pair<int,int>> loops2) {
  return ! (loops1 == loops2);
}

bool operator==(frame_summary_type fs1, frame_summary_type fs2) {
  if (fs1.fp != fs2.fp) {
    return false;
  }
  if (fs1.link_type != fs2.link_type) {
    return false;
  }
  if (fs1.loops != fs2.loops) {
    return false;
  }
  if (fs1.loop_link_type != fs2.loop_link_type) {
    return false;
  }
  return true;
}

bool operator!=(frame_summary_type fs1, frame_summary_type fs2) {
  return ! (fs1 == fs2);
}

void print_frame_summary(frame_summary_type fs) {
  printf("fs.fp=%p \t", fs.fp);
  if (fs.link_type == cactus_stack::plus::Call_link_async) {
    printf("fs.link_type = async \t");
  } else {
    printf("fs.link_type = sync \t");
  }
  if (fs.loop_link_type == cactus_stack::plus::Loop_link_child) {
    printf("fs.loop_link_type = child \t");
  } else {
    printf("fs.loop_link_type = none \t");
  }  
  printf("|fs.loops|=%lld", fs.loops.size());
}

void print_frame_summaries(std::vector<frame_summary_type> vfs) {
  for (auto fs : vfs) {
    print_frame_summary(fs);
    printf("\n");
  }
}

frame_summary_type summarize_frame(cactus_stack::plus::frame_header_type* fp) {
  frame_summary_type result;
  result.fp = fp;
  result.link_type = fp->ext.clt;
  result.loop_link_type = fp->ext.llt;
  char* ar = frame_data(fp);
  private_activation_record* par = get_private_frame_pointer<private_activation_record>(ar);
  result.loops = par->loop_activation_records();
  return result;    
}

bool has_parallel_loop(frame_summary_type s) {
  for (auto r : s.loops) {
    if ((r.second - r.first) >= 2) {
      return true;
    }
  }
  return false;
}

bool is_mark(frame_summary_type s) {
  if (s.link_type == cactus_stack::plus::Call_link_async) {
    return true;
  }
  return has_parallel_loop(s);
}

std::vector<frame_summary_type> summarize_stack(stack_type s) {
  std::vector<frame_summary_type> result;
  auto start = s.begin();
  for (auto it = start; it != s.end(); it++) {
    auto& fh = *it;
    frame_summary_type fs = summarize_frame(&fh);
    result.push_back(fs);
  }
  return result;
}

std::vector<frame_summary_type> summarize_stack_marks(stack_type s) {
  std::vector<frame_summary_type> result;
  auto start = s.begin();
  for (auto it = start; it != s.end(); it++) {
    auto& fh = *it;
    frame_summary_type fs = summarize_frame(&fh);
    if (! is_mark(fs)) {
      continue;
    }
    result.push_back(fs);
  }
  return result;
}

std::vector<frame_summary_type> summarize_mark_stack(stack_type s, bool& broke_invar) {
  std::vector<frame_summary_type> result;
  auto start = s.begin_mark();
  for (auto it = start; it != s.end_mark(); it++) {
    auto& fh = *it;
    frame_summary_type fs = summarize_frame(&fh);
    if (! is_mark(fs)) {
      if (fh.ext.llt != cactus_stack::plus::Loop_link_child) {
        broke_invar = true;
      } else {
        continue;
      }
    }
    result.push_back(fs);
  }
  return result;
}

void check_stack(stack_type s) {
  std::vector<frame_summary_type> vfs = summarize_stack(s);
  std::vector<frame_summary_type> vfs_marks = summarize_stack_marks(s);
  bool broke_invar = false;
  std::vector<frame_summary_type> mark_vfs = summarize_mark_stack(s, broke_invar);
  if (broke_invar || (vfs_marks != mark_vfs)) {    
    printf("-------------------\n");
    printf("stack (size=%lld)\n", vfs.size());
    print_frame_summaries(vfs);
    printf("\n");
    printf("stack marks (size=%lld)\n", vfs_marks.size());
    print_frame_summaries(vfs_marks);
    printf("\n");
    printf("mark stack (size=%lld)\n", mark_vfs.size());
    print_frame_summaries(mark_vfs);
    printf("\n");
    printf("-------------------\n");
    printf("\n");
    assert(false);
  }
}

#endif
  
bool never_promote = false;
  
class interpreter : public sched::vertex {
public:

  stack_type stack;
  
  interpreter() {
    stack = cactus::create_stack();
  }
  
  interpreter(stack_type stack)
  : stack(stack) { }
  
  int nb_strands() {
    if (empty_stack(stack)) {
      return 0;
    } else if (empty_mark_stack(stack)) {
      return 1;
    } else {
      return peek_marked_private_frame<private_activation_record>(stack).nb_strands();
    }
  }
  
  fuel::check_type run() {
    fuel::check_type f = fuel::check_no_promote;
    {
      stack_type s = stack;
      while ((! empty_stack(s)) && (f == fuel::check_no_promote)) {
        auto r = peek_newest_shared_frame<shared_activation_record>(s).run(s);
        s = r.first;
        f = r.second;
      }
      stack = cactus::update_mark_stack(s, [&] (char* _ar) {
        return pcfg::is_splittable(_ar);
      });
    }
#ifdef DEBUG_HEARTBEAT_STACK
    check_stack(stack);
#endif
    if (nb_strands() == 0) {
      assert(! is_suspended);
      assert(f != fuel::check_suspend);
      return fuel::check_no_promote;
    }
    if (f == fuel::check_suspend) {
      is_suspended = true;
      f = fuel::check_yes_promote;
    }
    if (f == fuel::check_no_promote) {
      return f;
    }
    if (never_promote) {
      schedule(this);
      return fuel::check_no_promote;
    }
    assert(f == fuel::check_yes_promote);
    auto r = peek_mark(stack);
    switch (r.tag) {
      case Peek_mark_none: {
        if (is_suspended) {
          is_suspended = false;
          auto& par = peek_newest_shared_frame<shared_activation_record>(stack);
          auto dep = par.get_dependency_of_join_minus(stack);
          sched::new_edge(dep, this);
        } else {
          schedule(this);
        }
        break;
      }
      case Peek_mark_fork: {
        r.sar->promote_mark(this, r.par);
        break;
      }
      case Peek_mark_loop_split: {
        auto r = split(nb_strands() / 2);
        schedule(r.v2);
        schedule(r.v1);
        if (r.v0 != nullptr) {
          schedule(r.v0);
        }
        stats::on_promotion();
        break;
      }
      default: {
        assert(false);
        break;
      }
    }
    return f;
  }
  
  vertex_split_type split(int nb) {
    return peek_marked_private_frame<private_activation_record>(stack).split(this, nb);
  }
  
};

template <class Shared_activation_record>
std::pair<stack_type, fuel::check_type> step(cfg_type<Shared_activation_record>& cfg, stack_type stack) {
  using private_activation_record = private_activation_record_of<Shared_activation_record>;
  assert(! empty_stack(stack));
  fuel::check_type f = fuel::check_no_promote;
  auto& sar = peek_newest_shared_frame<Shared_activation_record>(stack);
  auto& par = peek_newest_private_frame<private_activation_record>(stack);
  basic_block_label_type pred = par.trampoline.succ;
  basic_block_label_type succ;
  if (pred == exit_block_label) {
    stack = pop_call<Shared_activation_record>(stack);
#ifdef DEBUG_HEARTBEAT_STACK
    check_stack(stack);
#endif
    return std::make_pair(stack, f);
  }
  auto start_time = cycles::now();
  assert(pred >= 0 && pred < cfg.nb_basic_blocks());
  bool possibly_updated_parallel_loop_range = false;
  auto& block = cfg.basic_blocks[pred];
  switch (block.tag) {
    case tag_unconditional_jump: {
      block.variant_unconditional_jump.code(sar, par);
      succ = block.variant_unconditional_jump.next;
      possibly_updated_parallel_loop_range = true;
      break;
    }
    case tag_conditional_jump: {
      succ = block.variant_conditional_jump.code(sar, par);
      possibly_updated_parallel_loop_range = true;      
      break;
    }
    case tag_spawn_join: {
      stack = block.variant_spawn_join.code(sar, par, cactus::Parent_link_sync, stack);
      succ = block.variant_spawn_join.next;
      break;
    }
    case tag_spawn2_join: {
      stack = block.variant_spawn2_join.code1(sar, par, cactus::Parent_link_async, stack);
      succ = block.variant_spawn2_join.next;
      break;
    }
    case tag_tail: {
      stack = pop_call<Shared_activation_record>(stack);
      stack = block.variant_tail.code(sar, par, cactus::Parent_link_sync, stack);
      succ = block.variant_tail.next;
      break;
    }
    case tag_join_plus: {
      *block.variant_join_plus.getter(sar, par) = nullptr;
      stack = block.variant_join_plus.code(sar, par, cactus::Parent_link_async, stack);
      succ = block.variant_join_plus.next;
      break;
    }
    case tag_spawn_minus: {
      stack = block.variant_spawn_minus.code(sar, par, cactus::Parent_link_async, stack);
      succ = block.variant_spawn_minus.next;
      break;
    }
    case tag_spawn_plus: {
      *block.variant_spawn_plus.getter(sar, par) = sched::future();
      stack = block.variant_spawn_plus.code(sar, par, cactus::Parent_link_async, stack);
      succ = block.variant_spawn_plus.next;
      break;
    }
    case tag_join_minus: {
      auto future = *block.variant_join_minus.getter(sar, par);
      if (future) {
        f = fuel::check_suspend;
      }
      succ = block.variant_join_minus.next;
      break;
    }
    default: {
      assert(false);
    }
  }
  par.trampoline.pred = pred;
  par.trampoline.succ = succ;
  auto end_time = cycles::now();
  f = (f == fuel::check_suspend) ? f : fuel::check(end_time);
  auto elapsed = cycles::diff(start_time, end_time);
  grain::callback(elapsed);
  if (possibly_updated_parallel_loop_range) {
    stack = cactus::update_mark_stack_just_for_loops(stack, [&] (char* _ar) {
              return pcfg::is_splittable(_ar);
            });
#ifdef DEBUG_HEARTBEAT_STACK
    check_stack(stack);
#endif
  }
#ifdef DEBUG_HEARTBEAT_STACK
  check_stack(stack);
#endif
  return std::make_pair(stack, f);
}

template <class Shared_activation_record, class Private_activation_record>
void promote_mark(cfg_type<Shared_activation_record>& cfg, interpreter* interp,
                  Shared_activation_record* sar, Private_activation_record* par) {
  basic_block_label_type pred = par->trampoline.pred;
  assert(pred != exit_block_label);
  assert(pred >= 0 && pred < cfg.nb_basic_blocks());
  auto& block = cfg.basic_blocks[pred];
  switch (block.tag) {
    case tag_spawn2_join: {
      interpreter* join = interp;
      auto stacks = fork_mark(interp->stack, [&] (char* _ar) {
        return is_splittable(_ar);
      });
      join->stack = stacks.first;
      interpreter* branch1 = new interpreter(stacks.second);
      interpreter* branch2 = new interpreter;
      branch1->get_outset()->make_unary();
      branch2->get_outset()->make_unary();
      basic_block_label_type pred = par->trampoline.succ;
      auto& spawn_join_block = cfg.basic_blocks[pred];
      assert(spawn_join_block.tag == tag_spawn_join);
      par->trampoline.pred = pred;
      par->trampoline.succ = spawn_join_block.variant_spawn_join.next;
      branch2->stack = spawn_join_block.variant_spawn_join.code(*sar, *par, cactus::Parent_link_sync, branch2->stack);
#ifdef DEBUG_HEARTBEAT_STACK
    check_stack(branch1->stack);
    check_stack(branch2->stack);
#endif
      sched::new_edge(branch2, join);
      sched::new_edge(branch1, join);
      release(branch2);
      release(branch1);
      logging::push_promote_spawn2_join(sar->get_name());
      break;
    }
    case tag_spawn_minus: {
      interpreter* continuation = interp;
      auto stacks = fork_mark(interp->stack, [&] (char* _ar) {
        return is_splittable(_ar);
      });
      continuation->stack = stacks.first;
      interpreter* branch = new interpreter(stacks.second);
      branch->get_outset()->make_unary();
      sched::incounter* incounter = *block.variant_spawn_minus.getter(*sar, *par);
      assert(incounter != nullptr);
      sched::new_edge(branch, incounter);
      schedule(continuation);
      release(branch);
      logging::push_promote_spawn_minus(sar->get_name());
      break;
    }
    case tag_spawn_plus: {
      interpreter* continuation = interp;
      auto stacks = fork_mark(interp->stack, [&] (char* _ar) {
        return is_splittable(_ar);
      });
      continuation->stack = stacks.first;
      interpreter* branch = new interpreter(stacks.second);
      auto branch_out = branch->get_outset();
      auto future = branch_out->make_chain_future();
      assert(! *block.variant_spawn_plus.getter(*sar, *par));
      *block.variant_spawn_plus.getter(*sar, *par) = future;
      schedule(continuation);
      release(branch);
      logging::push_promote_spawn_plus(sar->get_name());
      break;
    }
    case tag_join_plus: {
      interpreter* continuation = interp;
      auto stacks = fork_mark(interp->stack, [&] (char* _ar) {
        return is_splittable(_ar);
      });
      continuation->stack = stacks.first;
      interpreter* branch = new interpreter(stacks.second);
      assert(*block.variant_join_plus.getter(*sar, *par) == nullptr);
      *block.variant_join_plus.getter(*sar, *par) = continuation->get_incounter();
      sched::new_edge(branch, continuation);
      release(branch);
      logging::push_promote_join_plus(sar->get_name());
      break;
    }
    default: {
      assert(false);
      return;
    }
  }
  stats::on_promotion();
}
  
template <class Shared_activation_record>
sched::future get_dependency_of_join_minus(cfg_type<Shared_activation_record>& cfg, stack_type stack) {
  using private_activation_record = private_activation_record_of<Shared_activation_record>;
  auto& sar = peek_newest_shared_frame<Shared_activation_record>(stack);
  auto& par = peek_newest_private_frame<private_activation_record>(stack);
  auto& block = cfg.basic_blocks.at(par.trampoline.pred);
  assert(block.tag == tag_join_minus);
  auto r = *block.variant_join_minus.getter(sar, par);
  assert(r);
  return r;
}
  
static constexpr parallel_loop_id_type not_a_parallel_loop_id = -1;
 
template <class Shared_activation_record, class Private_activation_record>
class parallel_loop_private_activation_record : public private_activation_record {
public:
  
  using sar_type = Shared_activation_record;
  using par_type = Private_activation_record;
  
  parallel_loop_activation_record* loop_activation_record_of(parallel_loop_id_type id) {
    auto p = (par_type*)this;
    return p->_heartbeat_loop_activation_record_of(id);
  }
  
  void initialize_descriptors() {
    for (parallel_loop_id_type id = 0; id < sar_type::cfg.nb_loops(); id++) {
      sar_type::cfg.loop_descriptors[id].initializer(*(par_type*)this, loop_activation_record_of(id));
    }
  }
  
  parallel_loop_id_type get_id_of_current_parallel_loop() {
    return sar_type::cfg.loop_of.at(trampoline.pred);
  }

#ifndef NDEBUG
  std::vector<std::pair<int,int>> loop_activation_records() {
    std::vector<std::pair<int,int>> result;
    parallel_loop_id_type current = get_id_of_current_parallel_loop();
    if (current == not_a_parallel_loop_id) {
      return result;
    }
    for (parallel_loop_id_type id : sar_type::cfg.loop_descriptors[current].parents) {
      assert(id != not_a_parallel_loop_id);
      auto loop_ar = loop_activation_record_of(id);
      result.push_back(loop_ar->loop_range());
    }
    return result;
  }
#endif
  
  parallel_loop_id_type get_id_of_oldest_nonempty() {
    parallel_loop_id_type current = get_id_of_current_parallel_loop();
    if (current == not_a_parallel_loop_id) {
      return not_a_parallel_loop_id;
    }
    for (parallel_loop_id_type id : sar_type::cfg.loop_descriptors[current].parents) {
      assert(id != not_a_parallel_loop_id);
      if (loop_activation_record_of(id)->nb_strands() >= 2) {
        return id;
      }
    }
    if (loop_activation_record_of(current)->nb_strands() >= 2) {
      return current;
    } else {
      return not_a_parallel_loop_id;
    }
  }
  
  parallel_loop_activation_record* get_oldest_nonempty() {
    auto id = get_id_of_oldest_nonempty();
    if (id == not_a_parallel_loop_id) {
      return nullptr;
    } else {
      return loop_activation_record_of(id);
    }
  }
  
  sched::vertex*& get_join(parallel_loop_id_type id) {
    return loop_activation_record_of(id)->get_join();
  }
  
  sched::vertex*& get_join() {
    return get_join(get_id_of_current_parallel_loop());
  }
  
  int nb_strands() {
    if (trampoline.pred == exit_block_label) {
      return 0;
    }
    auto d = get_oldest_nonempty();
    if (d == nullptr) {
      return 1;
    }
    return std::max(1, d->nb_strands());
  }
  
  vertex_split_type split_join_trivial(interpreter* interp0, par_type* par0, sar_type* sar0,
                                       parallel_loop_id_type id, int nb) {
    assert(nb < interp0->nb_strands());
    assert(nb > 0);
    parallel_loop_descriptor_type<par_type>& lpdescr = sar_type::cfg.loop_descriptors[id];
    interpreter* interp1 = nullptr;
    interpreter* interp2 = nullptr;
    auto lpar0 = par0->loop_activation_record_of(id);
    par_type* par1 = par0;
    sched::vertex* join = lpar0->get_join();
    interpreter* interp00 = nullptr;
    if (join == nullptr) {
      join = interp0;
      auto stacks = split_stack(interp0->stack);
      interp0->stack = stacks.first;
      interpreter* interp01 = new interpreter(create_stack(sar0, par0));
      interp01->get_outset()->make_unary();
      par_type* par01 = &peek_newest_private_frame<par_type>(interp01->stack);
      par01->initialize_descriptors();
      auto lpar01 = par01->loop_activation_record_of(id);
      lpar0->split(lpar01, lpar0->nb_strands());
      interp1 = new interpreter(create_stack(sar0, par01));
      interp1->get_outset()->make_unary();
      par1 = &peek_newest_private_frame<par_type>(interp1->stack);
      par1->initialize_descriptors();
      auto lpar1 = par1->loop_activation_record_of(id);
      lpar01->split(lpar1, lpar01->nb_strands() - 1);
      interp0->stack = cactus::update_mark_stack(interp0->stack, [&] (char* _ar) {
        return pcfg::is_splittable(_ar);
      });
      interp01->stack = cactus::update_mark_stack(interp01->stack, [&] (char* _ar) {
        return pcfg::is_splittable(_ar);
      });
      nb--;
      std::swap(interp0->is_suspended, interp01->is_suspended);
      par0->trampoline = lpdescr.exit;
      par1->trampoline = lpdescr.entry;
      lpar1->get_join() = join;
      lpar01->get_join() = join;
      lpar0->get_join() = nullptr;
      interp1->stack = cactus::update_mark_stack(interp1->stack, [&] (char* _ar) {
        return pcfg::is_splittable(_ar);
      });
      sched::new_edge(interp01, join);
      sched::new_edge(interp1, join);
      if (empty_stack(stacks.second)) {
        interp00 = interp01;
      } else {
        interp00 = new interpreter(stacks.second);
        interp00->get_outset()->make_unary();
        sched::new_edge(interp00, interp01);
        release(interp01);
      }
      interp00->make_ready();
      interp1->make_ready();
    } else {
      interp1 = interp0;
    }
    interp2 = new interpreter(create_stack(sar0, par1));
    interp2->get_outset()->make_unary();
    par_type* par2 = &peek_newest_private_frame<par_type>(interp2->stack);
    par2->initialize_descriptors();
    auto lpar2 = par2->loop_activation_record_of(id);
    par1->loop_activation_record_of(id)->split(lpar2, nb);
    lpar2->get_join() = join;
    par2->trampoline = lpdescr.entry;
    interp1->stack = cactus::update_mark_stack(interp1->stack, [&] (char* _ar) {
      return pcfg::is_splittable(_ar);
    });
    interp2->stack = cactus::update_mark_stack(interp2->stack, [&] (char* _ar) {
      return pcfg::is_splittable(_ar);
    });
    sched::new_edge(interp2, join);
    interp2->make_ready();
    logging::push_promote_loop_split_join_trivial(sar0->get_name());
#ifdef DEBUG_HEARTBEAT_STACK
    check_stack(interp1->stack);
    check_stack(interp2->stack);
#endif
    return sched::make_vertex_split(interp00, interp1, interp2);
  }
  
  vertex_split_type split_join_associative_combine(interpreter* interp0, par_type* par0, sar_type* sar0,
                                                   parallel_loop_id_type id, int nb) {
    parallel_loop_descriptor_type<par_type>& lpdescr = sar_type::cfg.loop_descriptors[id];
    interpreter* interp1 = interp0;
    interpreter* interp2 = nullptr;
    par_type* par1 = par0;
    sar_type* sar1 = sar0;
    auto lpar1 = par1->loop_activation_record_of(id);
    interp2 = new interpreter(create_stack(sar1, par1));
    auto interp2_out = interp2->get_outset();
    auto interp2_future = interp2_out->make_chain_future();
    par_type* par2 = &peek_newest_private_frame<par_type>(interp2->stack);
    par2->initialize_descriptors();
    auto lpar2 = par2->loop_activation_record_of(id);
    lpar1->split(lpar2, nb);
    par2->trampoline = lpdescr.entry;
    auto& children = lpar1->get_children();
    if (! children) {
      children.reset(new children_record);
    }
    private_activation_record* destination = new par_type;
    children->futures.push_back(std::make_pair(interp2_future, destination));
    lpar2->get_destination() = destination;
    interp1->stack = cactus::update_mark_stack(interp1->stack, [&] (char* _ar) {
      return pcfg::is_splittable(_ar);
    });
    interp2->stack = cactus::update_mark_stack(interp2->stack, [&] (char* _ar) {
      return pcfg::is_splittable(_ar);
    });
    interp2->make_ready();
#ifdef DEBUG_HEARTBEAT_STACK
    check_stack(interp1->stack);
    check_stack(interp2->stack);
#endif
    logging::push_promote_loop_split_join_associative_combine(sar0->get_name());
    return sched::make_vertex_split(interp1, interp2);
  }
  
  vertex_split_type split(interpreter* interp0, int nb) {
    vertex_split_type result;
    par_type* par0 = &peek_marked_private_frame<par_type>(interp0->stack);
    sar_type* sar0 = &peek_marked_shared_frame<sar_type>(interp0->stack);
    parallel_loop_id_type id = par0->get_id_of_oldest_nonempty();
    switch (sar_type::cfg.loop_descriptors[id].join) {
      case parallel_loop_descriptor_type<par_type>::join_trivial: {
        result = split_join_trivial(interp0, par0, sar0, id, nb);
        break;
      }
      case parallel_loop_descriptor_type<par_type>::join_binary_associative_combine: {
        result = split_join_associative_combine(interp0, par0, sar0, id, nb);
        break;
      }
      default: {
        assert(false);
      }
    }
    return result;
  }
  
};

class parallel_for_activation_record : public parallel_loop_activation_record {
public:
  
  parallel_for_activation_record()
  : lo(nullptr), hi(nullptr), join(nullptr) { }
  
  parallel_for_activation_record(int& lo, int& hi)
  : lo(&lo), hi(&hi), join(nullptr) { }
  
  int* lo;
  
  int* hi;
  
  sched::vertex* join = nullptr;
  
  int nb_strands() {
    return *hi - *lo;
  }
  
  void split(parallel_loop_activation_record* _dest, int nb) {
    parallel_for_activation_record* dest = (parallel_for_activation_record*)_dest;
    int orig = nb_strands();
    assert(nb >= 0 && nb <= orig);
    int mid = (orig - nb) + *lo;
    *(dest->hi) = *hi;
    *hi = mid;
    *(dest->lo) = mid;
    assert((dest->nb_strands() == nb) && (nb_strands() + nb == orig));
  }
  
  sched::vertex*& get_join() {
    return join;
  }
  
  std::unique_ptr<children_record>& get_children() {
    assert(false); // impossible
    return dummy_children;
  }

  private_activation_record*& get_destination() {
    assert(false); // impossible
    return dummy_destination;
  }

#ifndef NDEBUG
  std::pair<int, int> loop_range() {
    return std::make_pair(*lo, *hi);
  }
#endif

};

class parallel_combine_activation_record : public parallel_loop_activation_record {
public:
  
  parallel_combine_activation_record()
  : lo(nullptr), hi(nullptr), destination(nullptr) { }
  
  parallel_combine_activation_record(int& lo, int& hi)
  : lo(&lo), hi(&hi), destination(nullptr) { }
  
  parallel_combine_activation_record(const parallel_combine_activation_record& other)
  : lo(other.lo), hi(other.hi), destination(nullptr) { }
  
  int* lo;
  
  int* hi;
  
  std::unique_ptr<children_record> children;
  
  private_activation_record* destination;
 
  int nb_strands() {
    return *hi - *lo;
  }
  
  void split(parallel_loop_activation_record* _dest, int nb) {
    parallel_combine_activation_record* dest = (parallel_combine_activation_record*)_dest;
    int orig = nb_strands();
    assert(nb >= 0 && nb <= orig);
    int mid = (orig - nb) + *lo;
    *(dest->hi) = *hi;
    *hi = mid;
    *(dest->lo) = mid;
    assert((dest->nb_strands() == nb) && (nb_strands() + nb == orig));
  }
  
  sched::vertex*& get_join() {
    assert(false); // impossible
    return dummy_join;
  }
  
  std::unique_ptr<children_record>& get_children() {
    return children;
  }
  
  private_activation_record*& get_destination() {
    return destination;
  }

#ifndef NDEBUG
  std::pair<int, int> loop_range() {
    return std::make_pair(*lo, *hi);
  }
#endif
  
};
  
} // end namespace
} // end namespace
} // end namespace

#endif /*! _HEARTBEAT_INTERPRETER_H_ */
