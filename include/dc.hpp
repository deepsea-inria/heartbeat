
#include <functional>
#include <vector>
#include <map>
#include <algorithm>
#include <memory>
#include <array>
#include <atomic>
#include <string>

#include "vertex.hpp"
#include "cactus-plus.hpp"
#include "stats.hpp"
#include "scheduler.hpp"
#include "interpreter.hpp"
#include "pcfg.hpp"
#include "cycles.hpp"
#include "grain.hpp"

#ifndef _HEARTBEAT_DC_H_
#define _HEARTBEAT_DC_H_

namespace heartbeat {
namespace edsl {
  
/*---------------------------------------------------------------------*/
/* Dag Calculus */
  
namespace dc {

static constexpr
logging::program_point_type dflt_ppt = logging::dflt_ppt;
    
using stmt_tag_type = enum {
  tag_stmt, tag_stmts, tag_cond, tag_exit_function, tag_exit_loop,
  tag_sequential_loop, tag_parallel_for_loop, tag_parallel_combine_loop,
  tag_spawn_join, tag_spawn2_join,
  tag_join_plus, tag_spawn_minus,
  tag_spawn_plus, tag_join_minus,
  tag_profile_statement,
  tag_none
};
  
template <class Shared_activation_record, class Private_activation_record>
class stmt_type {
public:

  using sar_type = Shared_activation_record;
  using par_type = Private_activation_record;
  
  using unconditional_jump_code_type = std::function<void(sar_type&, par_type&)>;
  using predicate_code_type = std::function<bool(sar_type&, par_type&)>;
  using procedure_call_code_type = std::function<pcfg::stack_type(sar_type&, par_type&, pcfg::cactus::parent_link_type, pcfg::stack_type)>;
  using incounter_getter_code_type = std::function<sched::incounter**(sar_type&, par_type&)>;
  using future_getter_code_type = std::function<sched::future*(sar_type&, par_type&)>;
  using conds_type = std::vector<std::pair<predicate_code_type, stmt_type>>;
  using parallel_loop_range_getter_type = std::function<std::pair<int*, int*>(par_type&)>;
  using parallel_loop_combine_initializer_type = std::function<void(sar_type&, par_type&)>;
  using parallel_combining_operator_type = std::function<void(sar_type&, par_type&, par_type&)>;
  using loop_range_getter_type = std::function<std::pair<int*, int*>(sar_type&, par_type&)>;
  using leaf_loop_body_type = std::function<void(sar_type&, par_type&, int, int)>;
  using profile_prefix_getter_type = std::function<std::pair<uint64_t*, uint64_t*>(sar_type&, par_type&)>;
  using profile_report_code_type = std::function<void(sar_type&, par_type&, uint64_t, uint64_t)>;
  
  stmt_tag_type tag;
  
  union {
    struct {
      unconditional_jump_code_type code;
    } variant_stmt;
    struct {
      std::vector<stmt_type> stmts;
    } variant_stmts;
    struct {
      conds_type conds;
      std::unique_ptr<stmt_type> otherwise;
    } variant_cond;
    struct {
      // nothing to put here
    } variant_exit_function;
    struct {
      // nothing to put here
    } variant_exit_loop;
    struct {
      predicate_code_type predicate;
      std::unique_ptr<stmt_type> body;
    } variant_sequential_loop;
    struct {
      predicate_code_type predicate;
      parallel_loop_range_getter_type getter;
      std::unique_ptr<stmt_type> body;
    } variant_parallel_for_loop;
    struct {
      predicate_code_type predicate;
      parallel_loop_range_getter_type getter;
      parallel_loop_combine_initializer_type initialize;
      parallel_combining_operator_type combine;
      std::unique_ptr<stmt_type> body;
    } variant_parallel_combine_loop;
    struct {
      procedure_call_code_type code;
    } variant_spawn_join;
    struct {
      procedure_call_code_type code1;
      procedure_call_code_type code2;
    } variant_spawn2_join;
    struct {
      procedure_call_code_type code;
      incounter_getter_code_type getter;
    } variant_join_plus;
    struct {
      procedure_call_code_type code;
      incounter_getter_code_type getter;
    } variant_spawn_minus;
    struct {
      procedure_call_code_type code;
      future_getter_code_type getter;
    } variant_spawn_plus;
    struct {
      future_getter_code_type getter;
    } variant_join_minus;
    struct {
      profile_prefix_getter_type getter;
      std::unique_ptr<stmt_type> body;
      profile_report_code_type reporter;
    } variant_profile_statement;
  };
  
  stmt_type() : tag(tag_none) { }
  
private:
  
  void copy_constructor(stmt_type const& other) {
    tag = other.tag;
    switch (tag) {
      case tag_stmt: {
        new (&variant_stmt.code) unconditional_jump_code_type(other.variant_stmt.code);
        break;
      }
      case tag_stmts: {
        new (&variant_stmts.stmts) std::vector<stmt_type>(other.variant_stmts.stmts);
        break;
      }
      case tag_cond: {
        new (&variant_cond.conds) conds_type(other.variant_cond.conds);
        auto p = new stmt_type;
        p->copy_constructor(*other.variant_cond.otherwise);
        new (&variant_cond.otherwise) std::unique_ptr<stmt_type>(p);
        break;
      }
      case tag_exit_function: {
        // nothing to put here
        break;
      }
      case tag_exit_loop: {
        // nothing to put here
        break;
      }
      case tag_sequential_loop: {
        new (&variant_sequential_loop.predicate) predicate_code_type(other.variant_sequential_loop.predicate);
        auto p = new stmt_type;
        p->copy_constructor(*other.variant_sequential_loop.body);
        new (&variant_sequential_loop.body) std::unique_ptr<stmt_type>(p);
        break;
      }
      case tag_parallel_for_loop: {
        new (&variant_parallel_for_loop.predicate) predicate_code_type(other.variant_parallel_for_loop.predicate);
        new (&variant_parallel_for_loop.getter) parallel_loop_range_getter_type(other.variant_parallel_for_loop.getter);
        auto p = new stmt_type;
        p->copy_constructor(*other.variant_parallel_for_loop.body);
        new (&variant_parallel_for_loop.body) std::unique_ptr<stmt_type>(p);
        break;
      }
      case tag_parallel_combine_loop: {
        new (&variant_parallel_combine_loop.predicate) predicate_code_type(other.variant_parallel_combine_loop.predicate);
        new (&variant_parallel_combine_loop.getter) parallel_loop_range_getter_type(other.variant_parallel_combine_loop.getter);
        new (&variant_parallel_combine_loop.initialize) parallel_loop_combine_initializer_type(other.variant_parallel_combine_loop.initialize);
        new (&variant_parallel_combine_loop.combine) parallel_combining_operator_type(other.variant_parallel_combine_loop.combine);
        auto p = new stmt_type;
        p->copy_constructor(*other.variant_parallel_combine_loop.body);
        new (&variant_parallel_combine_loop.body) std::unique_ptr<stmt_type>(p);
        break;
      }
      case tag_spawn_join: {
        new (&variant_spawn_join.code) procedure_call_code_type(other.variant_spawn_join.code);
        break;
      }
      case tag_spawn2_join: {
        new (&variant_spawn2_join.code1) procedure_call_code_type(other.variant_spawn2_join.code1);
        new (&variant_spawn2_join.code2) procedure_call_code_type(other.variant_spawn2_join.code2);
        break;
      }
      case tag_join_plus: {
        new (&variant_join_plus.code) procedure_call_code_type(other.variant_join_plus.code);
        new (&variant_join_plus.getter) incounter_getter_code_type(other.variant_join_plus.getter);
        break;
      }
      case tag_spawn_minus: {
        new (&variant_spawn_minus.code) procedure_call_code_type(other.variant_spawn_minus.code);
        new (&variant_spawn_minus.getter) incounter_getter_code_type(other.variant_spawn_minus.getter);
        break;
      }
      case tag_spawn_plus: {
        new (&variant_spawn_plus.code) procedure_call_code_type(other.variant_spawn_plus.code);
        new (&variant_spawn_plus.getter) future_getter_code_type(other.variant_spawn_plus.getter);
        break;
      }
      case tag_join_minus: {
        new (&variant_join_minus.getter) future_getter_code_type(other.variant_join_minus.getter);
        break;
      }
      case tag_profile_statement: {
        new (&variant_profile_statement.getter) profile_prefix_getter_type(other.variant_profile_statement.getter);
        auto p = new stmt_type;
        p->copy_constructor(*other.variant_profile_statement.body);
        new (&variant_profile_statement.body) std::unique_ptr<stmt_type>(p);
        new (&variant_profile_statement.reporter) profile_report_code_type(other.variant_profile_statement.reporter);
        break;
      }
      default: {
        break;
      }
    }
  }
  
public:
  
  stmt_type(stmt_type const& other) {
    copy_constructor(other);
  }
  
private:
  
  void move_constructor(stmt_type&& other) {
    tag = other.tag;
    other.tag = tag_none;
    switch (tag) {
      case tag_stmt: {
        new (&variant_stmt.code) unconditional_jump_code_type;
        variant_stmt.code = std::move(other.variant_stmt.code);
        break;
      }
      case tag_stmts: {
        new (&variant_stmts.stmts) std::vector<stmt_type>;
        variant_stmts.stmts = std::move(other.variant_stmts.stmts);
        break;
      }
      case tag_cond: {
        new (&variant_cond.conds) conds_type;
        variant_cond.conds = std::move(other.variant_cond.conds);
        new (&variant_cond.otherwise) std::unique_ptr<stmt_type>;
        variant_cond.otherwise = std::move(other.variant_cond.otherwise);
        break;
      }
      case tag_exit_function: {
        // nothing to put here
        break;
      }
      case tag_exit_loop: {
        // nothing to put here
        break;
      }
      case tag_sequential_loop: {
        new (&variant_sequential_loop.predicate) predicate_code_type;
        variant_sequential_loop.predicate = std::move(other.variant_sequential_loop.predicate);
        new (&variant_sequential_loop.body) std::unique_ptr<stmt_type>;
        variant_sequential_loop.body = std::move(other.variant_sequential_loop.body);
        break;
      }
      case tag_parallel_for_loop: {
        new (&variant_parallel_for_loop.predicate) predicate_code_type;
        new (&variant_parallel_for_loop.getter) parallel_loop_range_getter_type;
        variant_parallel_for_loop.getter = std::move(other.variant_parallel_for_loop.getter);
        new (&variant_parallel_for_loop.body) std::unique_ptr<stmt_type>;
        variant_parallel_for_loop.body = std::move(other.variant_parallel_for_loop.body);
        break;
      }
      case tag_parallel_combine_loop: {
        new (&variant_parallel_combine_loop.predicate) predicate_code_type;
        new (&variant_parallel_combine_loop.getter) parallel_loop_range_getter_type;
        variant_parallel_combine_loop.getter = std::move(other.variant_parallel_combine_loop.getter);
        new (&variant_parallel_combine_loop.initialize) parallel_loop_combine_initializer_type;
        variant_parallel_combine_loop.initialize = std::move(other.variant_parallel_combine_loop.initialize);
        new (&variant_parallel_combine_loop.combine) parallel_combining_operator_type;
        variant_parallel_combine_loop.combine = std::move(other.variant_parallel_combine_loop.combine);
        new (&variant_parallel_combine_loop.body) std::unique_ptr<stmt_type>;
        variant_parallel_combine_loop.body = std::move(other.variant_parallel_combine_loop.body);
        break;
      }
      case tag_spawn_join: {
        new (&variant_spawn_join.code) procedure_call_code_type;
        variant_spawn_join.code = std::move(other.variant_spawn_join.code);
        break;
      }
      case tag_spawn2_join: {
        new (&variant_spawn2_join.code1) procedure_call_code_type;
        variant_spawn2_join.code1 = std::move(other.variant_spawn2_join.code1);
        new (&variant_spawn2_join.code2) procedure_call_code_type;
        variant_spawn2_join.code2 = std::move(other.variant_spawn2_join.code2);
        break;
      }
      case tag_join_plus: {
        new (&variant_join_plus.code) procedure_call_code_type;
        variant_join_plus.code = std::move(other.variant_join_plus.code);
        new (&variant_join_plus.getter) incounter_getter_code_type;
        variant_join_plus.getter = std::move(other.variant_join_plus.getter);
        break;
      }
      case tag_spawn_minus: {
        new (&variant_spawn_minus.code) procedure_call_code_type;
        variant_spawn_minus.code = std::move(other.variant_spawn_minus.code);
        new (&variant_spawn_minus.getter) incounter_getter_code_type;
        variant_spawn_minus.getter = std::move(other.variant_spawn_minus.getter);
        break;
      }
      case tag_spawn_plus: {
        new (&variant_spawn_plus.code) procedure_call_code_type;
        variant_spawn_plus.code = std::move(other.variant_spawn_plus.code);
        new (&variant_spawn_plus.getter) future_getter_code_type;
        variant_spawn_plus.getter = std::move(other.variant_spawn_plus.getter);
        break;
      }
      case tag_join_minus: {
        new (&variant_join_minus.getter) future_getter_code_type;
        variant_join_minus.getter = std::move(other.variant_join_minus.getter);
        break;
      }
      case tag_profile_statement: { 
        new (&variant_profile_statement.getter) profile_prefix_getter_type;
        variant_profile_statement.getter = std::move(other.variant_profile_statement.getter);
        new (&variant_profile_statement.body) std::unique_ptr<stmt_type>;
        variant_profile_statement.body = std::move(other.variant_profile_statement.body);
        new (&variant_profile_statement.reporter) profile_report_code_type;
        variant_profile_statement.reporter = std::move(other.variant_profile_statement.reporter);
        break;
      }
      default: {
        break;
      }
    }
  }
  
public:
  
  stmt_type(stmt_type&& other) {
    move_constructor(std::move(other));
  }
  
  ~stmt_type() {
    switch (tag) {
      case tag_stmt: {
        variant_stmt.code.~unconditional_jump_code_type();
        break;
      }
      case tag_stmts: {
        using st = std::vector<stmt_type>;
        variant_stmts.stmts.~st();
        break;
      }
      case tag_cond: {
        variant_cond.conds.~conds_type();
        using st = std::unique_ptr<stmt_type>;
        variant_cond.otherwise.~st();
        break;
      }
      case tag_exit_function: {
        // nothing to put here
        break;
      }
      case tag_exit_loop: {
        // nothing to put here
        break;
      }
      case tag_sequential_loop: {
        variant_sequential_loop.predicate.~predicate_code_type();
        using st = std::unique_ptr<stmt_type>;
        variant_sequential_loop.body.~st();
        break;
      }
      case tag_parallel_for_loop: {
        variant_parallel_for_loop.predicate.~predicate_code_type();
        variant_parallel_for_loop.getter.~parallel_loop_range_getter_type();
        using st = std::unique_ptr<stmt_type>;
        variant_parallel_for_loop.body.~st();
        break;
      }
      case tag_parallel_combine_loop: {
        variant_parallel_combine_loop.predicate.~predicate_code_type();
        variant_parallel_combine_loop.getter.~parallel_loop_range_getter_type();
        variant_parallel_combine_loop.initialize.~parallel_loop_combine_initializer_type();
        variant_parallel_combine_loop.combine.~parallel_combining_operator_type();
        using st = std::unique_ptr<stmt_type>;
        variant_parallel_combine_loop.body.~st();
        break;
      }
      case tag_spawn_join: {
        variant_spawn_join.code.~procedure_call_code_type();
        break;
      }
      case tag_spawn2_join: {
        variant_spawn2_join.code1.~procedure_call_code_type();
        variant_spawn2_join.code2.~procedure_call_code_type();
        break;
      }
      case tag_join_plus: {
        variant_join_plus.code.~procedure_call_code_type();
        variant_join_plus.getter.~incounter_getter_code_type();
        break;
      }
      case tag_spawn_minus: {
        variant_spawn_minus.code.~procedure_call_code_type();
        variant_spawn_minus.getter.~incounter_getter_code_type();
        break;
      }
      case tag_spawn_plus: {
        variant_spawn_plus.code.~procedure_call_code_type();
        variant_spawn_plus.getter.~future_getter_code_type();
        break;
      }
      case tag_join_minus: {
        variant_join_minus.getter.~future_getter_code_type();
        break;
      }
      case tag_profile_statement: { 
        variant_profile_statement.getter.~profile_prefix_getter_type();
        using st = std::unique_ptr<stmt_type>;
        variant_profile_statement.body.~st();
        variant_profile_statement.reporter.~profile_report_code_type();
        break;
      }
      default: {
        break;
      }
    }
  }
  
  stmt_type& operator=(stmt_type const& other) {
    copy_constructor(other);
    return *this;
  }
  
  stmt_type& operator=(stmt_type&& other) {
    move_constructor(std::move(other));
    return *this;
  }
  
  static
  stmt_type stmt(unconditional_jump_code_type code) {
    stmt_type s;
    s.tag = tag_stmt;
    new (&s.variant_stmt.code) unconditional_jump_code_type(code);
    return s;
  }
  
  static
  stmt_type stmts(std::vector<stmt_type> stmts) {
    stmt_type s;
    s.tag = tag_stmts;
    new (&s.variant_stmts.stmts) std::vector<stmt_type>(stmts);
    return s;
  }
  
  static
  stmt_type cond(conds_type conds, stmt_type otherwise) {
    stmt_type s;
    s.tag = tag_cond;
    new (&s.variant_cond.conds) conds_type(conds);
    new (&s.variant_cond.otherwise) std::unique_ptr<stmt_type>(new stmt_type(otherwise));
    return s;
  }
  
  static
  stmt_type exit_function() {
    stmt_type s;
    s.tag = tag_exit_function;
    return s;
  }

  static
  stmt_type exit_loop() {
    stmt_type s;
    s.tag = tag_exit_loop;
    return s;
  }

  static
  stmt_type sequential_loop(predicate_code_type predicate, stmt_type body) {
    stmt_type s;
    s.tag = tag_sequential_loop;
    new (&s.variant_sequential_loop.predicate) predicate_code_type(predicate);
    new (&s.variant_sequential_loop.body) std::unique_ptr<stmt_type>(new stmt_type(body));
    return s;
  }
  
  static
  stmt_type parallel_for_loop(predicate_code_type predicate, parallel_loop_range_getter_type getter, stmt_type body) {
    stmt_type s;
    s.tag = tag_parallel_for_loop;
    new (&s.variant_parallel_for_loop.predicate) predicate_code_type(predicate);
    new (&s.variant_parallel_for_loop.getter) parallel_loop_range_getter_type(getter);
    new (&s.variant_parallel_for_loop.body) std::unique_ptr<stmt_type>(new stmt_type(body));
    return s;
  }
  
  static
  stmt_type parallel_combine_loop(predicate_code_type predicate,
                                  parallel_loop_range_getter_type getter,
                                  parallel_loop_combine_initializer_type initialize,
                                  parallel_combining_operator_type combine,
                                  stmt_type body) {
    stmt_type s;
    s.tag = tag_parallel_combine_loop;
    new (&s.variant_parallel_combine_loop.predicate) predicate_code_type(predicate);
    new (&s.variant_parallel_combine_loop.getter) parallel_loop_range_getter_type(getter);
    new (&s.variant_parallel_combine_loop.initialize) parallel_loop_combine_initializer_type(initialize);
    new (&s.variant_parallel_combine_loop.combine) parallel_combining_operator_type(combine);
    new (&s.variant_parallel_combine_loop.body) std::unique_ptr<stmt_type>(new stmt_type(body));
    return s;
  }
  
  static
  stmt_type spawn_join(procedure_call_code_type code) {
    stmt_type s;
    s.tag = tag_spawn_join;
    new (&s.variant_spawn_join) procedure_call_code_type(code);
    return s;
  }
  
  static
  stmt_type spawn2_join(procedure_call_code_type code1, procedure_call_code_type code2) {
    stmt_type s;
    s.tag = tag_spawn2_join;
    new (&s.variant_spawn2_join.code1) procedure_call_code_type(code1);
    new (&s.variant_spawn2_join.code2) procedure_call_code_type(code2);
    return s;
  }
  
  static
  stmt_type join_plus(procedure_call_code_type code, incounter_getter_code_type getter) {
    stmt_type s;
    s.tag = tag_join_plus;
    new (&s.variant_join_plus.code) procedure_call_code_type(code);
    new (&s.variant_join_plus.getter) incounter_getter_code_type(getter);
    return s;
  }
  
  static
  stmt_type spawn_minus(procedure_call_code_type code, incounter_getter_code_type getter) {
    stmt_type s;
    s.tag = tag_spawn_minus;
    new (&s.variant_spawn_minus.code) procedure_call_code_type(code);
    new (&s.variant_spawn_minus.getter) incounter_getter_code_type(getter);
    return s;
  }
  
  static
  stmt_type spawn_plus(procedure_call_code_type code, future_getter_code_type getter) {
    stmt_type s;
    s.tag = tag_spawn_plus;
    new (&s.variant_spawn_plus.code) procedure_call_code_type(code);
    new (&s.variant_spawn_plus.getter) future_getter_code_type(getter);
    return s;
  }
  
  static
  stmt_type join_minus(future_getter_code_type getter) {
    stmt_type s;
    s.tag = tag_join_minus;
    new (&s.variant_join_minus.getter) future_getter_code_type(getter);
    return s;
  }

  static
  stmt_type profile_statement(profile_prefix_getter_type getter,
                              stmt_type body,
                              profile_report_code_type reporter) {
    stmt_type s;
    s.tag = tag_profile_statement;
    new (&s.variant_profile_statement.getter) profile_prefix_getter_type(getter);
    new (&s.variant_profile_statement.body) std::unique_ptr<stmt_type>(new stmt_type(body));
    new (&s.variant_profile_statement.reporter) profile_report_code_type(reporter);
    return s;
  }
  
  static
  stmt_type mk_if(predicate_code_type pred, stmt_type branch1, stmt_type branch2) {
    return cond({ std::make_pair(pred, branch1) }, branch2);
  }
  
  static
  stmt_type mk_if(predicate_code_type pred, stmt_type branch1) {
    return mk_if(pred, branch1, stmt([] (sar_type&, par_type&) { }));
  }

  template <int threshold=grain::automatic, class Unconditional_jump_code_type>
  static
  stmt_type sequential_loop(predicate_code_type predicate, Unconditional_jump_code_type body,
                            int line_nb=dflt_ppt.line_nb, const char* source_fname=dflt_ppt.source_fname) {
    using controller_type = grain::controller<threshold, Unconditional_jump_code_type>;
    controller_type::set_ppt(line_nb, source_fname);
    return sequential_loop(predicate, stmt([=] (sar_type& s, par_type& p) {
      auto lg_lt = controller_type::predict_lg_nb_iterations();
      auto lt = controller_type::predict_nb_iterations(lg_lt);
      int i = 0;
      for (; i < lt; i++) {
        if (predicate(s, p)) {
          body(s, p);
        } else {
          break;
        }
      }
      controller_type::register_callback(lg_lt, i);
    }));
  }

  using loop_direction_type = enum { forward_loop, backward_loop };

  template <int threshold=grain::automatic, class Leaf_loop_body_type>
  static
  stmt_type sequential_loop(unconditional_jump_code_type initializer,
                            loop_range_getter_type getter,
                            Leaf_loop_body_type body,
                            loop_direction_type direction = forward_loop,
                            int line_nb=dflt_ppt.line_nb, const char* source_fname=dflt_ppt.source_fname) {
    using controller_type = grain::controller<threshold, Leaf_loop_body_type>;
    controller_type::set_ppt(line_nb, source_fname);
    return stmts({
      stmt(initializer),
      sequential_loop([=] (sar_type& s, par_type& p) {
        auto rng = getter(s, p);
        return *rng.first != *rng.second;
      }, stmt([=] (sar_type& s, par_type& p) {
        auto rng = getter(s, p);
        auto lo = *rng.first;
        auto hi = *rng.second;
        auto lo2 = 0;
        auto hi2 = 0;
        auto lg_lt = controller_type::predict_lg_nb_iterations();
        auto lt = controller_type::predict_nb_iterations(lg_lt);
        if (direction == forward_loop) {
          auto mid = std::min(lo + lt, hi);
          *rng.first = mid;
          lo2 = lo;
          hi2 = mid;
        } else { // backward_loop
          auto mid = std::max(lo, hi - lt);
          *rng.second = mid;
          lo2 = mid;
          hi2 = hi;
        }
        body(s, p, lo2, hi2);
        controller_type::register_callback(lg_lt, hi2 - lo2);
      }))
    });
  }

  static
  stmt_type parallel_for_loop(parallel_loop_range_getter_type getter, stmt_type body) {
    return parallel_for_loop([=] (sar_type& s, par_type& p) {
      auto rng = getter(p);
      return *rng.first != *rng.second;
    }, getter, body);
  }
  
  template <int threshold=grain::automatic, class Leaf_loop_body_type>
  static
  stmt_type parallel_for_loop(unconditional_jump_code_type initializer,
                              parallel_loop_range_getter_type getter,
                              Leaf_loop_body_type body,
                              int line_nb=dflt_ppt.line_nb, const char* source_fname=dflt_ppt.source_fname) {
    using controller_type = grain::controller<threshold, Leaf_loop_body_type>;
    controller_type::set_ppt(line_nb, source_fname);
    return stmts({
      stmt(initializer),
      parallel_for_loop(getter, stmt([=] (sar_type& s, par_type& p) {
        auto lg_lt = controller_type::predict_lg_nb_iterations();
        auto lt = controller_type::predict_nb_iterations(lg_lt);
        auto rng = getter(p);
        auto lo = *rng.first;
        auto mid = std::min(lo + lt, *rng.second);
        *rng.first = mid;
        body(s, p, lo, mid);
        controller_type::register_callback(lg_lt, mid - lo);
      }))
    });
  }

  static
  stmt_type parallel_combine_loop(parallel_loop_range_getter_type getter,
                                  parallel_loop_combine_initializer_type initialize,
                                  parallel_combining_operator_type combine,
                                  stmt_type body) {
    return parallel_combine_loop([=] (sar_type& s, par_type& p) {
      auto rng = getter(p);
      return *rng.first != *rng.second;
    }, getter, initialize, combine, body);
  }

  template <int threshold=grain::automatic, class Leaf_loop_body_type>
  static
  stmt_type parallel_combine_loop(unconditional_jump_code_type initializer,
                                  parallel_loop_range_getter_type getter,
                                  parallel_loop_combine_initializer_type combine_initializer,
                                  parallel_combining_operator_type combine,
                                  Leaf_loop_body_type body,
                                  int line_nb=dflt_ppt.line_nb, const char* source_fname=dflt_ppt.source_fname) {
    using controller_type = grain::controller<threshold, Leaf_loop_body_type>;
    controller_type::set_ppt(line_nb, source_fname);
    return stmts({
      stmt(initializer),
      parallel_combine_loop(getter, combine_initializer, combine,
                            stmt([=] (sar_type& s, par_type& p) {
        auto rng = getter(p);
        auto lo = *rng.first;
        auto lg_lt = controller_type::predict_lg_nb_iterations();
        auto lt = controller_type::predict_nb_iterations(lg_lt);
        auto mid = std::min(lo + lt, *rng.second);
        *rng.first = mid;
        body(s, p, lo, mid);
        controller_type::register_callback(lg_lt, mid - lo);
      }))
    });    
  }
  
};
  
template <class Shared_activation_record, class Private_activation_record>
class linearize {
private:

  using sar = Shared_activation_record;
  using par = Private_activation_record;
  using stmt_type = stmt_type<sar, par>;
  
  using basic_block_label_type = pcfg::basic_block_label_type;
  using basic_block_type = pcfg::basic_block_type<sar, par>;
  using bbt = basic_block_type;
  using lt = basic_block_label_type;
  using block_map_type = std::map<basic_block_label_type, basic_block_type>;
  
  using parallel_loop_id_type = pcfg::parallel_loop_id_type;
  using parallel_loop_of_basic_block_type = std::map<basic_block_label_type, parallel_loop_id_type>;
  using parallel_loop_descriptor_type = pcfg::parallel_loop_descriptor_type<Private_activation_record>;
  using parallel_loop_descriptor_table_type = std::map<parallel_loop_id_type, parallel_loop_descriptor_type>;
  
  static
  lt next_block_label;
  
  static
  parallel_loop_id_type next_parallel_loop_id;
  
  using configuration_type = struct {
    block_map_type blocks;
    parallel_loop_of_basic_block_type loop_of_basic_block;
    parallel_loop_descriptor_table_type descriptor_of_loop;
  };
  
  static
  configuration_type transform(stmt_type stmt,
                               lt entry, lt exit,
                               lt loop_exit_block,
                               std::vector<parallel_loop_id_type> loop_scope,
                               configuration_type config) {
    configuration_type result = config;
    auto new_label = [&] {
      return next_block_label++;
    };
    auto add_block = [&] (lt label, bbt block) {
      assert(result.blocks.find(label) == result.blocks.end());
      result.blocks.insert(std::make_pair(label, block));
      if (! loop_scope.empty()) {
        result.loop_of_basic_block.insert(std::make_pair(label, loop_scope.back()));
      }
    };
    auto new_parallel_loop_label = [&] {
      return next_parallel_loop_id++;
    };
    auto add_parallel_loop = [&] (parallel_loop_id_type label, parallel_loop_descriptor_type d) {
      assert(result.descriptor_of_loop.find(label) == result.descriptor_of_loop.end());
      result.descriptor_of_loop.insert(std::make_pair(label, d));
    };
    switch (stmt.tag) {
      case tag_stmt: {
        add_block(entry, bbt::unconditional_jump(stmt.variant_stmt.code, exit));
        break;
      }
      case tag_stmts: {
        std::vector<lt> entries;
        auto stmts = stmt.variant_stmts.stmts;
        for (int i = 0; i < stmts.size(); i++) {
          entries.push_back((i == 0) ? entry : new_label());
        }
        int i = 0;
        for (auto& s : stmts) {
          lt entry = entries.at(i);
          lt exit2 = (i + 1 == entries.size()) ? exit : entries.at(i + 1);
          result = transform(s, entry, exit2, loop_exit_block, loop_scope, result);
          i++;
        }
        break;
      }
      case tag_cond: {
        auto nb_branches = stmt.variant_cond.conds.size();
        std::vector<typename stmt_type::predicate_code_type> predicates(nb_branches);
        std::vector<lt> entries(nb_branches);
        std::vector<stmt_type> stmts(nb_branches);
        int i = 0;
        for (auto& c : stmt.variant_cond.conds) {
          predicates[i] = c.first;
          stmts[i] = c.second;
          entries[i] = new_label();
          i++;
        }
        auto otherwise_label = new_label();
        i = 0;
        for (auto& s : stmts) {
          result = transform(s, entries.at(i), exit, loop_exit_block, loop_scope, result);
          i++;
        }
        entries.push_back(otherwise_label);
        result = transform(*stmt.variant_cond.otherwise, otherwise_label, exit, loop_exit_block, loop_scope, result);
        auto selector = [predicates, entries] (sar& s, par& p) {
          int i = 0;
          for (auto& pred : predicates) {
            if (pred(s, p)) {
              break;
            }
            i++;
          }
          return entries[i];
        };
        add_block(entry, bbt::conditional_jump(selector));
        break;
      }
      case tag_exit_function: {
        add_block(entry, bbt::unconditional_jump([] (sar&, par&) { }, pcfg::exit_block_label));
        break;
      }
      case tag_exit_loop: {
        if (loop_exit_block == -1) {
          atomic::die("bogus loop exit\n");
        }
        add_block(entry, bbt::unconditional_jump([] (sar&, par&) { }, loop_exit_block));
        break;
      }
      case tag_sequential_loop: {
        auto header_label = entry;
        auto body_label = new_label();
        auto predicate = stmt.variant_sequential_loop.predicate;
        auto selector = [predicate, body_label, exit] (sar& s, par& p) {
          return predicate(s, p) ? body_label : exit;
        };
        add_block(header_label, bbt::conditional_jump(selector));
        result = transform(*stmt.variant_sequential_loop.body, body_label, header_label, exit, loop_scope, result);
        break;
      }
      case tag_parallel_for_loop: {
        parallel_loop_id_type loop_label = new_parallel_loop_label();
        loop_scope.push_back(loop_label);
        parallel_loop_descriptor_type descriptor;
        descriptor.join = parallel_loop_descriptor_type::join_trivial;
        descriptor.entry = { .pred=entry, .succ=entry };
        descriptor.exit = { .pred=exit, .succ=exit };
        descriptor.parents = loop_scope;
        auto getter = stmt.variant_parallel_for_loop.getter;
        descriptor.initializer = [getter] (Private_activation_record& p, pcfg::parallel_loop_activation_record* _ar) {
          pcfg::parallel_for_activation_record& ar = *((pcfg::parallel_for_activation_record*)_ar);
          new (&ar) pcfg::parallel_for_activation_record;
          std::pair<int*, int*> range = getter(p);
          ar.lo = range.first;
          ar.hi = range.second;
        };
        add_parallel_loop(loop_label, descriptor);
        auto header_label = entry;
        auto body_label = new_label();
        auto footer_label = new_label();
        auto predicate = stmt.variant_parallel_for_loop.predicate;
        auto selector = [predicate, body_label, footer_label] (sar& s, par& p) {
          return predicate(s, p) ? body_label : footer_label;
        };
        add_block(header_label, bbt::conditional_jump(selector));
        auto footer_selector = [loop_label, exit] (sar&, par& p) {
          return p.get_join(loop_label) == nullptr ? exit : pcfg::exit_block_label;
        };
        add_block(footer_label, bbt::conditional_jump(footer_selector));
        result = transform(*stmt.variant_parallel_for_loop.body, body_label, header_label, loop_exit_block, loop_scope, result);
        break;
      }
      case tag_parallel_combine_loop: {
        parallel_loop_id_type loop_label = new_parallel_loop_label();
        loop_scope.push_back(loop_label);
        parallel_loop_descriptor_type descriptor;
        descriptor.join = parallel_loop_descriptor_type::join_binary_associative_combine;
        descriptor.entry = { .pred=entry, .succ=entry };
        descriptor.exit = { .pred=exit, .succ=exit };
        descriptor.parents = loop_scope;
        auto getter = stmt.variant_parallel_combine_loop.getter;
        descriptor.initializer = [getter] (Private_activation_record& p, pcfg::parallel_loop_activation_record* _ar) {
          pcfg::parallel_combine_activation_record& ar = *((pcfg::parallel_combine_activation_record*)_ar);
          new (&ar) pcfg::parallel_combine_activation_record;
          std::pair<int*, int*> range = getter(p);
          ar.lo = range.first;
          ar.hi = range.second;
        };
        add_parallel_loop(loop_label, descriptor);
        auto initialize_label = entry;
        auto header_label = new_label();
        auto body_label = new_label();
        auto children_loop_finalize_label = new_label();
        auto children_loop_header_label = new_label();
        auto children_loop_body0_label = new_label();
        auto children_loop_body1_label = new_label();
        auto parent_check_label = new_label();
        add_block(initialize_label, bbt::unconditional_jump(stmt.variant_parallel_combine_loop.initialize, header_label));
        auto predicate = stmt.variant_parallel_combine_loop.predicate;
        auto selector = [predicate, body_label, children_loop_finalize_label] (sar& s, par& p) {
          return predicate(s, p) ? body_label : children_loop_finalize_label;
        };
        add_block(header_label, bbt::conditional_jump(selector));
        add_block(children_loop_finalize_label, bbt::unconditional_jump([loop_label] (sar&, par& p) {
          auto& children = p.loop_activation_record_of(loop_label)->get_children();
          if (children) {
            children->nb = (int)children->futures.size();
          }
        }, children_loop_header_label));
        auto children_loop_header = [loop_label, parent_check_label, children_loop_body0_label] (sar&, par& p) {
          auto& children = p.loop_activation_record_of(loop_label)->get_children();
          if (children) {
            if (children->nb <= 0) {
              children.reset();
              return parent_check_label;
            }
            return children_loop_body0_label;
          } else {
            return parent_check_label;
          }
        };
        add_block(children_loop_header_label, bbt::conditional_jump(children_loop_header));
        auto children_loop_body0 = [loop_label] (sar&, par& p) {
          auto& children = p.loop_activation_record_of(loop_label)->get_children();
          assert(children);
          auto i = children->nb - 1;
          return &(children->futures.at(i).first);
        };
        add_block(children_loop_body0_label, bbt::join_minus(children_loop_body0, children_loop_body1_label));
        auto combine = stmt.variant_parallel_combine_loop.combine;
        auto children_loop_body1 = [loop_label, combine] (sar& s, par& p) {
          auto& children = p.loop_activation_record_of(loop_label)->get_children();
          assert(children);
          auto i = children->nb - 1;
          par* c = (par*)children->futures.at(i).second;
          combine(s, *c, p);
          delete c;
          children->futures.at(i).second = nullptr; // to avoid the dangling pointer
          children->nb--;
        };
        add_block(children_loop_body1_label, bbt::unconditional_jump(children_loop_body1, children_loop_header_label));
        auto parent_check = [loop_label, exit] (sar& s, par& p) {
          auto destination = (par*)p.loop_activation_record_of(loop_label)->get_destination();
          if (destination != nullptr) {
            new (destination) par(p);
            return pcfg::exit_block_label;
          }
          return exit;
        };
        add_block(parent_check_label, bbt::conditional_jump(parent_check));
        result = transform(*stmt.variant_parallel_combine_loop.body, body_label, header_label, loop_exit_block, loop_scope, result);
        break;
      }
      case tag_spawn_join: {
        add_block(entry, bbt::spawn_join(stmt.variant_spawn_join.code, exit));
        break;
      }
      case tag_spawn2_join: {
        auto branch1_label = entry;
        auto branch2_label = new_label();
        add_block(branch1_label, bbt::spawn2_join(stmt.variant_spawn2_join.code1, branch2_label));
        add_block(branch2_label, bbt::spawn_join(stmt.variant_spawn2_join.code2, exit));
        break;
      }
      case tag_join_plus: {
        add_block(entry, bbt::join_plus(stmt.variant_join_plus.code, stmt.variant_join_plus.getter, exit));
        break;
      }
      case tag_spawn_minus: {
        add_block(entry, bbt::spawn_minus(stmt.variant_spawn_minus.code, stmt.variant_spawn_minus.getter, exit));
        break;
      }
      case tag_spawn_plus: {
        add_block(entry, bbt::spawn_plus(stmt.variant_spawn_plus.code, stmt.variant_spawn_plus.getter, exit));
        break;
      }
      case tag_join_minus: {
        add_block(entry, bbt::join_minus(stmt.variant_join_minus.getter, exit));
        break;
      }
      case tag_profile_statement: {
#if 0 //defined(HEARTBEAT_ENABLE_LOGGING)
        auto start_profiling_label = entry;
        auto body_label = new_label();
        auto end_profiling_label = new_label();
        auto getter = stmt.variant_profile_statement.getter;
        auto start_profiling = [getter] (sar& s, par& p) {
          auto ws = getter(s, p);
          *(ws.first) = s.pc.profile.work.load();
          *(ws.second) = s.pc.profile.span.load();
        };
        add_block(start_profiling_label, bbt::unconditional_jump(start_profiling, body_label));
        result = transform(*stmt.variant_profile_statement.body, body_label, end_profiling_label, loop_exit_block, loop_scope, result);
        auto reporter = stmt.variant_profile_statement.reporter;
        auto end_profiling = [getter, reporter] (sar& s, par& p) {
          auto ws = getter(s, p);
          auto work = s.pc.profile.work.load() - *(ws.first);
          auto span = s.pc.profile.span.load() - *(ws.second);
          reporter(s, p, work, span);
        };
        add_block(end_profiling_label, bbt::unconditional_jump(end_profiling, exit));
#else
        result = transform(*stmt.variant_profile_statement.body, entry, exit, loop_exit_block, loop_scope, result);
#endif
        break;
      }
      default: {
        break;
      }
    }
    return result;
  }
  
public:
  
  static
  pcfg::cfg_type<Shared_activation_record> transform(stmt_type stmt) {
    lt entry = pcfg::entry_block_label;
    lt exit = pcfg::exit_block_label;
    configuration_type result = transform(stmt, entry, exit, -1, { }, { });
    pcfg::cfg_type<Shared_activation_record> cfg(result.blocks.size());
    for (int i = 0; i < cfg.nb_basic_blocks(); i++) {
      cfg[i] = result.blocks[i];
    }
    auto nb_loops = result.descriptor_of_loop.size();
    assert(nb_loops == next_parallel_loop_id);
    cfg.loop_descriptors.resize(nb_loops);
    for (auto& p : result.descriptor_of_loop) {
      cfg.loop_descriptors[p.first] = p.second;
    }
    cfg.loop_of.resize(cfg.nb_basic_blocks());
    for (int i = 0; i < cfg.nb_basic_blocks(); i++) {
      parallel_loop_id_type loop_id = pcfg::not_a_parallel_loop_id;
      auto it = result.loop_of_basic_block.find(i);
      if (it != result.loop_of_basic_block.end()) {
        loop_id = it->second;
      }
      cfg.loop_of[i] = loop_id;
    }
    return cfg;
  }
  
};
  
template <class Sar, class Par>
typename linearize<Sar, Par>::linearize::lt linearize<Sar, Par>::next_block_label = pcfg::entry_block_label + 1;
  
template <class Sar, class Par>
pcfg::parallel_loop_id_type linearize<Sar, Par>::next_parallel_loop_id = 0;
  
} // end namespace  
} // end namespace
} // end namespace

#endif /*! _HEARTBEAT_DC_H_ */
