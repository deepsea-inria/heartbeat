#include <typeinfo>

#include "pcfg.hpp"
#include "interpreter.hpp"
#include "dc.hpp"

#ifndef _HEARTBEAT_EDSL_H_
#define _HEARTBEAT_EDSL_H_

/*---------------------------------------------------------------------*/
/* Macros */

#define heartbeat_pcfg_default_private_activation_record(edsl) \
class private_activation_record : public edsl::private_activation_record { };

#define heartbeat_pcfg_declare(edsl, name, sar, par, bb) \
std::pair<edsl::pcfg::stack_type, heartbeat::fuel::check_type> run(edsl::pcfg::stack_type stack) const { \
  return edsl::pcfg::step(cfg, stack); \
} \
\
void promote_mark(edsl::pcfg::interpreter* interp, edsl::pcfg::private_activation_record* p) { \
  edsl::pcfg::promote_mark(cfg, interp, this, (par*)p); \
} \
\
heartbeat::sched::future get_dependency_of_join_minus(edsl::pcfg::stack_type  stack) { \
  return edsl::pcfg::get_dependency_of_join_minus(cfg, stack); \
} \
\
using cfg_type = edsl::pcfg::cfg_type<name>; \
using sar = name; \
using par = private_activation_record; \
using stt = edsl::pcfg::stack_type; \
using bb = typename cfg_type::basic_block_type; \
using plt = edsl::pcfg::cactus::parent_link_type; \
\
template <class Activation_record, class ...Args> \
static stt heartbeat_call(stt s, plt p, Args... args) { \
  return edsl::pcfg::push_call<Activation_record>(s, p, args...); \
} \
\
static \
cfg_type cfg; \

#define heartbeat_get_name(__name) \
const char* get_name() { \
  return typeid(__name).name(); \
} \

#define heartbeat_dc_declare(edsl, __name, sar, par, dc, get_dc) \
heartbeat_pcfg_default_private_activation_record(edsl::pcfg) \
heartbeat_pcfg_declare(edsl, __name, sar, par, _bb0123) \
using dc = edsl::dc::stmt_type<sar, par>; \
\
static \
cfg_type get_cfg() { \
  return edsl::dc::linearize<sar, par>::transform(get_dc()); \
} \
\
heartbeat_get_name(__name) \

#define heartbeat_pcfg_allocate(name, get_cfg) \
name::cfg_type name::cfg = name::get_cfg(); \

#define heartbeat_dc_loop_declare(edsl, __name, sar, par, dc, get_dc) \
heartbeat_pcfg_declare(edsl, __name, sar, par, _bb0123) \
using dc = edsl::dc::stmt_type<sar, par>; \
\
static \
cfg_type get_cfg() { \
return edsl::dc::linearize<sar, par>::transform(get_dc()); \
} \
\
heartbeat_get_name(__name) \

#define heartbeat_pcfg_allocate(name, get_cfg) \
name::cfg_type name::cfg = name::get_cfg(); \

#define heartbeat_parallel_loop_alloc_default(edsl, nb_loops) \
private_activation_record() { \
  private_activation_record::initialize_descriptors(); \
} \
\
std::array<edsl::pcfg::parallel_for_activation_record, nb_loops> _ars1; \
edsl::pcfg::parallel_loop_activation_record* _heartbeat_loop_activation_record_of(edsl::pcfg::parallel_loop_id_type id) { \
  return &_ars1[id]; \
}

#define heartbeat_private_activation_record_begin(edsl, name, nb_loops) \
class private_activation_record \
: public edsl::pcfg::parallel_loop_private_activation_record<name,private_activation_record> { \
public: \
  heartbeat_parallel_loop_alloc_default(edsl, nb_loops)

#define heartbeat_private_activation_record_end(edsl, name, sar, par, dc, get_dc) \
}; \
heartbeat_dc_loop_declare(edsl, name, sar, par, dc, get_dc)

#endif /*! _HEARTBEAT_EDSL_H_ */
