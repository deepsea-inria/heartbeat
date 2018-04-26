
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "parallel.h"

/* every item in the knapsack has a weight and a value */
#define MAX_ITEMS 256

struct item {
     int value;
     int weight;
};

int best_so_far = INT_MIN;

int compare(struct item *a, struct item *b)
{
     double c = ((double) a->value / a->weight) -
     ((double) b->value / b->weight);

     if (c > 0)
	  return -1;
     if (c < 0)
	  return 1;
     return 0;
}

int read_input(const char *filename, struct item *items, int *capacity, int *n)
{
     int i;
     FILE *f;

     if (filename == NULL)
	  filename = "\0";
     f = fopen(filename, "r");
     if (f == NULL) {
	  fprintf(stderr, "open_input(\"%s\") failed\n", filename);
	  return -1;
     }
     /* format of the input: #items capacity value1 weight1 ... */
     fscanf(f, "%d", n);
     fscanf(f, "%d", capacity);

     for (i = 0; i < *n; ++i)
	  fscanf(f, "%d %d", &items[i].value, &items[i].weight);

     fclose(f);

     /* sort the items on decreasing order of value/weight */
     /* cilk2c is fascist in dealing with pointers, whence the ugly cast */
     qsort(items, *n, sizeof(struct item),
	    (int (*)(const void *, const void *)) compare);

     return 0;
}

namespace pbbs {
  
/* 
 * return the optimal solution for n items (first is e) and
 * capacity c. Value so far is v.
 */
int knapsack(struct item *e, int c, int n, int v)
{
  int with, without, best;
  double ub;
  
  /* base case: full knapsack or no items */
  if (c < 0)
	  return INT_MIN;
  
  if (n == 0 || c == 0)
	  return v;		/* feasible solution, with value v */
  
  ub = (double) v + c * e->value / e->weight;
  
  if (ub < best_so_far) {
	  /* prune ! */
	  return INT_MIN;
  }
  /* 
   * compute the best solution without the current item in the knapsack 
   */
  without = cilk_spawn knapsack(e + 1, c, n - 1, v);
  
  /* compute the best solution with the current item in the knapsack */
  with = knapsack(e + 1, c - e->weight, n - 1, v + e->value);
  
  cilk_sync;
  
  best = with > without ? with : without;
  
  /* 
   * notice the race condition here. The program is still
   * correct, in the sense that the best solution so far
   * is at least best_so_far. Moreover best_so_far gets updated
   * when returning, so eventually it should get the right
   * value. The program is highly non-deterministic.
   */
  if (best > best_so_far)
	  best_so_far = best;
  
  return best;
}

} // end namespace

#include "heartbeatbench.hpp"

namespace heartbeatbench {

class knapsack : public heartbeat::edsl::pcfg::shared_activation_record {
public:

  struct item *e; int c; int n; int v; int* dest;
  int with, without, best;
  double ub;

  knapsack(struct item *e, int c, int n, int v, int* dest)
    : e(e), c(c), n(n), v(v), dest(dest) { }

  heartbeat_dc_declare(heartbeat::edsl, knapsack, sar, par, dc, get_dc)

  static
  dc get_dc() {
    return dc::stmts({
      dc::mk_if([] (sar& s, par& p) { return s.c < 0; }, dc::stmts({
        /* base case: full knapsack or no items */
        dc::stmt([] (sar& s, par& p) {
          *s.dest = INT_MIN;
        }),
        dc::exit_function()
      })),
      dc::mk_if([] (sar& s, par& p) { return (s.n == 0) || (s.c == 0); }, dc::stmts({
        /* feasible solution, with value v */
        dc::stmt([] (sar& s, par& p) {
          *s.dest = s.v;
        }),
        dc::exit_function()
      })),
      dc::stmt([] (sar& s, par& p) {
        s.ub = (double) s.v + s.c * s.e->value / s.e->weight;
      }),
      dc::mk_if([] (sar& s, par& p) { return s.ub < best_so_far; }, dc::stmts({
        /* prune ! */
        dc::stmt([] (sar& s, par& p) {
          *s.dest = INT_MIN;
        }),
        dc::exit_function()
      })),
      dc::spawn2_join(
        [] (sar& s, par&, plt pt, stt st) {
          return heartbeat_call<knapsack>(st, pt, s.e + 1, s.c, s.n - 1, s.v, &s.without);
        },
        [] (sar& s, par&, plt pt, stt st) {
          return heartbeat_call<knapsack>(st, pt, s.e + 1, s.c - s.e->weight, s.n - 1, s.v + s.e->value, &s.with);
      }),
      dc::stmt([] (sar& s, par& p) {
        s.best = s.with > s.without ? s.with : s.without;

        /* 
         * notice the race condition here. The program is still
         * correct, in the sense that the best solution so far
         * is at least best_so_far. Moreover best_so_far gets updated
         * when returning, so eventually it should get the right
         * value. The program is highly non-deterministic.
         */
        if (s.best > best_so_far)
          best_so_far = s.best;

        *s.dest = s.best;
      })
    });
  }
  
};

heartbeat_pcfg_allocate(knapsack, get_cfg)
  
} // end namespace

namespace cmdline = deepsea::cmdline;

void benchmark(std::string infile) {
  auto filename = infile.c_str();
  struct item items[MAX_ITEMS];	/* array of items */
  int n, capacity, sol;
  if (read_input(filename, items, &capacity, &n)) {
	  return;
  }
  deepsea::cmdline::dispatcher d;
  d.add("heartbeat", [&] {
    heartbeat::launch_interpreter<heartbeatbench::knapsack>(items, capacity, n, 0, &sol);
  });
  d.add("pbbs", [&] {
    heartbeatbench::run_and_report_elapsed_time([&] {
      sol = pbbs::knapsack(items, capacity, n, 0);
    });
  });
  d.dispatch("algorithm");
  std::cout << "solution " << sol << std::endl; 
}

int main(int argc, char** argv) {
  heartbeatbench::initialize(argc, argv);
  std::string infile = deepsea::cmdline::parse_or_default_string("infile", "");
  if (infile == "") {
    std::cout << "bogus input file name" << std::endl;
    exit(0);
  }
  benchmark(infile);  
  return 0;
}
