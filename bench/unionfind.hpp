// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2010 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef _HEARTBEAT_A_UNIONFIND_INCLUDED
#define _HEARTBEAT_A_UNIONFIND_INCLUDED

namespace heartbeatbench {

struct unionFind {
  int* parents;

  unionFind() { }
  // initialize with all roots marked with -1
  unionFind(int n) {
    parents = (int*) malloc(sizeof(int) * n);
    //    cilk_for(int i = 0; i < n; i++) parents[i] = -1;
  }

  void del() {free(parents);}

  intT find(intT i) {
    if (parents[i] < 0) return i;
    intT j = parents[i];     
    if (parents[j] < 0) return j;
    do j = parents[j]; 
    while (parents[j] >= 0);
    intT tmp;
    while ((tmp = parents[i]) != j) { 
      parents[i] = j;
      i = tmp;
    }
    return j;
  }

  void link(intT u, intT v) { 
    parents[u] = v;
  }
};

} // end namespace

#endif
