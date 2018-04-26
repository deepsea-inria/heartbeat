heartbeat
======

A C++ library for prototyping the Heartbeat Scheduling algorithm.

Try it out
----------

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
$ cd heartbeat/example
$ make fib.dbg
$ fib.dbg -algorithm dc -n 20
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Example output:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
exectime 0.367
nb_promotions 5978
nb_steals 0
nb_stacklet_allocations 0
nb_stacklet_deallocations 0
launch_duration 0.367092
utilization 1
exectime 0.367
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Another example
---------------

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
$ fib.dbg -algorithm dc -n 20 -proc 40
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Example output:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
exectime 0.076
nb_promotions 6144
nb_steals 1220
nb_stacklet_allocations 0
nb_stacklet_deallocations 0
launch_duration 0.0771156
utilization 0.8211
exectime 0.078
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now, a real example
-------------------

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
$ cd heartbeat/bench
$ make merge.log
$ merge.log -algorithm heartbeat -n 10000000 -proc 40
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Example output:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
exectime 0.020
nb_promotions 1590
nb_steals 335
nb_stacklet_allocations 0
nb_stacklet_deallocations 0
launch_duration 0.0215029
utilization 0.702188
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Another example
---------------

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
$ convexhull.cilk -algorithm pbbs -proc 40 -type 2d -infile _data/array_point2d_in_circle_large.bin 
$ convexhull.heartbeat -algorithm heartbeat -proc 40 -type 2d -infile _data/array_point2d_in_circle_large.bin 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For viewing the execution
-------------------------

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
git clone https://github.com/deepsea-inria/pview.git
cd pview
make
alias pview=<your path to pview>/pview/pview
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Then:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
$ cd heartbeat/bench 
$ merge.log -algorithm heartbeat -n 10000000 -proc 40 --pview
$ pview
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~