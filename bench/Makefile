####################################
# Example settings.sh 
#
# USE_HWLOC=1
# USE_CILK=1
# CUSTOM_MALLOC_PREFIX=-ltcmalloc 
# CILK_EXTRAS_PREFIX=-L ../../../../cilk-plus-rts/lib -I ../../../../cilk-plus-rts/include -ldl -DCILK_RUNTIME_WITH_STATS
####################################

####################################################################

HWLOC_PREFIX_ALL=-DHAVE_HWLOC $(HWLOC_FLAGS) $(HWLOC_LIBS)

HWLOC_PREFIX=
ifeq ($(strip $(USE_HWLOC)),1)
   HWLOC_PREFIX=$(HWLOC_PREFIX_ALL)
endif

INCLUDE_FILES=$(wildcard $(CHUNKEDSEQ_PATH)/*.hpp) $(wildcard $(HEARTBEAT_INCLUDE_PATH)/*.hpp) $(wildcard $(CMDLINE_PATH)/*.hpp) $(wildcard $(PBBS_INCLUDE_PATH)/*.hpp) $(wildcard $(SPTL_PATH)/*.hpp) $(wildcard $(CACTUS_PATH)/*.hpp) $(wildcard *.hpp)

INCLUDE_PREFIX= -I $(CHUNKEDSEQ_PATH) -I $(HEARTBEAT_INCLUDE_PATH) -I $(CMDLINE_PATH) -I $(PBBS_INCLUDE_PATH) -I $(SPTL_PATH) -I $(CACTUS_PATH)

CILK_PREFIX=
ifeq ($(strip $(USE_CILK)),1)
	CILK_PREFIX=$(CILK_EXTRAS_PREFIX)
	CILK_PREFIX+=-fcilkplus -DUSE_CILK_PLUS -DUSE_CILK_PLUS_RUNTIME -DHEARTBEAT_BENCHMARK
	LD_FLAGS += -lcilkrts -ldl
endif

COMMON_PREFIX= -w -DHAVE_GCC_TLS -fpermissive -pthread -O2 -march=native -DNDEBUG -std=gnu++1y $(HWLOC_PREFIX) -DTARGET_LINUX -m64 -DTARGET_X86_64 -DMANUAL_CONTROL #-DTIME_MEASURE
LOGGING_PREFIX=-DHEARTBEAT_ENABLE_LOGGING
DEBUG_PREFIX= -w -fpermissive -pthread -O0 -std=gnu++11 -g $(LOGGING_PREFIX) -DTARGET_LINUX -DDEBUG_HEARTBEAT_STACK -DHEARTBEAT_RANDOMIZE_SCHEDULE #-DHEARTBEAT_SEQUENCE_USE_PBBS_VERSIONS

%.heartbeat: %.cpp $(INCLUDE_FILES)
	$(CXX) -DHEARTBEAT_ENABLE_STATS $(COMMON_PREFIX) $(INCLUDE_PREFIX) -o $@ $< $(LD_FLAGS) $(CUSTOM_MALLOC_PREFIX)

%.log: %.cpp $(INCLUDE_FILES)
	$(CXX) -DHEARTBEAT_ENABLE_STATS $(COMMON_PREFIX) $(INCLUDE_PREFIX) $(LOGGING_PREFIX) -o $@ $< $(LD_FLAGS) $(CUSTOM_MALLOC_PREFIX)

%.dbg: %.cpp $(INCLUDE_FILES)
	$(CXX) -DHEARTBEAT_ENABLE_STATS $(DEBUG_PREFIX) $(INCLUDE_PREFIX) -o $@ $<

%.cilk: %.cpp $(INCLUDE_FILES) 
	$(CXX) $(CILK_PREFIX) $(COMMON_PREFIX) $(INCLUDE_PREFIX) -o $@ $< $(CUSTOM_MALLOC_PREFIX)

%.cilk_debug: %.cpp $(INCLUDE_FILES) $(PBBS_INCLUDE_PATH)
	$(CXX) $(CILK_PREFIX) $(DEBUG_PREFIX) $(INCLUDE_PREFIX) -o $@ $<

%.cilk_elision: %.cpp $(INCLUDE_FILES) $(PBBS_INCLUDE_PATH)
	$(CXX) $(COMMON_PREFIX) $(INCLUDE_PREFIX) -o $@ $< $(CUSTOM_MALLOC_PREFIX) 

BENCHMARKS=\
	radixsort \
	samplesort \
	suffixarray \
	removeduplicates \
	convexhull \
	removeduplicates \
	nearestneighbors \
	delaunay \
	refine \
	raycast \
	mis \
	mst \
	spanning

dummy:
	echo hi

all: $(BENCHMARKS:%=%.heartbeat) $(BENCHMARKS:%=%.cilk) $(BENCHMARKS:%=%.cilk_elision)

install:
	cp *.heartbeat *.cilk *.cilk_elision $(INSTALL_FOLDER)

clean:
	rm -f *.heartbeat *.cilk *.dbg *.log *.cilk_elision *.cilk_debug
