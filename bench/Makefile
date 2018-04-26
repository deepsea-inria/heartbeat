
PACKAGE_PATH=../../

PBENCH_PATH=$(PACKAGE_PATH)/pbench
CMDLINE_PATH=$(PACKAGE_PATH)/cmdline/include
CHUNKEDSEQ_PATH=$(PACKAGE_PATH)/chunkedseq/include
SPTL_PATH=$(PACKAGE_PATH)/sptl/include
PBBS_INCLUDE_PATH=$(PACKAGE_PATH)/pbbs-include/
PBBS_SPTL_PATH=$(PACKAGE_PATH)/pbbs-sptl/include
PBBS_SPTL_BENCH_PATH=$(PACKAGE_PATH)/pbbs-sptl/bench
CACTUS_PATH=$(PACKAGE_PATH)/cactus-stack/include
HEARTBEAT_INCLUDE_PATH=../include/

CUSTOM_MALLOC_PREFIX=

HWLOC_FLAGS=`pkg-config --cflags hwloc`
HWLOC_LIBS=`pkg-config --libs hwloc`
HWLOC_PREFIX_ALL=-DHAVE_HWLOC $(HWLOC_FLAGS) $(HWLOC_LI

####################################################################
# Makefile options

# Create a file called "settings.sh" in this folder if you want to
# configure particular options. See section below for options.

-include settings.sh

####################################
# Example settings.sh 
#
# USE_HWLOC=1
# USE_CILK=1
# CUSTOM_MALLOC_PREFIX=-ltcmalloc 
# CILK_EXTRAS_PREFIX=-L ../../../../cilk-plus-rts/lib -I ../../../../cilk-plus-rts/include -ldl -DCILK_RUNTIME_WITH_STATS
####################################

####################################################################

include $(PBENCH_PATH)/Makefile_common

HWLOC_FLAGS=`pkg-config --cflags hwloc`
HWLOC_LIBS=`pkg-config --libs hwloc`
ifeq ($(strip $(USE_MANUAL_HWLOC_PATH)),1)
	HWLOC_FLAGS=$(MY_HWLOC_FLAGS)
	HWLOC_LIBS=$(MY_HWLOC_LIBS)
endif

HWLOC_PREFIX_ALL=-DHAVE_HWLOC $(HWLOC_FLAGS) $(HWLOC_LIBS)

HWLOC_PREFIX=
ifeq ($(strip $(USE_HWLOC)),1)
   HWLOC_PREFIX=$(HWLOC_PREFIX_ALL)
endif

INCLUDE_FILES=$(wildcard $(CHUNKEDSEQ_PATH)/*.hpp) $(wildcard $(HEARTBEAT_INCLUDE_PATH)/*.hpp) $(wildcard $(CMDLINE_PATH)/*.hpp) $(wildcard $(PBBS_INCLUDE_PATH)/*.hpp) $(wildcard $(SPTL_PATH)/*.hpp) $(wildcard $(PBBS_SPTL_PATH)/*.hpp) $(wildcard $(PBBS_SPTL_BENCH_PATH)/*.hpp) $(wildcard $(CACTUS_PATH)/*.hpp) $(wildcard *.hpp)

INCLUDE_PREFIX= -I $(CHUNKEDSEQ_PATH) -I $(HEARTBEAT_INCLUDE_PATH) -I $(CMDLINE_PATH) -I $(PBBS_INCLUDE_PATH) -I $(SPTL_PATH) -I $(PBBS_SPTL_PATH) -I $(PBBS_SPTL_BENCH_PATH) -I $(CACTUS_PATH)

CILK_PREFIX=
ifeq ($(strip $(USE_CILK)),1)
	CILK_PREFIX= -fcilkplus -lcilkrts -ldl -DUSE_CILK_PLUS -DUSE_CILK_PLUS_RUNTIME -DHEARTBEAT_BENCHMARK
	CILK_PREFIX+=$(CILK_EXTRAS_PREFIX)
endif

COMMON_PREFIX= -w -DHAVE_GCC_TLS -fpermissive -pthread -O2 -march=native -DNDEBUG -std=gnu++1y $(HWLOC_PREFIX) -DTARGET_LINUX -m64 -DTARGET_X86_64 -DMANUAL_CONTROL #-DTIME_MEASURE
LOGGING_PREFIX=-DHEARTBEAT_ENABLE_LOGGING
DEBUG_PREFIX= -w -fpermissive -pthread -O0 -std=gnu++11 -g $(LOGGING_PREFIX) -DTARGET_LINUX -DDEBUG_HEARTBEAT_STACK -DHEARTBEAT_RANDOMIZE_SCHEDULE #-DHEARTBEAT_SEQUENCE_USE_PBBS_VERSIONS

%.heartbeat: %.cpp $(INCLUDE_FILES)
	g++ -DHEARTBEAT_ENABLE_STATS $(COMMON_PREFIX) $(CUSTOM_MALLOC_PREFIX) $(INCLUDE_PREFIX) -o $@ $<

%.log: %.cpp $(INCLUDE_FILES)
	g++ -DHEARTBEAT_ENABLE_STATS $(COMMON_PREFIX) $(CUSTOM_MALLOC_PREFIX) $(INCLUDE_PREFIX) $(LOGGING_PREFIX) -o $@ $<

%.dbg: %.cpp $(INCLUDE_FILES)
	g++ -DHEARTBEAT_ENABLE_STATS $(DEBUG_PREFIX) $(INCLUDE_PREFIX) -o $@ $<

%.cilk: %.cpp $(INCLUDE_FILES) $(PBBS_INCLUDE_PATH)
	g++ $(CILK_PREFIX) $(COMMON_PREFIX) $(CUSTOM_MALLOC_PREFIX) $(INCLUDE_PREFIX) -o $@ $<

%.cilk_debug: %.cpp $(INCLUDE_FILES) $(PBBS_INCLUDE_PATH)
	g++ $(CILK_PREFIX) $(DEBUG_PREFIX) $(CUSTOM_MALLOC_PREFIX) $(INCLUDE_PREFIX) -o $@ $<

%.cilk_elision: %.cpp $(INCLUDE_FILES) $(PBBS_INCLUDE_PATH)
	g++ $(COMMON_PREFIX) $(CUSTOM_MALLOC_PREFIX) $(INCLUDE_PREFIX) -o $@ $<

clean: pbench_clean
	rm -f *.heartbeat *.cilk *.dbg *.log *.cilk_elision *.cilk_debug
