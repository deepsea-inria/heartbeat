TARGET_EXEC ?= fib.dbg

ROOT_DIR=../
EXTERN_DIR=$(ROOT_DIR)/../

BUILD_DBG_DIR ?= ./build_dbg
SRC_DIRS ?= $(ROOT_DIR)/include $(EXTERN_DIR)/chunkedseq/include $(EXTERN_DIR)/cmdline/include $(EXTERN_DIR)/cactus-stack/include

SRCS := $(shell find $(SRC_DIRS) -name '*.cpp' -or -name '*.c' -or -name '*.s')
OBJS_DBG := $(SRCS:%=$(BUILD_DBG_DIR)/%.o)
DEPS_DBG := $(OBJS_DBG:.o=.d)

INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))
INCS := $(shell find $(SRC_DIRS) -name '*.hpp' -or -name '*.h')

CPPFLAGS_SHARED=$(INC_FLAGS) -w -I../bench/ -MMD -MP -pthread -std=gnu++11 -fpermissive -DHEARTBEAT_ENABLE_STATS -DTARGET_LINUX

CPPFLAGS_DBG ?= $(CPPFLAGS_SHARED) -O0 -g -DHEARTBEAT_ENABLE_LOGGING -DDEBUG_HEARTBEAT_STACK -DHEARTBEAT_RANDOMIZE_SCHEDULE

CC_DBG=g++ $(CPPFLAGS_DBG)

ALLDEPS_DBG=$(OBJS_DBG) $(INCS)

CPPFLAGS_OPT ?= $(CPPFLAGS_SHARED) -O2 -march=native -DNDEBUG

CC_OPT=g++ $(CPPFLAGS_OPT)

ALLDEPS_OPT=$(OBJS_OPT) $(INCS)

# assembly
$(BUILD_DBG_DIR)/%.s.o: %.s
	$(MKDIR_P) $(dir $@)
	$(AS) $(ASFLAGS) -c $< -o $@

# c source
$(BUILD_DBG_DIR)/%.c.o: %.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CPPFLAGS_DBG) $(CFLAGS_DBG) -c $< -o $@

# c++ source
$(BUILD_DBG_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CC) $(CPPFLAGS_DBG) -c $< -o $@

%.dbg: %.cpp $(ALLDEPS_DBG)
	$(MKDIR_P) $(dir $@)
	$(CC_DBG) $(OBJS_DBG) $< -o $@ $(LDFLAGS_DBG)

%.opt: %.cpp $(ALLDEPS_OPT)
	$(MKDIR_P) $(dir $@)
	$(CC_OPT) $(OBJS_OPT) $< -o $@ $(LDFLAGS_OPT)

.PHONY: clean

clean:
	$(RM) -r $(BUILD_DBG_DIR) *.dbg

-include $(DEPS_DBG)

MKDIR_P ?= mkdir -p

