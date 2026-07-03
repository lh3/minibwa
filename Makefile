CC ?= gcc
AR ?= ar

BUILD_DIR ?= build
OBJ_DIR := $(BUILD_DIR)/obj

PROG := minibwa
LIB_TARGET := libminibwa.a

CFLAGS ?= -std=c99 -g -Wall -O3
CPPFLAGS ?=
LDFLAGS ?=
LDLIBS ?= -lpthread -lz -lm

INCLUDES := -Iinclude -Isrc -Ithird_party/libsais
ARCH := $(shell uname -m)
omp ?= $(shell printf '\043include <omp.h>\nint main(){return 0;}' | $(CC) -x c -fopenmp -o /dev/null - 2>/dev/null && echo "1" || echo "0")

LIB_SRCS := \
	src/kommon.c \
	src/kalloc.c \
	src/bwt.c \
	src/l2bit.c \
	src/options.c \
	src/seed.c \
	src/map-algo.c \
	src/lchain.c \
	src/align.c \
	src/pe.c \
	src/cs.c \
	src/format.c \
	src/ksw2_extz2_sse.c \
	src/ksw2_extd2_sse.c \
	src/ksw2_ll_sse.c

APP_SRCS := \
	src/kthread.c \
	third_party/libsais/libsais.c \
	third_party/libsais/libsais64.c \
	src/index.c \
	src/bseq.c \
	src/map-main.c \
	src/fastmap.c

GPL_SRCS := \
	third_party/bwtgen/QSufSort.c \
	third_party/bwtgen/bwtgen.c

MAIN_SRC := src/main.c
MIMALLOC_SRC := third_party/mimalloc/static.c

LIB_OBJS := $(patsubst %.c,$(OBJ_DIR)/%.o,$(LIB_SRCS))
APP_OBJS := $(patsubst %.c,$(OBJ_DIR)/%.o,$(APP_SRCS))
MAIN_OBJ := $(patsubst %.c,$(OBJ_DIR)/%.o,$(MAIN_SRC))
GPL_OBJS := $(patsubst %.c,$(OBJ_DIR)/%.o,$(GPL_SRCS))
MIMALLOC_OBJ := $(patsubst %.c,$(OBJ_DIR)/%.o,$(MIMALLOC_SRC))

ifneq ($(asan),)
	CFLAGS += -fsanitize=address
	LDFLAGS += -fsanitize=address
	LDLIBS += -ldl
endif

ifeq ($(omp),1)
	CPPFLAGS += -DLIBSAIS_OPENMP
	CFLAGS += -fopenmp
	LDLIBS += -fopenmp
endif

ifneq ($(gpl),0)
	APP_OBJS += $(GPL_OBJS)
	CPPFLAGS += -DUSE_GPL
endif

ifeq ($(mimalloc),0)
	MIMALLOC_OBJ :=
	CPPFLAGS += -DHAVE_KALLOC
endif

ifeq ($(ARCH),x86_64)
	CFLAGS += -msse4.2 -mpopcnt
endif

DEPS := $(LIB_OBJS:.o=.d) $(APP_OBJS:.o=.d) $(MAIN_OBJ:.o=.d) $(MIMALLOC_OBJ:.o=.d)

.PHONY: all clean examples

all: $(PROG)

$(LIB_TARGET): $(LIB_OBJS)
	$(AR) rcs $@ $^

$(PROG): $(LIB_TARGET) $(MIMALLOC_OBJ) $(APP_OBJS) $(MAIN_OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) $(MIMALLOC_OBJ) $(APP_OBJS) $(MAIN_OBJ) $(LIB_TARGET) -o $@ $(LDLIBS)

examples: $(LIB_TARGET)
	$(MAKE) -C examples

$(OBJ_DIR)/third_party/mimalloc/static.o: third_party/mimalloc/static.c
	@mkdir -p $(dir $@)
	$(CC) -c -std=gnu11 -O3 -Wall -Wextra -DNDEBUG -DMI_MALLOC_OVERRIDE -DMI_OSX_INTERPOSE=1 -DMI_OSX_ZONE=1 -MMD -MP -Ithird_party/mimalloc $< -o $@

$(OBJ_DIR)/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) -MMD -MP $< -o $@

clean:
	rm -rf $(BUILD_DIR) $(PROG) $(LIB_TARGET) a.out *.o *.a *.dSYM *~
	$(MAKE) -C examples clean

-include $(DEPS)
