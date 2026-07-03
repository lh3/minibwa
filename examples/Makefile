CC ?= gcc
ROOT := ..

CFLAGS ?= -std=gnu99 -Wall -O3
CPPFLAGS ?= -I$(ROOT)/include -I$(ROOT)/src
LDLIBS ?= -lz -lm

EXE := mbmap-one mbmap-batch
LIBMINIBWA := $(ROOT)/libminibwa.a

.PHONY: all clean

all: $(EXE)

$(LIBMINIBWA):
	$(MAKE) -C $(ROOT) libminibwa.a

mbmap-one: ex-one.c $(LIBMINIBWA)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)

mbmap-batch: ex-batch.c $(LIBMINIBWA)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -f $(EXE)
