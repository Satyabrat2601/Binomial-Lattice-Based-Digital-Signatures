# DUET-Sign Reference Implementation -- Makefile
# Usage: make MODE=1   (or 3, or 5)
#        make          (builds all three modes)

CC     ?= gcc
CFLAGS  = -std=c11 -Wall -Wextra -O3 -fomit-frame-pointer
MODE   ?= 1

ifeq ($(MODE),1)
  MDEF := DUET_MODE1
  MSTR := mode1
else ifeq ($(MODE),3)
  MDEF := DUET_MODE3
  MSTR := mode3
else ifeq ($(MODE),5)
  MDEF := DUET_MODE5
  MSTR := mode5
else
  $(error MODE must be 1, 3, or 5)
endif

CFLAGS += -DDUET_CONFIG_MODE=$(MDEF) -Iinclude

BUILDDIR := build/$(MSTR)

SRCS_CORE := \
    src/fips202.c \
    src/reduce.c \
    src/ntt.c \
    src/poly.c \
    src/polyvec.c \
    src/polymat.c \
    src/packing.c \
    src/symmetric-shake.c \
    src/randombytes.c \
    src/encoding.c \
    src/sign.c

OBJS_CORE := $(patsubst src/%.c, $(BUILDDIR)/%.o, $(SRCS_CORE))

TARGET_TEST  := $(BUILDDIR)/test_$(MSTR)
TARGET_BENCH := $(BUILDDIR)/bench_$(MSTR)

.PHONY: all test bench clean

all:
	$(MAKE) MODE=1 test
	$(MAKE) MODE=3 test
	$(MAKE) MODE=5 test

$(BUILDDIR):
	@mkdir -p $@

$(BUILDDIR)/%.o: src/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(TARGET_TEST): $(OBJS_CORE) test/main.c | $(BUILDDIR)
	$(CC) $(CFLAGS) $(OBJS_CORE) test/main.c -o $@

$(TARGET_BENCH): $(OBJS_CORE) benchmark/speed.c benchmark/cpucycles.c \
                  benchmark/speed_print.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -Ibenchmark \
	      $(OBJS_CORE) benchmark/speed.c benchmark/cpucycles.c \
	      benchmark/speed_print.c -o $@

test: $(TARGET_TEST)
	@echo ""
	@echo "=== DUET-Sign $(MSTR) tests ==="
	./$(TARGET_TEST)

bench: $(TARGET_BENCH)
	@echo ""
	@echo "=== DUET-Sign $(MSTR) benchmark ==="
	./$(TARGET_BENCH)

clean:
	rm -rf build/
