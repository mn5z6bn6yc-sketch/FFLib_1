# Makefile for fflib — finite-field library, demo, and integrity tests
#
# Targets:
#   make        build ./demo
#   make test   build and run ./test_fflib (18 mathematical integrity tests)
#   make clean  remove build artifacts
#
# GMP is required for ff_bigint.c. On Apple Silicon Homebrew installs it
# under /opt/homebrew — the GMP_* lines below handle that. If GMP is in
# your compiler's default paths (Linux, Intel Mac), delete both GMP lines.

CC         ?= gcc
CFLAGS     ?= -Wall -Wextra -O3 -std=c11 -I.
GMP_CFLAGS  = -I/opt/homebrew/include
GMP_LIBS    = -L/opt/homebrew/lib -lgmp

# Library sources (post-refactor layout)
SRCS = ff_lib.c \
       ff_poly.c ff_poly_arith.c ff_poly_factor.c ff_poly_charpoly.c \
       ff_krylov.c ff_frobenius.c ff_mat_util.c \
       ff_snf.c \
       ff_jordan_ext2.c ff_staircase.c ff_ratrec.c \
       ff_extk.c ff_mat_extk.c ff_jordan_extk.c \
       ff_bigint.c \
       poly_random_impl.c

OBJS     = $(SRCS:.c=.o)
HDRS     = $(wildcard *.h) ff_poly_internal.h ff_extk_internal.h
TARGET   = demo
TEST_SRC = tests/test_fflib.c
TEST_OBJ = tests/test_fflib.o
TEST_BIN = test_fflib

all: $(TARGET)

$(TARGET): $(OBJS) demo.o
	$(CC) $(CFLAGS) $(GMP_CFLAGS) -o $@ $^ $(GMP_LIBS)

$(TEST_BIN): $(OBJS) $(TEST_OBJ)
	$(CC) $(CFLAGS) $(GMP_CFLAGS) -o $@ $^ $(GMP_LIBS)

# Rebuild objects when any header changes (cheap at this project size)
%.o: %.c $(HDRS)
	$(CC) $(CFLAGS) $(GMP_CFLAGS) -c $< -o $@

test: $(TEST_BIN)
	./$(TEST_BIN)

clean:
	rm -f $(OBJS) demo.o $(TEST_OBJ) $(TARGET) $(TEST_BIN)

.PHONY: all test clean
