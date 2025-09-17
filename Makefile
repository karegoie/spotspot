# Simple Makefile for spotspot

CC      := gcc
CFLAGS  := -Wall -Wextra -std=c99 -O3
LDFLAGS := -lm

TARGET  := spotspot
SRCS    := matrix.c rng.c stats.c mcmc.c io.c model.c rjmcmc.c spotspot.c
OBJS    := $(SRCS:.c=.o)
DEPS    := $(SRCS:.c=.d)

.PHONY: all clean run

all: $(TARGET)
	@rm -f $(OBJS) $(DEPS)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

-include $(DEPS)

clean:
	rm -f $(OBJS) $(DEPS) $(TARGET)

.PHONY: all clean run build-clean

build-clean: $(TARGET)
	@echo "Built $(TARGET) — removing intermediate files..."
	rm -f $(OBJS) $(DEPS)

# Example run: make run ARGS="testMatrix.csv 100 20 1234 2"
run: $(TARGET)
	./$(TARGET) $(ARGS)
