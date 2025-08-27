# Simple Makefile for mq_forensics
# Usage:
#   make         # build
#   make clean   # remove binary and objects

CC      ?= gcc
CFLAGS  ?= -O3 -std=c11 -Wall -Wextra -Wno-unused-parameter
HTS_CFLAGS := $(shell pkg-config --cflags htslib 2>/dev/null)
HTS_LDLIBS := $(shell pkg-config --libs htslib 2>/dev/null)

# Fallback if pkg-config isn't set up; comment/uncomment as needed:
# HTS_CFLAGS := -I$(HOME)/opt/include
# HTS_LDLIBS := -L$(HOME)/opt/lib -lhts -lpthread -lz -lbz2 -llzma

LDLIBS  := $(HTS_LDLIBS) -lm

.PHONY: all clean

all: mq_forensics

mq_forensics: mq_forensics.c
	$(CC) $(CFLAGS) $(HTS_CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f mq_forensics *.o
