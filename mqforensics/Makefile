CC      ?= gcc
CFLAGS  ?= -O3 -std=c11 -Wall -Wextra -Wno-unused-parameter
HTS_CFLAGS := $(shell pkg-config --cflags htslib 2>/dev/null)
HTS_LDLIBS := $(shell pkg-config --libs htslib 2>/dev/null)
LDLIBS  := $(HTS_LDLIBS) -lm

SRC = src/main.c src/analyze.c src/summarize.c src/bam_utils.c src/stats.c src/hist.c
OBJ = $(SRC:.c=.o)
INCLUDES = -Iinclude

all: mq_forensics

mq_forensics: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LDLIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(HTS_CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJ) mq_forensics
