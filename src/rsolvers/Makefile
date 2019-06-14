CC = gcc
CFLAGS = -Wall -g -pg
LIBS = -lm
OBJS = hllc_sr.o \
  

SRC = $(OBJS:.o=.c)

.c.o:
	${CC} -c $< ${CFLAGS}

.PHONY: clean compile all

all: compile

compile: ${OBJS}

depend: Makedepend
Makedepend:
	${CC} -MM ${CFLAGS} ${SRC} > Makedepend

clean:
	rm -f $(OBJS) Makedepend
include Makedepend
