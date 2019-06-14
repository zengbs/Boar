CC = gcc
CFLAGS = -Wall -g -pg
LIBS = -lm
OBJS = integrate.o \
       integrate_1d_vl.o \
       integrate_2d_vl.o \
       integrate_3d_vl.o

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
