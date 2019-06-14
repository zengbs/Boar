# object files-------------------------------------------------------

CORE_OBJ = convert_var.o \
           cons2prim.o
#           prim2cons.o    

ALL_OBJ = ${CORE_OBJ} ${INTEGRATOR_OBJ} ${RECONSTRUCTION_OBJ} ${RSOLVERS_OBJ}

# macro ddefinitions------------------------------------------------
CC = gcc

CFLAGS = -Wall

DLIB = /usr/local/myprog/gsl-2.4/lib/
DINC = /usr/local/myprog/gsl-2.4/include/

EXEDIR =

BIN = ${EXEDIR}a.out

SRC = $(ALL_OBJ:.o=.c)

# implicit rules-----------------------------------------------------
%.o:%.c
	${CC} -I${DINC} -c $< ${CFLAGS}

# targets------------------------------------------------------------

.PHONY: clean all compile

all: compile

compile: ${BIN}

${BIN}: ${ALL_OBJ}
#	${CC} -o ${EXEDIR}a.out ${ALL_OBJ} -lm ${DLIB}libgslcblas.so.0.0.0 ${DLIB}libgsl.so.23.0.0 ${CFLAGS}   # shared library
	${CC} -o ${EXEDIR}a.out ${ALL_OBJ}  -lm -L${DLIB} -lgsl -lgslcblas  ${CFLAGS}   # shared library
clean:
	rm -f *.o Makedepend $(BIN)

depend: Makedepend

Makedepend:
	${CC} -MM ${CFLAGS} ${SRC} > Makedepend

include Makedepend
