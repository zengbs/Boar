# object files-------------------------------------------------------

CORE_OBJ = array.o \
           bvals_hydro.o \
           init_grid.o \
           problem.o \
           convert_var.o \
           cc_pos.o \
           output.o \
           main.o

INTEGRATOR_OBJ = integrator/integrate.o \
		 integrator/integrate_1d_vl.o \
		 integrator/integrate_2d_vl.o \
		 integrator/integrate_3d_vl.o

RECONSTRUCTION_OBJ = reconstruction/lr_states_prim2.o

RSOLVERS_OBJ = rsolvers/hllc_sr.o

ALL_OBJ = ${CORE_OBJ} ${INTEGRATOR_OBJ} ${RECONSTRUCTION_OBJ} ${RSOLVERS_OBJ}

# macro ddefinitions------------------------------------------------
CC = gcc

CFLAGS = -Wall -g -DSR_DEBUG

DLIB = /usr/local/myprog/gsl-2.4/lib/
DINC = /usr/local/myprog/gsl-2.4/include/

EXEDIR = ../bin/

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
